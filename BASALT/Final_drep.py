#!/usr/bin/env python

"""
Final dereplication and reporting utilities for BASALT.

This module performs PCA-based dereplication on bin sets and produces
summary statistics for downstream evaluation.
"""

from Bio import SeqIO
import sys, os, threading, copy
from multiprocessing import Pool
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from basalt_runtime import depth_sample_count

def Contigs_aligner(Contigs_assembly1, num_threads):
    print('Using BLAST to align '+Contigs_assembly1+' to '+Contigs_assembly1)
    try:
        blast_output=str(Contigs_assembly1)+'_vs_'+str(Contigs_assembly1)+'.txt'
        os.system('makeblastdb -in '+str(Contigs_assembly1)+' -dbtype nucl -hash_index -parse_seqids -logfile temp_db.txt')
        os.system('blastn -query '+str(Contigs_assembly1)+' -db '+str(Contigs_assembly1)+' -evalue 1e-20 -outfmt 6 -num_threads '+str(num_threads)+' -out '+str(blast_output))
    except:
        print('Alignment error! Please check whether BLAST+ is installed in your system')
    
    blast_output2=open('Filtrated_'+str(Contigs_assembly1)+'_vs_'+str(Contigs_assembly1)+'.txt','w')
    record_result = {}
    for line in open(blast_output,'r'):
        simi=str(line).strip().split('\t')[2]
        if float(simi) >= 99:
            length=eval(str(line).strip().split('\t')[3])
            contig1=str(line).strip().split('\t')[0]
            contig2=str(line).strip().split('\t')[1]

            if contig1 != contig2:
                if int(length) >= 500:
                    blast_output2.write(str(line))
            else:
                if '||' in contig1:
                    if str(line) not in record_result.keys():
                        blast_output2.write(str(line))
                        record_result[str(line)]=0
    blast_output2.close()
    return 'Filtrated_'+str(Contigs_assembly1)+'_vs_'+str(Contigs_assembly1)+'.txt'

def Sequence_length_recorder(ORF):
    print('Reading '+ORF+' Contigs length')
    seq_len={}
    for record in SeqIO.parse(ORF,'fasta'):
        seq_len[str(record.id)]=len(record.seq)
    print('---------------------')
    return seq_len

def PCA_slector(data_array, num_contig):
    pca = PCA(n_components=1)
    pca.fit(data_array)
    explained_variance_ratio=pca.explained_variance_ratio_
    # print(explained_variance_ratio)
    # num=len(explained_variance_ratio)
    newData=pca.fit_transform(data_array)
    newData2=newData.reshape((1,num_contig))
    # print('Shape', num, num_contig
    newData_list=newData2.tolist()
    n=0
    for item in newData_list:
        n+=1
        if n == 1:
            newData_list_item=item
    return newData_list_item, explained_variance_ratio

def splitting_kmer_file(bin_name, assembly, bin_contigs_mock, pwd):
    """Write the assembly k-mer rows belonging to one candidate bin."""
    source_path = os.path.join(pwd, str(assembly) + '.kmer.txt')
    output_path = os.path.join(os.getcwd(), str(bin_name) + '_kmer.txt')
    retained = bin_contigs_mock[bin_name]
    with open(source_path, 'r') as source, open(output_path, 'w') as output:
        for line in source:
            contig_id = line.split('\t', 1)[0].strip()
            if contig_id in retained:
                output.write(line)


def outlier_remover(
    bin_id,
    contig_ids,
    thresholds,
    item_data,
    explained_variance_ratio,
    bin_outliers,
    feature_type,
):
    """Record per-threshold IQR outliers for the legacy TNF refinement helper."""
    del explained_variance_ratio, feature_type
    summary = pd.Series(item_data).describe()
    q1 = summary['25%']
    q3 = summary['75%']
    iqr = q3 - q1
    for threshold in thresholds:
        lower = q1 - float(threshold) * iqr
        upper = q3 + float(threshold) * iqr
        bin_outliers[bin_id][threshold] = [
            contig_ids[index]
            for index, value in enumerate(item_data)
            if value < lower or value > upper
        ]
    return bin_outliers


def core_contigs_filtration(bin_contig_cov, bin_contig, contig_cov, folder_binset, contig_file):
    pwd=os.getcwd()
    print('Filtrating core contigs')

    core_contigs, core_contigs_IQR, core_contigs_IQR_coverage, core_contigs_IQR_ave_coverage ={}, {}, {}, {}
    for bin in bin_contig_cov.keys():
        print('Processing '+str(bin))
        core_contigs[bin]={}
        core_contigs_IQR[bin]={}
        core_contigs_IQR_coverage[bin]={}
        core_contigs_IQR_ave_coverage[bin]={}
        # contig_num=len(bin_contig_cov[bin])

        ### Filtration of contigs
        print(str(len(bin_contig_cov[bin]))+' contigs')
        cov_index_total_num, cov_index={}, {}
        for item in bin_contig_cov[bin].keys():  
            num_coverage=len(bin_contig_cov[bin][item])
            for i in range(0, num_coverage):
                try:
                    cov_index_total_num[i]+=float(bin_contig_cov[bin][item][i+1])
                    cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))
                except:
                    cov_index_total_num[i]=float(bin_contig_cov[bin][item][i+1])
                    cov_index[i]=[]
                    cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))

        upper, lower, upper_CC, lower_CC = {}, {}, {}, {}
        for item in cov_index.keys():
            cov_index[item].sort() ### sort low to high
            list_key_num=len(cov_index[item])
            p25=int(0.25*list_key_num)
            p75=int(0.75*list_key_num)

            n=0
            for i in cov_index[item]:
                n+=1
                if n == p25:
                    Q1=i
                elif n == p75:
                    Q3=i
                else:
                    continue

            IQR = Q3 - Q1
            upper[item+1] = Q3 + 1.5*IQR
            lower[item+1] = Q1 - 1.5*IQR
            upper_CC[item+1] = Q3
            lower_CC[item+1] = Q1
            print('Coverage:'+str(item)+'; 75%: '+str(Q3)+'; 25%: '+str(Q1))

        for i in range(1, num_coverage+1):
            core_contigs_IQR_coverage[bin][i]=0
            core_contigs_IQR_ave_coverage[bin][i]=0

        for i in range(1, num_coverage+1):
            for item in bin_contig_cov[bin].keys():
                if contig_cov[item][i] <= upper[i] and contig_cov[item][i] >= lower[i]:
                    try:
                        core_contigs[bin][item]+=1
                    except:
                        core_contigs[bin][item]=1
                    
                    if contig_cov[item][i] <= upper_CC[i] and contig_cov[item][i] >= lower_CC[i]:
                        # core_contigs_IQR_coverage[bin][i]+=contig_cov[item][i]
                        try:
                            core_contigs_IQR[bin][item]+=1
                        except:
                            core_contigs_IQR[bin][item]=1   
        print(str(len(core_contigs_IQR[bin]))+' contigs between 75% and 25%')             
        
    core_contigs2={}
    core_contigs2=copy.deepcopy(core_contigs)
    core_contigs_IQR2={}
    core_contigs_IQR2=copy.deepcopy(core_contigs_IQR)  

    ftest=open('core_contigs_test_'+str(folder_binset)+'.txt','w')
    for bin_id in core_contigs2.keys():
        for contigs in core_contigs2[bin_id].keys():
            ftest.write(str(bin_id)+'\t'+str(contigs)+'\t'+str(core_contigs2[bin_id][contigs])+'\n')
            if core_contigs2[bin_id][contigs] != num_coverage:
                del core_contigs[bin_id][contigs]
    ftest.close()

    for bin_id in core_contigs_IQR2.keys():
        for contigs in core_contigs_IQR2[bin_id].keys():
            if core_contigs_IQR2[bin_id][contigs] != num_coverage:
                del core_contigs_IQR[bin_id][contigs]

    core_contigs2, core_contigs_IQR2={}, {}
    for bin_id in bin_contig_cov.keys():
        print('Processing '+str(bin_id)+' coverage matrix: '+str(len(bin_contig_cov[bin_id]))+' contigs.')
        core_contigs2[bin_id]={}
        core_contigs_IQR2[bin_id]={}
        coverage_data, contigs_ids, coverage_list, num_contig={}, [], [], 0
        for contig in bin_contig_cov[bin_id].keys():
            num_contig+=1
            contigs_ids.append(contig)
            num_coverage=len(bin_contig_cov[bin_id][contig])
            for i in range(1, num_coverage+1):
                if i not in coverage_data.keys():
                    coverage_data[i]=[]
                coverage_data[i].append(bin_contig_cov[bin_id][contig][i])
                coverage_list.append(bin_contig_cov[bin_id][contig][i])

        try:
            coverage_array=np.array(coverage_list).reshape((num_contig,num_coverage))

            A=PCA_slector(coverage_array, num_contig)
            newData=A[0]
            # explained_variance_ratio=A[1]
            four = pd.Series(newData).describe()
            Q1 = four['25%']
            Q3 = four['75%']
            IQR = Q3 - Q1
            upper1 = Q3 + 1.5*IQR
            lower1 = Q1 - 1.5*IQR
            upper_CC = Q3
            lower_CC = Q1
            for i in range(0, len(newData)):
                if newData[i] <= float(upper1) and newData[i] >= float(lower1):
                    passed_contig=str(contigs_ids[i])
                    # core_contigs2[bin_id][passed_contig]=0
                    if passed_contig in core_contigs[bin_id].keys():
                        core_contigs2[bin_id][passed_contig]=0
                    if newData[i] <= float(upper_CC) and newData[i] >= float(lower_CC):
                        core_contigs_IQR2[bin_id][passed_contig]=0
        except:
            print(str(bin_id)+' PCA error. Adding all contigs as core contigs')
            for contig in bin_contig_cov[bin_id].keys():
                core_contigs_IQR2[bin_id][contig]=0

    for bin_id in core_contigs_IQR2.keys():
        for contigs in core_contigs_IQR2[bin_id].keys():
            for i in range(1, num_coverage+1):
                core_contigs_IQR_coverage[bin_id][i]+=contig_cov[contigs][i]
    
    for bin_id in core_contigs_IQR_coverage.keys():
        for i in core_contigs_IQR_coverage[bin_id].keys():
            if len(core_contigs_IQR2[bin_id]) != 0:
                core_contigs_IQR_ave_coverage[bin_id][i]=core_contigs_IQR_coverage[bin_id][i]/len(core_contigs_IQR2[bin_id])
            else:
                core_contigs_IQR_ave_coverage[bin_id][i]=0
    
    f=open(str(folder_binset)+'_bins_coverage.txt','w')
    i=1
    a='Bin'
    while i <= int(num_coverage):
        a+='\t'+'Coverage'+str(i)
        i+=1
    f.write(str(a)+'\n')

    for bin_id in core_contigs_IQR_ave_coverage.keys():
        a=str(bin_id)
        for i in core_contigs_IQR_ave_coverage[bin_id].keys():
            a+='\t'+str(core_contigs_IQR_ave_coverage[bin_id][i])
        f.write(str(a)+'\n')
    f.close()

    print('Coverage filtrating core contigs done.')

    os.system('mkdir CC_'+folder_binset)
    os.chdir('CC_'+folder_binset)
    new_contigs={}
    for bin_id in core_contigs2.keys():
        f=open(bin_id, 'w')
        for contigs in core_contigs2[bin_id].keys():
            f.write('>'+str(contigs)+'\n'+str(bin_contig[bin_id][contigs])+'\n')
            new_contigs[contigs]=str(bin_contig[bin_id][contigs])
        f.close()
    os.chdir(pwd)

    fnc=open('CC_'+contig_file,'w')
    for contigs in new_contigs.keys():
        fnc.write('>'+str(contigs)+'\n'+str(new_contigs[contigs])+'\n')
    fnc.close()

    os.chdir(folder_binset)
    for root, dirs, files in os.walk(pwd+'/'+folder_binset):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                os.system('cp '+file+' '+pwd+'/CC_'+folder_binset)
    os.chdir(pwd)
    return core_contigs2, core_contigs_IQR_ave_coverage

def bin_depth_normalization(target_bin, bin_dict, binset_coverage_total, contigs_coverage, pwd, num_threads):
    # try:
    #     f=open(target_bin+'_bins_depth_comparison.txt','a')
    #     f.close()
    # except:
    f=open(target_bin+'_bins_depth_comparison.txt','w')
    f.close()

    for blast_output in bin_dict.keys():
        bin1=str(blast_output).split('_vs_')[1].split('.txt')[0]
        bin2=str(blast_output).split('_vs_')[0].split('Filtrated_')[1]
        target_bin=bin2
        print('Processing '+str(bin1)+' and '+str(bin2))

        x=0
        for contigs in contigs_coverage.keys():
            x+=1
            if x == 1:
                num=len(contigs_coverage[contigs])

        bin1_cov_ave, bin2_cov_ave, contig_num = {}, {}, 0
        for line in open(blast_output,'r'):
            contig_num+=1
            # try:
            #     contig2=str(line).strip().split('\t')[0].split('|')[0]
            # except:
            contig2=str(line).strip().split('\t')[0]

            # try:
            #     contig1=str(line).strip().split('\t')[1].split('|')[0]
            # except:
            contig1=str(line).strip().split('\t')[1]

            for i in range(1, num+1):
                try:
                    bin1_cov_ave[i]+=float(contigs_coverage[contig1][i])
                    bin2_cov_ave[i]+=float(contigs_coverage[contig2][i])
                except:
                    bin1_cov_ave[i]=float(contigs_coverage[contig1][i])
                    bin2_cov_ave[i]=float(contigs_coverage[contig2][i])

        normalization={}
        for i in bin1_cov_ave.keys():
            bin1_cov_ave[i]=float(bin1_cov_ave[i]/contig_num)
            bin2_cov_ave[i]=float(bin2_cov_ave[i]/contig_num)
            if bin2_cov_ave[i] != 0:
                normalization[i]=bin1_cov_ave[i]/bin2_cov_ave[i]
            else:
                normalization[i]=1
                fabs=open('Abnor_normalization_record.txt','a')
                fabs.write(str(bin1)+'\t'+str(bin2)+'\t'+str(i)+'\t'+str(bin1_cov_ave[i])+'\t'+str(bin2_cov_ave[i])+'\n')
                fabs.close()
        
        if len(normalization) == 0:
            os.system('rm Filtrated_'+str(bin2)+'_vs_'+str(bin1)+'.txt')
        else:
            normalized_bin2_depth={}
            f=open(target_bin+'_bins_depth_comparison.txt','a')
            for i in binset_coverage_total[bin2].keys():
                normalized_bin2_depth[i]=binset_coverage_total[bin2][i]*float(normalization[i])
            f.write(str(bin1)+'\t'+str(binset_coverage_total[bin1])+'\t'+str(bin2)+'\t'+str(normalized_bin2_depth)+'\n')
            f.close()
            os.system('mv Filtrated_'+str(bin2)+'_vs_'+str(bin1)+'.txt '+pwd+'/bin_comparison_folder')

def TNFs_refiner(binset, assembly, coverage_refined_folder, threshold, pwd, num_threads):
    os.system('mkdir '+binset+'_TNFs_outliner')
    os.system('mkdir '+binset+'_TNFs')
    print('Calculating TNFs of '+assembly)
    os.system('calc.kmerfreq.pl -i '+str(assembly)+' -o '+str(assembly)+'.kmer.txt')
    
    os.chdir(pwd+'/'+coverage_refined_folder)
    bin_contigs, bin_contigs_mock, bin_outliner={}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+coverage_refined_folder):
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                if '.fa' in hz:
                    bin_outliner[file]={}
                    bin_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][record.id]=record.seq
                        bin_contigs_mock[file][record.id]=1
    
    os.chdir(pwd+'/'+binset+'_TNFs')
    Bins_TNFs, bin_contig_list={}, {}
    pool=Pool(processes=num_threads)
    for item in bin_contigs.keys():
        pool.apply_async(splitting_kmer_file, args=(item, assembly, bin_contigs_mock, pwd))
    pool.close()
    pool.join()

    for item in bin_contigs.keys():
        Bins_TNFs[item], bin_contig_list[item] = [], []
        for line in open(pwd+'/'+binset+'_TNFs/'+item+'_kmer.txt', 'r'):
            contig=str(line).strip().split('\t')[0]
            lis=str(line).strip().split('\t')
            bin_contig_list[item].append(contig)
            for i in range(1, len(lis)):
                Bins_TNFs[item].append(lis[i])
    # print(str(Bins_TNFs))

    ftnf=open('Bin_TNFs_evaluation.txt','w')
    for item in Bins_TNFs.keys():
        try:
            print('Processing TNFs of bin '+item)
            ftnf.write(str(item))
            num_contig=len(bin_contigs[item])
            ftnf.write('\t'+str(num_contig)+'\t'+str(Bins_TNFs[item]))
            TNF_array=np.array(Bins_TNFs[item]).reshape((num_contig, 256))
            os.chdir(pwd+'/'+binset+'_TNFs_outliner')
            A=PCA_slector(TNF_array, num_contig)
            newData=A[0]
            explained_variance_ratio=A[1]
            bin_outliner=outlier_remover(item, bin_contig_list[item], threshold, newData, explained_variance_ratio, bin_outliner, 'TNF')
            ftnf.write('\t'+'Passed'+'\n')
            os.chdir(pwd)
        except:
            ftnf.write('\t'+'Failed'+'\n')
    
    os.system('mkdir '+str(binset)+'_TNFs_refined')
    os.chdir(str(binset)+'_TNFs_refined')
    outliner_sum=open('Summary_TNFs_outliners.txt','w')
    outliner_sum.write('Bin'+'\t'+'Threshold'+'\t'+'Outliners'+'\n')
    for item in bin_outliner.keys():
        for th in bin_outliner[item].keys():
            if len(bin_outliner[item][th]) != 0:
                outliner_sum.write(str(item)+'\t'+str(th)+'\t'+str(bin_outliner[item][th])+'\n')
    outliner_sum.close()

    bin_outliner2={}
    bin_outliner2=copy.deepcopy(bin_outliner)
    for item in bin_outliner2.keys():
        dereplicate={}
        for th in bin_outliner2[item].keys():
            if len(bin_outliner2[item][th]) not in dereplicate.keys():
                dereplicate[len(bin_outliner2[item][th])]=float(th)
            else:
                if float(th) > dereplicate[len(bin_outliner2[item][th])]:
                    del bin_outliner[item][th]

    for item in bin_contigs.keys():
        for th in bin_outliner[item].keys():
            f=open(item+'_TNFs_'+str(th)+'.fa','w')
            for contig in bin_contigs[item].keys():
                if contig not in bin_outliner[item][th]:
                    f.write('>'+str(contig)+'\n'+str(bin_contigs[item][contig])+'\n')
            f.close()
    os.chdir(pwd)

def genome_contigs_recorder(binset, binset_record, binset_genome_size, coverage_matrix):
    binset_record[str(binset)]={}
    binset_genome_size[str(binset)]={}
    bins_coverage, bin_num_contigs, bin_total_length, bin_GC, bin_GC_ratio, binset_coverage, all_bins={}, {}, {}, {}, {}, {}, {}
    pwd=os.getcwd()

    print('Parsing '+binset)

    n1=0
    for line in open(coverage_matrix, 'r'):
        n1+=1
        if n1 == 1:
            num=str(line).strip().count('drange')
        else:
            ids=str(line).strip().split('\t')[0]
            bins_coverage[str(ids)]={}
            i=1
            while i <= int(num):
                bins_coverage[str(ids)][i]=str(line).strip().split('\t')[3*i+1]
                i+=1

    for root, dirs, files in os.walk(binset):
        os.chdir(pwd+'/'+binset)
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    binset_coverage[file]={}
                    bin_total_length[file]=0
                    bin_GC[file]=0
                    n=0
                    i=1
                    while i <= int(num):
                        binset_coverage[file][i]=0
                        i+=1

                    for record in SeqIO.parse(file,'fasta'):        
                        n+=1
                        bin_num_contigs[file]=n
                        bin_total_length[file]+=len(record.seq)
                        bin_GC[file]+=int(str(record.seq).count('G'))
                        bin_GC[file]+=int(str(record.seq).count('C'))
                        if n == 1:
                            binset_genome_size[str(binset)][str(file)]=len(record.seq)
                        else:
                            binset_genome_size[str(binset)][str(file)]+=len(record.seq)

                        if str(record.id) not in binset_record[str(binset)].keys():
                            binset_record[str(binset)][str(record.id)]=[]
                            binset_record[str(binset)][str(record.id)].append(str(file))
                        else:
                            binset_record[str(binset)][str(record.id)].append(str(file))
                    
                        i=1
                        while i <= int(num):
                            binset_coverage[file][i]+=float(bins_coverage[str(record.id)][i])
                            i+=1
                
                    if 'noclass' not in file and 'unbined' not in file:
                        all_bins[file]=1
            else:
                hz=file.split('.')[-1]
                if hz == 'fasta' or hz == 'fa':
                    #print(str(file))### test
                    binset_coverage[file]={}
                    bin_total_length[file]=0
                    bin_GC[file]=0
                    n=0
                    i=1
                    while i <= int(num):
                        binset_coverage[file][i]=0
                        i+=1

                    for record in SeqIO.parse(file,'fasta'):        
                        n+=1
                        bin_num_contigs[file]=n
                        bin_total_length[file]+=len(record.seq)
                        bin_GC[file]+=int(str(record.seq).count('G'))
                        bin_GC[file]+=int(str(record.seq).count('C'))
                        if n == 1:
                            binset_genome_size[str(binset)][str(file)]=len(record.seq)
                        else:
                            binset_genome_size[str(binset)][str(file)]+=len(record.seq)

                        if str(record.id) not in binset_record[str(binset)].keys():
                            binset_record[str(binset)][str(record.id)]=[]
                            binset_record[str(binset)][str(record.id)].append(str(file))
                        else:
                            binset_record[str(binset)][str(record.id)].append(str(file))
                    
                        i=1
                        while i <= int(num):
                            binset_coverage[file][i]+=float(bins_coverage[str(record.id)][i])
                            i+=1
                
                    if 'noclass' not in file and 'unbined' not in file:
                        all_bins[file]=1
    os.chdir(pwd)

    # print(binset_coverage
    binset_coverage_avg={}
    f=open(str(binset)+'_bins_coverage.txt','w')
    i=1
    a='Bin'
    while i <= int(num):
        a+='\t'+'Coverage'+str(i)
        i+=1
    f.write(str(a)+'\n')

    for item in binset_coverage.keys():
        binset_coverage_avg[item]={}
        a=str(item)
        i=1
        while i <= int(num):
            avg_coverage=float(binset_coverage[item][i])/int(bin_num_contigs[item])
            binset_coverage_avg[item][i]=avg_coverage
            a+='\t'+str(avg_coverage)
            i+=1
        f.write(str(a)+'\n')
    f.close()

    f_t=open(binset+'_gc.txt', 'w')
    for item in bin_GC.keys():
        bin_GC_ratio[item]=round(100*float(bin_GC[item])/float(bin_total_length[item]),1)
        f_t.write(str(item)+'\t'+str(bin_GC_ratio[item])+'\n')
    f_t.close()
    print('---------------------')
    
    return binset_record, binset_genome_size, binset_coverage_avg, bin_GC_ratio, all_bins

def coverage_GC_comparitor(binset_coverage_avg1, binset_coverage_avg2, binset_GC_ratio1, binset_GC_ratio2, iteration_num):
    print('Comparing coverages')
    bins_score, bins_score_total, bins_score_delta, bin_gc={}, {}, {}, {}
    for item in binset_coverage_avg1.keys():
        num=len(binset_coverage_avg1[item])
        i=1
        bins_score[item], bins_score_total[item], bins_score_delta[item]={}, {}, {}
        while i <= num:
            a=str(binset_coverage_avg1[item][i])
            for item2 in binset_coverage_avg2.keys():
                b=str(binset_coverage_avg2[item2][i])
                ave_cov=(float(a)+float(b))/2
                delta_coverage=abs(float(b)-float(a))
                if float(ave_cov) != 0:
                    delta_coverage_vari=round(100*float(delta_coverage)/float(ave_cov), 2) # lower the delta_coverage is, closer the two bins
                else:
                    delta_coverage_vari=round(100*float(delta_coverage)/0.001, 2) # lower the delta_coverage is, closer the two bins
                
                if item2 not in bins_score_total[item].keys():
                    bins_score_total[item][item2]=str(a)+'\t'+str(b)+'\t'+str(delta_coverage_vari)
                else:
                    bins_score_total[item][item2]+='\t'+str(a)+'\t'+str(b)+'\t'+str(delta_coverage_vari)
                    
                if float(a) <= 10 or float(b) <= 10:
                    if item2 not in bins_score[item].keys():
                        bins_score[item][item2]=1
                        bins_score_delta[item][item2]=float(delta_coverage_vari) 
                    else:
                        bins_score[item][item2]+=1
                        bins_score_delta[item][item2]+=float(delta_coverage_vari)
                # elif ave_cov > 1 and ave_cov <= 5:
                #     if delta_coverage_vari <= 60: ### variable percentage
                #         if item2 not in bins_score[item].keys():
                #             bins_score[item][item2]=1
                #             bins_score_delta[item][item2]=float(delta_coverage_vari) 
                #         else:
                #             bins_score[item][item2]+=1
                #             bins_score_delta[item][item2]+=float(delta_coverage_vari)
                # elif ave_cov > 5 and ave_cov <= 15:
                #     if delta_coverage_vari <= 30: ### variable percentage
                #         if item2 not in bins_score[item].keys():
                #             bins_score[item][item2]=1
                #             bins_score_delta[item][item2]=float(delta_coverage_vari) 
                #         else:
                #             bins_score[item][item2]+=1
                #             bins_score_delta[item][item2]+=float(delta_coverage_vari)
                # elif ave_cov > 15 and ave_cov <= 30:
                #     if delta_coverage_vari <= 20: ### variable percentage
                #         if item2 not in bins_score[item].keys():
                #             bins_score[item][item2]=1
                #             bins_score_delta[item][item2]=float(delta_coverage_vari) 
                #         else:
                #             bins_score[item][item2]+=1
                #             bins_score_delta[item][item2]+=float(delta_coverage_vari)
                else:
                    if delta_coverage_vari <= 100: ### variable percentage
                        if item2 not in bins_score[item].keys():
                            bins_score[item][item2]=1
                            bins_score_delta[item][item2]=float(delta_coverage_vari) 
                        else:
                            bins_score[item][item2]+=1
                            bins_score_delta[item][item2]+=float(delta_coverage_vari)
            i+=1
        
        del_item=[]
        for item2 in bins_score[item].keys():
            if int(bins_score[item][item2]) < num:
                del_item.append(item2)
        
        if len(del_item) != 0:
            for item2 in del_item:
                del bins_score[item][item2]
                del bins_score_delta[item][item2]

    binset_coverage_gc={}
    for item in binset_GC_ratio1.keys():
        for item2 in binset_GC_ratio2.keys():
            delta_GC_ratio=round(100*abs(float(binset_GC_ratio2[item2])-float(binset_GC_ratio1[item]))/float(binset_GC_ratio1[item]),2)
            # delta_GC_ratio=round(100*float(delta_GC)/float(binset_GC_ratio1[item]),1)
            if delta_GC_ratio <= 5: ## +- 3% total 6% var
                if item2 in bins_score[item].keys():
                    if str(item) not in binset_coverage_gc.keys():
                        binset_coverage_gc[str(item)]=str(item2)+':'+str(bins_score_delta[item][item2])+':'+str(delta_GC_ratio)
                    else:
                        binset_coverage_gc[str(item)]+='\t'+str(item2)+':'+str(bins_score_delta[item][item2])+':'+str(delta_GC_ratio)

                if item not in bin_gc.keys():
                    bin_gc[item]=str(item)+':'+str(binset_GC_ratio1[item])+'\t'+str(item2)+':'+str(binset_GC_ratio2[item2])+':'+str(delta_GC_ratio)
                else:
                    bin_gc[item]+='\t'+str(item2)+':'+str(binset_GC_ratio2[item2])+':'+str(delta_GC_ratio)

    bins_coverage_score={}
    f_bin_coverage=open('Bins_similar_coverage_'+str(iteration_num)+'.txt', 'w')
    f_bin_coverage.write('Target Bin'+'\t'+'Similar Bin'+'\t'+'Target Bin Coverage'+'\t'+'Similar Bin Coverage'+'\t'+'Average Coverage Variation'+'\n')
    for item in bins_score.keys():
        for item2 in bins_score_delta[item].keys():
            average_coverage_var=round(float(bins_score_delta[item][item2]/num),2)
            # if item in 
            f_bin_coverage.write(str(item)+'\t'+str(item2)+'\t'+str(average_coverage_var)+'\n')
            bins_coverage_score[str(item)+'\t'+str(item2)]=str(average_coverage_var)
    f_bin_coverage.close()

    f_bin_coverage=open('Total_bins_similar_coverage_'+str(iteration_num)+'.txt', 'w')
    f_bin_coverage.write('Target Bin'+'\t'+'Similar Bin'+'\t'+'Target Bin Coverage'+'\t'+'Similar Bin Coverage'+'\t'+'Average Coverage Variation'+'\n')
    for item in bins_score_total.keys():
        for item2 in bins_score_total[item].keys():
            f_bin_coverage.write(str(item)+'\t'+str(item2)+'\t'+str(bins_score_total[item][item2])+'\n')
    f_bin_coverage.close()

    f_gc=open('Bins_similar_GC_'+str(iteration_num)+'.txt', 'w')
    f_gc.write('Group'+'\t'+'Bin:GC%:deltaGC%'+'\n')
    n=0
    for item in bin_gc.keys():
        n+=1
        f_gc.write(str(n)+'\t'+str(bin_gc[item])+'\n')
    f_gc.close()

    f_co=open('Bins_similar_co_GC_coverage_'+str(iteration_num)+'.txt', 'w')
    f_co.write('Target Bin'+'\t'+'Similar:Coverage_var(%):GC_var(%)'+'\n')
    n=0
    for item in binset_coverage_gc.keys():
        n+=1
        f_co.write(str(item)+'\t'+str(binset_coverage_gc[item])+'\n')
    f_co.close()
    print('-----------------------')
    return bins_coverage_score, bin_gc, binset_coverage_gc

def seq_comparitor(blast_output, binset1, binset2, seq_len1, seq_len2, binset_record1, binset_record2, binset_genome_size1, binset_genome_size2, bins_coverage_score):
    print('Comparing blast output of '+binset1+' and '+binset2)
    print('-----------------------')
    contig_contig_score={}
    contig_contig_precentage={}
    contig_bin_genome={}
    contig_contig_score_mock={}

    binset_record1_mock={}
    binset_record1_mock[str(binset1)]={}
    binset_record2_mock={}
    binset_record2_mock[str(binset2)]={}
    for item in binset_record1[str(binset1)].keys():
        binset_record1_mock[str(binset1)][item]=0
    for item in binset_record2[str(binset2)].keys():
        binset_record2_mock[str(binset2)][item]=0

    # ft=open('Test_Raw_Bins_Comparison_'+blast_output, 'w')
    # print(binset_record1
    print('Parsing '+blast_output)
    for line in open(blast_output,'r'):
        simi=str(line).strip().split('\t')[2]
        length=eval(str(line).strip().split('\t')[3])
        # if float(simi) >= 99 and int(length) >= 500:
            # ID1_ORFs=str(line).strip().split('\t')[1] ### Be aware of the difference between ID1 and ID2. Seq1 has been set as database. 
            # ID1=ID1_ORFs.split('_')[0]
            # ID2_ORFs=str(line).strip().split('\t')[0]
            # ID2=ID2_ORFs.split('_')[0]
        ID1=str(line).strip().split('\t')[1]
        ID2=str(line).strip().split('\t')[0]
            
        Alignment_length=int(int(length)*float(simi)/100)
        ID1_aligned_percentage=round(float(int(length)*float(simi)/seq_len1[ID1]),2)
        ID2_aligned_percentage=round(float(int(length)*float(simi)/seq_len2[ID2]),2)
        num=len(contig_contig_score_mock)
        contig_contig_score_mock[str(ID1)+'\t'+str(ID2)]=1
        if len(contig_contig_score_mock) == num+1:
            contig_contig_score[str(ID1)+'\t'+str(ID2)]=Alignment_length
            contig_contig_precentage[str(ID1)+'\t'+str(ID2)]=str(ID1_aligned_percentage)+'\t'+str(ID2_aligned_percentage)
        else:
            contig_contig_score[str(ID1)+'\t'+str(ID2)]+=Alignment_length
            ID1_aligned_percentage_2=float(contig_contig_precentage[str(ID1)+'\t'+str(ID2)].split('\t')[0])+ID1_aligned_percentage
            ID2_aligned_percentage_2=float(contig_contig_precentage[str(ID1)+'\t'+str(ID2)].split('\t')[1])+ID2_aligned_percentage
            contig_contig_precentage[str(ID1)+'\t'+str(ID2)]=str(ID1_aligned_percentage_2)+'\t'+str(ID2_aligned_percentage_2)

        bin_name=[]
        m1=len(binset_record1_mock[str(binset1)])
        binset_record1_mock[str(binset1)][str(ID1)]=0
        if len(binset_record1_mock[str(binset1)]) == m1:
            m2=len(binset_record2_mock[str(binset2)])
            binset_record2_mock[str(binset2)][str(ID2)]=0
            if len(binset_record2_mock[str(binset2)]) == m2:
                bin1_list=binset_record1[str(binset1)][str(ID1)]
                bin2_list=binset_record2[str(binset2)][str(ID2)]
                for item in bin1_list:
                    xx=0
                    while xx < len(bin2_list):
                        bin_name.append(str(item)+'\t'+bin2_list[xx])
                        xx+=1

                    # bin1=str(binset_record1[str(binset1)][str(ID1)])
                    # bin2=str(binset_record2[str(binset2)][str(ID2)])
                    # if ',' not in bin1 and ',' not in bin2:
                    #     bin_name.append(bin1+'\t'+bin2)
                    # elif ',' in bin1 and ',' not in bin2:
                    #     lis=str(bin1).split(',')
                    #     xx=0
                    #     while xx < len(lis):
                    #         bin_name.append(str(lis[xx])+'\t'+bin2)
                    #         xx+=1
                    # elif ',' not in bin1 and ',' in bin2:
                    #     lis=bin2.split(',')
                    #     xx=0
                    #     while xx < len(lis):
                    #         bin_name.append(bin1+'\t'+bin2.split(',')[xx])
                    #         xx+=1
                    # else:
                    #     lis1=bin1.split(',')
                    #     lis2=bin2.split(',')
                    #     for item in lis1:
                    #         xx=0
                    #         while xx < len(lis2):
                    #             bin_name.append(str(item)+'\t'+bin2.split(',')[xx])
                    #             xx+=1
                    
                for item in bin_name:
                    # ft.write(str(item)+'\n')
                    if str(item) not in contig_bin_genome.keys():
                        contig_bin_genome[str(item)]=Alignment_length
                    else:
                        contig_bin_genome[str(item)]+=Alignment_length
            else:
                del binset_record2_mock[str(binset2)][str(ID2)]
        else:
            del binset_record1_mock[str(binset1)][str(ID1)]

   # ft.close()

    print('Sorting contigs score')
    contig_contig_score=sorted(contig_contig_score.items(), key=lambda d:d[0])

    f=open('Contigs_scoring_'+blast_output, 'w')
    f.write('Contig1'+'\t'+'Contig2'+'\t'+'Alignment_length'+'\t'+'Contig1_aligned_percentage'+'\t'+'Contig2_aligned_percentage'+'\n')
    for item in contig_contig_score:
        f.write(str(item[0])+'\t'+str(item[1])+'\t'+str(contig_contig_precentage[item[0]])+'\n')
    f.close()

    bin_seq_coverage_filtrated={}
    f=open('Raw_Bins_Comparison_'+blast_output, 'w')
    f3=open('Filtrated_Bins_Seq_Similarity_Coverage_'+blast_output, 'w')
    title='GenomeA'+'\t'+'GnomesB'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'
    f.write(title+'\n')
    f3.write(title+'\t'+'Average_Coverage_Variation'+'\n')
    for item in contig_bin_genome.keys():
        # bin1=item.split('\t')[0]
        # bin2=item.split('\t')[1]
        num=contig_bin_genome[item]
        A_percentage=round(100*float(num)/float(binset_genome_size1[str(binset1)][item.split('\t')[0]]),3)
        B_percentage=round(100*float(num)/float(binset_genome_size2[str(binset2)][item.split('\t')[1]]),3)
        contig_bin_genome[item]=str(num)+'\t'+str(A_percentage)+'\t'+str(B_percentage)
        f.write(str(item)+'\t'+str(num)+'\t'+str(A_percentage)+'\t'+str(B_percentage)+'\n')
        
        Sum_sim=float(A_percentage)+float(B_percentage)
        if 'noclass' not in item and 'unbined' not in item and num >= 500000:
            if Sum_sim >= 80:
                if float(A_percentage) >= 50 or float(B_percentage) >= 50:
                    if str(item) in bins_coverage_score.keys():
                        # f2.write(str(item[0])+'\t'+str(item[1])+'\n')
                        bin_seq_coverage_filtrated[str(item)]=str(num)+'\t'+str(A_percentage)+'\t'+str(B_percentage)+'\t'+str(bins_coverage_score[str(item)])
                        f3.write(str(item)+'\t'+str(num)+'\t'+str(A_percentage)+'\t'+str(B_percentage)+'\t'+str(bins_coverage_score[str(item)])+'\n')
                    else:
                        continue
    f.close()
    f3.close()

    # contig_bin_genome=sorted(contig_bin_genome.items(), key=lambda d:d[0])

    # bin_seq_coverage_filtrated={}
    # f=open('Raw_Bins_Comparison_'+blast_output, 'w')
    # # f2=open('Filtrated_Bins_Comparison_'+blast_output, 'w')
    # f3=open('Filtrated_Bins_Seq_Similarity_Coverage_'+blast_output, 'w')
    # title='GenomeA'+'\t'+'GnomesB'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'
    # f.write(title+'\n')
    # # f2.write(title+'\n')
    # f3.write(title+'\t'+'Average_Coverage_Variation'+'\n')
    # for item in contig_bin_genome:
    #     f.write(str(item[0])+'\t'+str(item[1])+'\n') ###
        # print(str(item[1])
    #     alignment_length=int(str(item[1].split('\t')[0]))
    #     A_sim=item[1].split('\t')[1]
    #     B_sim=item[1].split('\t')[2]
    #     Sum_sim=float(A_sim)+float(B_sim)
    #     # if 'noclass' not in item[0] and 'unbined' not in item[0] and alignment_length >= 500000 and Sum_sim >= 100: ### Similar sequence which is larger than 100 Kbp will be kept to consider
    #     if 'noclass' not in item[0] and 'unbined' not in item[0] and alignment_length >= 500000:
    #         if Sum_sim >= 80:
    #             if float(A_sim) >= 50 or float(B_sim) >= 50:
    #                 if str(item[0]) in bins_coverage_score.keys():
    #                     # f2.write(str(item[0])+'\t'+str(item[1])+'\n')
    #                     bin_seq_coverage_filtrated[str(item[0])]=str(item[1])+'\t'+str(bins_coverage_score[str(item[0])])
    #                     f3.write(str(item[0])+'\t'+str(item[1])+'\t'+str(bins_coverage_score[str(item[0])])+'\n')
    #                 else:
    #                     continue
    # f.close()
    # # f2.close()
    # f3.close()
    return contig_contig_score, contig_bin_genome, bin_seq_coverage_filtrated

def checkm_connections(binset):
    print('Reading checkm output of '+binset)
    pwd=os.getcwd()
    binset_checkm_connection={}
    for root, dirs, files in os.walk(binset):
        os.chdir(pwd+'/'+binset)
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    bins=str(line).strip().split('\t')[0]
                    bin_id=str(bins)+'.fa'
                    bin_id2=str(bins)+'.fasta'
                    binset_checkm_connection[str(bin_id)]={}
                    binset_checkm_connection[str(bin_id)]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\'')[0]
                    binset_checkm_connection[str(bin_id)]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id)]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
                    binset_checkm_connection[str(bin_id)]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id2)]={}
                    binset_checkm_connection[str(bin_id2)]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id2)]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\'')[0]
                    binset_checkm_connection[str(bin_id2)]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id2)]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0])
                    binset_checkm_connection[str(bin_id2)]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
                    binset_checkm_connection[str(bin_id2)]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0])
    os.chdir(pwd)
    print('Done of reading checkm output of '+binset)
    print('-----------------------------------------')
    return binset_checkm_connection

def bin_comparitor(bin_seq_coverage_filtrated, binset_checkm_connection_1, binset_checkm_connection_2, num):
    print('Comparing bins')
    selected_bins, eliminated_bins={}, {}
    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}
    f=open('Selected_bins_'+str(num)+'.txt', 'w')
    f2=open('Removed_bins_'+str(num)+'.txt', 'w')
    try:
        f3=open('Bins_wth_sameQua_difSize.txt', 'a')
    except:
        f3=open('Bins_wth_sameQua_difSize.txt', 'w')
    f.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Bin1'+'\t'+'Bin1 Cov.'+'\t'+'Bin2'+'\t'+'Bin2 Cov.'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    f2.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Bin1'+'\t'+'Bin1 Cov.'+'\t'+'Bin2'+'\t'+'Bin2 Cov.'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')

    del_filtration={}
    for item in bin_seq_coverage_filtrated.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]

        bin_score_1=marker_score[binset_checkm_connection_1[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[binset_checkm_connection_2[ID2]['marker lineage'].split('__')[0]]

        delta_1=binset_checkm_connection_1[ID1]['Completeness']-binset_checkm_connection_1[ID1]['Contamination']
        delta_2=binset_checkm_connection_2[ID2]['Completeness']-binset_checkm_connection_2[ID2]['Contamination']

        if binset_checkm_connection_1[ID1]['Contamination'] >= 12 and float(delta_1) <= 80:
            eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
            del_filtration[item]=0
        
        if binset_checkm_connection_2[ID2]['Contamination'] >= 12 and float(delta_2) <= 80:
            eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
            del_filtration[item]=0

    for item in del_filtration.keys():
        del bin_seq_coverage_filtrated[item]

    for item in bin_seq_coverage_filtrated.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]
        bin_score_1=marker_score[binset_checkm_connection_1[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[binset_checkm_connection_2[ID2]['marker lineage'].split('__')[0]]
        CPN_CTN_1=binset_checkm_connection_1[ID1]['Completeness']-5*binset_checkm_connection_1[ID1]['Contamination']
        CPN_CTN_2=binset_checkm_connection_2[ID2]['Completeness']-5*binset_checkm_connection_2[ID2]['Contamination']
        bin1_size=binset_checkm_connection_1[ID1]['Genome size']
        bin2_size=binset_checkm_connection_2[ID2]['Genome size']
        bin1_mean_len=binset_checkm_connection_1[ID1]['Mean scaffold length']
        bin2_mean_len=binset_checkm_connection_2[ID2]['Mean scaffold length']       

        if bin_score_1 == bin_score_2:    
            if CPN_CTN_1 == CPN_CTN_2:
                if int(bin1_size) == int(bin2_size) and float(bin1_mean_len) == float(bin2_mean_len):
                    eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                    if ID1 in eliminated_bins.keys():
                        del eliminated_bins[ID1]
                    selected_bins[ID1]=binset_checkm_connection_1[ID1]
                    f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
                else:
                    ### Record bins with the same quality value but different size or medium length
                    f3.write(ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
                
            elif CPN_CTN_1 > CPN_CTN_2:
                eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                if ID1 in eliminated_bins.keys():
                    del eliminated_bins[ID1]
                selected_bins[ID1]=binset_checkm_connection_1[ID1]
                f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
            else:
                eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
                if ID2 in eliminated_bins.keys():
                    del eliminated_bins[ID2]
                selected_bins[ID2]=binset_checkm_connection_2[ID2]
                f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
        else:
            if float(bin_score_1)*float(CPN_CTN_1) >= float(bin_score_2)*float(CPN_CTN_2):
                eliminated_bins[ID2]=binset_checkm_connection_2[ID2]
                if ID1 in eliminated_bins.keys():
                    del eliminated_bins[ID1]
                selected_bins[ID1]=binset_checkm_connection_1[ID1]
                f.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n')
            else:
                eliminated_bins[ID1]=binset_checkm_connection_1[ID1]
                if ID2 in eliminated_bins.keys():
                    del eliminated_bins[ID2]
                selected_bins[ID2]=binset_checkm_connection_2[ID2]   
                f.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(bin_seq_coverage_filtrated[item])+'\t'+str(binset_checkm_connection_1[ID1])+'\t'+str(binset_checkm_connection_2[ID2])+'\n') 
    f.close()
    f3.close()

    # for item in selected_bins.keys():
    #     if item in eliminated_bins.keys():
    #         del selected_bins[item]
    for item in eliminated_bins.keys():
        f2.write(item+'\t'+str(eliminated_bins[item])+'\n')
    f2.close()

    print('Comparison done!')
    print('----------------')
    return selected_bins, eliminated_bins

def new_selected_bins_generator(selected_bins, eliminated_bins, all_bins_1, all_bins_2, binset_checkm_connection_1, binset_checkm_connection_2, iteration_num, binset1, binset2, Contigs_assembly1, Contigs_assembly2, coverage_matrix1, coverage_matrix2, binset_coverage_avg1, binset_coverage_avg2):
    pwd=os.getcwd()
    extract_bins_1, extract_bins_2, coverage_maxtrix={}, {}, {}

    nx=0
    for line in open(coverage_matrix1, 'r'):
        nx+=1
        if nx == 1:
            coverage_title=str(line).strip()
        else:
            coverage_maxtrix[str(line).strip().split('\t')[0]]=str(line).strip()
    
    nx=0
    for line in open(coverage_matrix2, 'r'):
        nx+=1
        if nx > 1:
            coverage_maxtrix[str(line).strip().split('\t')[0]]=str(line).strip()

    for item in all_bins_1.keys():
        if item not in selected_bins.keys() and item not in eliminated_bins.keys():
            delta_1=binset_checkm_connection_1[item]['Completeness']-binset_checkm_connection_1[item]['Contamination']
            if binset_checkm_connection_1[item]['Contamination'] >= 12 and float(delta_1) <= 80:
                eliminated_bins[item]=binset_checkm_connection_1[item]
            else:
                extract_bins_1[item]=binset_checkm_connection_1[item]
    
    for item in all_bins_2.keys():
        if item not in selected_bins.keys() and item not in eliminated_bins.keys():
            delta_2=binset_checkm_connection_2[item]['Completeness']-binset_checkm_connection_2[item]['Contamination']
            if binset_checkm_connection_2[item]['Contamination'] >= 12 and float(delta_2) <= 80:
                eliminated_bins[item]=binset_checkm_connection_2[item]
            else:
                extract_bins_2[item]=binset_checkm_connection_2[item]

    f=open('Extract_bins_'+str(iteration_num)+'.txt','w')
    f.write('Extract bins from binset 1'+'\n')
    for item in extract_bins_1.keys():
        f.write(str(item)+'\t'+str(extract_bins_1[item])+'\n')

    f.write('Extract bins from binset 2'+'\n')
    for item in extract_bins_2.keys():
        f.write(str(item)+'\t'+str(extract_bins_2[item])+'\n')
    f.close()

    try:
        os.mkdir('Iteration_'+str(iteration_num)+'_genomes')
    except:
        os.system('rm -rf Iteration_'+str(iteration_num)+'_genomes')
        os.mkdir('Iteration_'+str(iteration_num)+'_genomes')
        print('Iteration_'+str(iteration_num)+'_genomes exist')
        print('Re-created folder of Iteration_'+str(iteration_num)+'_genomes')

    new_Contigs, new_Contig_mock, bin_avg_cov_new={}, {}, {}
    for root, dirs, files in os.walk(binset1):
        os.chdir(binset1)
        for file in files:
            if file in extract_bins_1.keys() or file in selected_bins.keys():
                if file not in eliminated_bins.keys():
                    print( binset_coverage_avg1[file] )

                    bin_avg_cov_new[file]=binset_coverage_avg1[file]
                    os.system('cp '+file+' '+pwd+'/Iteration_'+str(iteration_num)+'_genomes')
                    for record in SeqIO.parse(file, 'fasta'):
                        new_Contigs[record.id]=str(record.seq)
                        new_Contig_mock[record.id]=0

    for root, dirs, files in os.walk(pwd+'/'+binset2):
        os.chdir(pwd+'/'+binset2)
        for file in files:
            if file in extract_bins_2.keys() or file in selected_bins.keys():
                if file not in eliminated_bins.keys():
                    print( binset_coverage_avg2[file] )

                    bin_avg_cov_new[file]=binset_coverage_avg2[file]
                    os.system('cp '+file+' '+pwd+'/Iteration_'+str(iteration_num)+'_genomes')
                    for record in SeqIO.parse(file, 'fasta'):
                        new_Contigs[record.id]=str(record.seq)
                        new_Contig_mock[record.id]=0

    os.chdir(pwd)                    
    new_contigs_name='Contigs_iteration_'+str(iteration_num)+'.fa'
    new_contigs=open(new_contigs_name, 'w')
    new_coverage_name='Coverage_matrix_for_binning_iteration_'+str(iteration_num)+'.txt'
    new_coverage=open(new_coverage_name,'w')
    new_coverage.write(str(coverage_title)+'\n')
    
    for item in new_Contigs.keys():
        new_contigs.write('>'+str(item)+'\n'+str(new_Contigs[item])+'\n')
        new_coverage.write(str(coverage_maxtrix[str(item)])+'\n')
    new_coverage.close()
    new_contigs.close()

    # new_ORFs=ORFs_predictor(str(new_contigs_name))

    # new_ORFs=open('ORFs_iteration_'+str(iteration_num)+'.fa','w')
    # id_num=len(new_Contig_mock)
    # for record in SeqIO.parse(ORFs_assembly1, 'fasta'):
    #     new_Contig_mock[str(record.id).split('_')[0]]=0
    #     if len(new_Contig_mock) == id_num:
    #     # if str(record.id).split('_')[0] in new_Contigs.keys():
    #        new_ORFs.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    #     else:
    #         del new_Contig_mock[str(record.id).split('_')[0]]
    
    # for record in SeqIO.parse(ORFs_assembly2, 'fasta'):
    #     new_Contig_mock[str(record.id).split('_')[0]]=0
    #     if len(new_Contig_mock) == id_num:
    #     # if str(record.id).split('_')[0] in new_Contigs.keys():
    #        new_ORFs.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    #     else:
    #         del new_Contig_mock[str(record.id).split('_')[0]]
    # new_ORFs.close()
    
    
    high_quality_bins,all_remain_bins={},{}
    os.chdir(pwd+'/Iteration_'+str(iteration_num)+'_genomes')
    f=open('Iteration_'+str(iteration_num)+'_bin_stats_ext.tsv','w')
    for item in selected_bins.keys():
        if item not in eliminated_bins.keys():
            all_remain_bins[str(item)]=str(selected_bins[item])
            #if '_metabat_genomes.' in str(item):
            #    bin_name=str(item).replace('.fa', '')
            #elif '_maxbin2_genomes.' in str(item):
            #    bin_name=str(item).replace('.fasta', '')  
            #else:
            #    continue
            item_list=item.split('.')
            item_list.remove(item_list[-1])
            bin_name='.'.join(item_list)
            f.write(str(bin_name)+'\t'+str(selected_bins[item])+'\n')

    for item in extract_bins_1.keys():
        all_remain_bins[str(item)]=str(extract_bins_1[item])
        item_list=item.split('.')
        item_list.remove(item_list[-1])
        bin_name='.'.join(item_list)
        f.write(str(bin_name)+'\t'+str(extract_bins_1[item])+'\n')

    for item in extract_bins_2.keys():
        all_remain_bins[str(item)]=str(extract_bins_2[item])
        item_list=item.split('.')
        item_list.remove(item_list[-1])
        bin_name='.'.join(item_list)
        f.write(str(bin_name)+'\t'+str(extract_bins_2[item])+'\n')
    f.close()

    f=open('Medium_quality_bins_CPN-5CTN_50_iteration_'+str(iteration_num)+'_.txt','w')
    f2=open('High_quality_bins_CPN-5CTN_70_iteration_'+str(iteration_num)+'_.txt','w')
    for item in all_remain_bins.keys():
        cpn=str(all_remain_bins[item]).split('Completeness\': ')[1].split(',')[0]
        # if str(all_remain_bins[item]).count(',') == 4:
        ctn=str(all_remain_bins[item]).split('Contamination\': ')[1].split('}')[0]
        quality=float(cpn)-5*float(ctn)
        if quality >= 50:
            f.write(str(item)+'\t'+str(all_remain_bins[str(item)])+'\n')
        else:
            continue
        
        if quality >= 70:
            f2.write(str(item)+'\t'+str(all_remain_bins[str(item)])+'\n')
    f.close()
    f2.close()
    os.chdir(pwd)
    return extract_bins_1, extract_bins_2, 'Contigs_iteration_'+str(iteration_num)+'.fa', 'Iteration_'+str(iteration_num)+'_genomes', new_coverage_name, bin_avg_cov_new

def final_iteration_mapping(contigs, datasets, num_threads, pwd):
    print('Building '+str(contigs)+' Bowtie2 index')
    os.system('bowtie2-build '+str(contigs)+' '+str(contigs))
    for i in range(1, len(datasets)+1):
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(contigs)+' -1 '+str(datasets[str(i)][0])+' -2 '+str(datasets[str(i)][1])+' -S '+str(contigs)+'-'+str(i)+'.sam -q --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(contigs)+'-'+str(i)+'.sam -o '+str(contigs)+'-'+str(i)+'.bam')
        # py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(contigs)+'-'+str(i)+'.bam '+str(contigs)+'-'+str(i)+'_sorted')
        try:
            with open(str(contigs)+'-'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(contigs)+'-'+str(i)+'_sorted.bam '+str(contigs)+'-'+str(i)+'.bam' )
    
        if i == 1:
            bam_sorted=str(contigs)+'-1_sorted.bam'
        else:
            bam_sorted+=' '+str(contigs)+'-'+str(i)+'_sorted.bam'
    
    os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(contigs)+'.depth.txt '+str(bam_sorted))
    # os.system(str(pwd)+'/jgi_summarize_bam_contig_depths --outputDepth '+str(contigs)+'.depth.txt '+str(bam_sorted))

    for i in range(1, len(datasets)+1):
        os.system('rm '+str(contigs)+'-'+str(i)+'.sam '+str(contigs)+'-'+str(i)+'_sorted.bam '+str(contigs)+'-'+str(i)+'.bam')

    bin_depth, bin_contigs_depth, bin_core_contigs, n={}, {}, {}, 0
    for line in open(str(contigs)+'.depth.txt','r'):
        n+=1
        if n >= 2:
            bin_id=str(line).strip().split('\t')[0].split('|')[0].strip()
            bin_depth[bin_id], bin_contigs_depth[bin_id], bin_core_contigs[bin_id] = {}, {}, {}

    n=0
    for line in open(str(contigs)+'.depth.txt','r'):
        n+=1
        if n == 1:
            coverage_num=depth_sample_count(line)
        else:
            bin_id=str(line).strip().split('\t')[0].split('|')[0].strip()
            contig_id_t=str(line).strip().split('\t')[0]
            if '||' not in contig_id_t:
                contig_id=contig_id_t.split('|')[1].strip()
            else:
                contig_id_ta=contig_id_t.replace('||','**')
                contig_id=contig_id_ta.split('|')[1].strip().replace('**','||')

            bin_contigs_depth[bin_id][contig_id]={}
            bin_core_contigs[bin_id][contig_id]={}
            # contig_len=int(str(line).strip().split('\t')[1])
            for i in range(0, coverage_num):
                depth_value=float(str(line).strip().split('\t')[i*2+3])
                i2=i+1
                bin_contigs_depth[bin_id][contig_id][i2]=depth_value
                bin_core_contigs[bin_id][contig_id][i2]=depth_value
                try:
                    bin_depth[bin_id][i2].append(depth_value)
                except:
                    bin_depth[bin_id][i2]=[]
                    bin_depth[bin_id][i2].append(depth_value)

    bin_core_total, bin_core_avg = {}, {}
    for bin_id in bin_depth.keys():
        bin_core_total[bin_id], bin_core_avg[bin_id] = {}, {}
        x25, x75, x25_dict, x75_dict = int(0.25*n), int(0.75*n), {}, {}
        if len(bin_contigs_depth[bin_id]) > 10:
            for i in range(1, len(bin_depth[bin_id])+1):
                depth_list=[y for y in bin_depth[bin_id][i]]
                depth_list.sort()
                x25_dict[i], x75_dict[i], x = 0, 0, 0
                for contig_value in depth_list:
                    x+=1
                    if x == x25:
                        x25_dict[i]=contig_value
                    elif x <= x75:
                        x75_dict[i]=contig_value

            for contig in bin_contigs_depth[bin_id].keys():
                for i in bin_contigs_depth[bin_id][contig].keys():
                    if float(bin_contigs_depth[bin_id][contig][i]) < x25_dict[i] and float(bin_contigs_depth[bin_id][contig][i]) > x75_dict[i]:
                        del bin_core_contigs[bin_id][contig]

        y =  {}
        for i in range(1, len(bin_depth[bin_id])+1):
            bin_core_total[bin_id][i]=0
            y[i]=0
            for contig in bin_core_contigs[bin_id].keys():
                y[i]+=1
                bin_core_total[bin_id][i]+=bin_core_contigs[bin_id][contig][i]

        for i in bin_core_total[bin_id].keys():
            bin_core_avg[bin_id][i]=bin_core_total[bin_id][i]/y[i]
    return bin_core_avg

def mapping(bin, datasets, num_threads, pwd):
    print('Building '+str(bin)+' Bowtie2 index')
    os.system('bowtie2-build '+str(bin)+' '+str(bin))
    for i in range(1, len(datasets)+1):
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(bin)+' -1 '+str(datasets[str(i)][0])+' -2 '+str(datasets[str(i)][1])+' -S '+str(bin)+'-'+str(i)+'.sam -q --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(bin)+'-'+str(i)+'.sam -o '+str(bin)+'-'+str(i)+'.bam')
        # py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(bin)+'-'+str(i)+'.bam '+str(bin)+'-'+str(i)+'_sorted')
        try:
            with open(str(bin)+'-'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(bin)+'-'+str(i)+'_sorted.bam '+str(bin)+'-'+str(i)+'.bam' )
    
        if i == 1:
            bam_sorted=str(bin)+'-1_sorted.bam'
        else:
            bam_sorted+=' '+str(bin)+'-'+str(i)+'_sorted.bam'
    
    os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(bin)+'.depth.txt '+str(bam_sorted))

    for i in range(1, len(datasets)+1):
        os.system('rm '+str(bin)+'-'+str(i)+'.sam '+str(bin)+'-'+str(i)+'_sorted.bam '+str(bin)+'-'+str(i)+'.bam')

    bin_depth, bin_contigs_depth, bin_core_contigs, n={}, {}, {}, 0    
    for line in open(str(bin)+'.depth.txt','r'):
        n+=1
        if n == 1:
            coverage_num=depth_sample_count(line)
        else:
            contig_id=str(line).strip().split('\t')[0]
            bin_contigs_depth[contig_id]={}
            bin_core_contigs[contig_id]={}
            # contig_len=int(str(line).strip().split('\t')[1])
            for i in range(0, coverage_num):
                depth_value=float(str(line).strip().split('\t')[i*2+3])
                i2=i+1
                bin_contigs_depth[contig_id][i2]=depth_value
                bin_core_contigs[contig_id][i2]=depth_value
                try:
                    bin_depth[i2].append(depth_value)
                except:
                    bin_depth[i2]=[]
                    bin_depth[i2].append(depth_value)
    
    x25, x75, x25_dict, x75_dict = int(0.25*n), int(0.75*n), {}, {}
    if n > 10:
        for i in range(1, len(bin_depth)+1):
            depth_list=[y for y in bin_depth[i]]
            depth_list.sort()
            x25_dict[i], x75_dict[i], x = 0, 0, 0
            for contig_value in depth_list:
                x+=1
                if x == x25:
                    x25_dict[i]=contig_value
                elif x <= x75:
                    x75_dict[i]=contig_value

        for contig in bin_contigs_depth.keys():
            for i in bin_contigs_depth[contig].keys():
                if float(bin_contigs_depth[contig][i]) < x25_dict[i] and float(bin_contigs_depth[contig][i]) > x75_dict[i]:
                    del bin_core_contigs[contig]

    bin_core_total, y =  {}, {}
    for i in range(1, len(bin_depth)+1):
        bin_core_total[i]=0
        y[i]=0
        for contig in bin_core_contigs.keys():
            y[i]+=1
            bin_core_total[i]+=bin_core_contigs[contig][i]
    
    bin_core_avg={}
    for i in bin_core_total:
        bin_core_avg[i]=bin_core_total[i]/y[i]
    return bin_core_avg

def parse_checkm(file):
    bin_checkm={}
    for line in open(file,'r'):
        bin_id=str(line).strip().split('\t')[0]
        bin_id_f=bin_id+'.fa'
        bin_id_f2=bin_id+'.fasta'
        bin_id_f3=bin_id.split('.')[0].split('_')[0]+'.fa'
        bin_id_f4=bin_id.split('.')[0].split('_')[0]+'.fasta'
        bin_checkm[bin_id_f]={}
        try:
            bin_checkm[bin_id_f]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f]['Connections']=0
        bin_checkm[bin_id_f]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
        bin_checkm[bin_id_f]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
        try:
            bin_checkm[bin_id_f]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f]['Genome size']=0

        try:
            bin_checkm[bin_id_f]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f]['Mean scaffold length']=0
        bin_checkm[bin_id_f]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())

        bin_checkm[bin_id_f2]={}
        try:
            bin_checkm[bin_id_f2]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f2]['Connections']=0
        bin_checkm[bin_id_f2]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
        bin_checkm[bin_id_f2]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
        try:
            bin_checkm[bin_id_f2]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f2]['Genome size']=0

        try:
            bin_checkm[bin_id_f2]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f2]['Mean scaffold length']=0
        bin_checkm[bin_id_f2]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())

        bin_checkm[bin_id_f3]={}
        try:
            bin_checkm[bin_id_f3]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f3]['Connections']=0
        bin_checkm[bin_id_f3]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
        bin_checkm[bin_id_f3]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
        try:
            bin_checkm[bin_id_f3]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f3]['Genome size']=0
        try:
            bin_checkm[bin_id_f3]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f3]['Mean scaffold length']=0
        bin_checkm[bin_id_f3]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())

        bin_checkm[bin_id_f4]={}
        try:
            bin_checkm[bin_id_f4]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f4]['Connections']=0
        bin_checkm[bin_id_f4]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
        bin_checkm[bin_id_f4]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
        try:
            bin_checkm[bin_id_f4]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f4]['Genome size']=0
        try:
            bin_checkm[bin_id_f4]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
        except:
            bin_checkm[bin_id_f4]['Mean scaffold length']=0
        bin_checkm[bin_id_f4]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
    return bin_checkm

def initial_drep_final_comparitor(final_iteration_folder, Coverage_maxtrix_list, datasets, num_threads, pwd, step):
    output_folder_name='BestBinset'
    try:
        os.system('mkdir BestBinset')
    except:
        os.system('mv -r BestBinset BestBinset_old')
        os.system('mkdir BestBinset')

    present_bins, bin_checkm={}, {}
    for root, dirs, files in os.walk(final_iteration_folder):
        os.chdir(pwd+'/'+final_iteration_folder)
        for file in files:
            if 'noclass' not in file and 'unbined' not in file:
                hz=file.split('.')[-1]
                if 'fasta' in hz or 'fa' in hz or 'fna' in hz:
                    present_bins[file]=0
            
                if 'bin_stats_ext.tsv' in file:
                    bin_checkm.update(parse_checkm(file))
    os.chdir(pwd)
    #print(str(bin_checkm))

    fx=open('Selected_bins_'+str(final_iteration_folder)+'.txt', 'w')
    f2x=open('Removed_bins_'+str(final_iteration_folder)+'.txt', 'w')
    fx.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    f2x.write('Removed_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')

    total_core_contigs={}
    try:
        for line in open(pwd+'/Core_contigs_total_bins.txt','r'):
            bin_id=str(line).strip().split('\t')[0]
            total_core_contigs[bin_id]={}
            bin_core_contigs_list=str(line).strip().split('\t')[1].replace('{','').replace('}','').replace('\'','').replace(': 0','').split(',')
            for cg in bin_core_contigs_list:
                total_core_contigs[bin_id][cg]=0
    except:
        for line in open(pwd+'/BestBinset_comparison_files/Core_contigs_total_bins.txt','r'):
            bin_id=str(line).strip().split('\t')[0]
            total_core_contigs[bin_id]={}
            bin_core_contigs_list=str(line).strip().split('\t')[1].replace('{','').replace('}','').replace('\'','').replace(': 0','').split(',')
            for cg in bin_core_contigs_list:
                total_core_contigs[bin_id][cg]=0
                
    # f=open('Total_coverage_matrix.txt','w')
    n, contigs_coverage=0, {}
    for item in Coverage_maxtrix_list:
        print('Processing '+str(item))
        n1=0
        for line in open(item, 'r'):
            n+=1
            n1+=1
            if n == 1:
                num=str(line).strip().count('drange')
                # f.write(str(line))
            if n1 >= 2:
                # f.write(str(line))
                ids=str(line).strip().split('\t')[0]
                contigs_coverage[str(ids)]={}
                for i in range(1, int(num)+1):
                    contigs_coverage[str(ids)][i]=float(str(line).strip().split('\t')[3*i+1])
    # f.close()

    print('Finding possible same contig')
    print('Parsing bins in '+str(final_iteration_folder))
    contigs, contig_bin_list, bin_total_length, bin_core_length, bin_contigs_num, bin_core_contigs_num, bin_contigs_coverage, bin_core_contigs_coverage = {}, {}, {}, {}, {}, {}, {}, {}
    contigs_len={}
    for root, dirs, files in os.walk(final_iteration_folder):
        os.chdir(pwd+'/'+final_iteration_folder)
        for file in files:
            if 'noclass' not in file and 'unbined' not in file:
                hz=file.split('.')[-1]

                if 'fasta' in hz or 'fa' in hz or 'fna' in hz:
                    bin_total_length[file], bin_core_length[file], bin_contigs_num[file], bin_core_contigs_num[file] = 0, 0, 0, 0
                    bin_contigs_coverage[file], bin_core_contigs_coverage[file] = {}, {}
                    for i in range(1, int(num)+1):
                        bin_contigs_coverage[file][i]=0
                        bin_core_contigs_coverage[file][i]=0

                    for record in SeqIO.parse(file,'fasta'):
                        contigs_len[record.id]=len(record.seq)
                        bin_total_length[file]+=len(record.seq)
                        contigs[str(record.id)]=str(record.seq)
                        try:
                            contig_bin_list[str(record.id)].append(str(file))
                        except:
                            contig_bin_list[str(record.id)]=[str(file)]
                        try:
                            for i in range(1, int(num)+1):
                                bin_contigs_coverage[file][i]+=contigs_coverage[str(record.id)][i]
                            bin_contigs_num[file]+=1
                        except:
                            xyzz=0

                        # if str(record.id) in total_core_contigs[file].keys():
                        #     bin_core_length[file]+=len(record.seq)
                        #     try:
                        #         for i in range(1, int(num)+1):
                        #             bin_core_contigs_coverage[file][i]+=contigs_coverage[str(record.id)][i]
                        #         bin_core_contigs_num[file]+=1
                        #     except:
                        #         xyzz=0
    os.chdir(pwd)

    bin_contigs_avg_coverage, bin_core_contigs_avg_coverage = {}, {}
    for bins in bin_contigs_coverage.keys():
        bin_contigs_avg_coverage[bins]={}
        if bin_contigs_num[bins] != 0:
            for i in range(1, int(num)+1):
                bin_contigs_avg_coverage[bins][i]=float(round(float(bin_contigs_coverage[bins][i]/bin_contigs_num[bins]),2))

    # for bins in bin_core_contigs_coverage.keys():
    #     bin_core_contigs_coverage[bins]={}
    #     if bin_core_contigs_num[bins] != 0:
    #         for i in range(1, int(num)+1):
    #             bin_core_contigs_avg_coverage[bins][i]=float(round(float(bin_core_contigs_coverage[bins][i]/bin_core_contigs_num[bins]),2))

    t_contigs=str(final_iteration_folder)+'_contigs.fa'
    total_contigs_file=open(t_contigs,'w')
    for ids in contigs.keys():
        total_contigs_file.write('>'+str(ids)+'\n'+str(contigs[ids])+'\n')
    total_contigs_file.close()

    blast_output=Contigs_aligner(t_contigs, num_threads)
    os.system('rm '+t_contigs)

    bin_similar, bin_core_similar, bin_similar_nor_sum, bin_similar_nor_item_num = {}, {}, {}, {}
    for line in open(blast_output, 'r'):
        contig1=str(line).strip().split('\t')[0]
        contig2=str(line).strip().split('\t')[1]
        simi=str(line).strip().split('\t')[2]
        leng=str(line).strip().split('\t')[3]
        aligned_length=float(simi)*int(leng)/100
        bin1s=contig_bin_list[str(contig1)]
        bin2s=contig_bin_list[str(contig2)]

        sim_bins = {}
        for item in bin1s:
            for item2 in bin2s:
                bins=[item, item2]
                bins.sort()
                sim_bins[bins[0]+'\t'+bins[1]]=aligned_length

        for item in sim_bins.keys():
            bin1, bin2 = item.split('\t')[0], item.split('\t')[1]
            if item not in bin_core_similar.keys():
                bin_core_similar[item]={}
                bin_similar_nor_sum[item]={}
                bin_similar_nor_item_num[item]={}
                for i in range(1, int(num)+1):
                    bin_similar_nor_sum[item][i]=0
                    bin_similar_nor_item_num[item][i]=0
                
            try:
                bin_similar[item]+=aligned_length
            except:
                bin_similar[item]=aligned_length

            try:
                bin_core_similar[item][bin1]+=aligned_length
            except:
                bin_core_similar[item][bin1]=aligned_length

            try:
                bin_core_similar[item][bin2]+=aligned_length
            except:
                bin_core_similar[item][bin2]=aligned_length

            if int(leng) >= 1000:
                contig1_ratio=100*int(leng)/int(contigs_len[contig1])
                contig2_ratio=100*int(leng)/int(contigs_len[contig2])
                if contig1_ratio >= 50 or contig2_ratio >= 50:
                    for i in range(1, int(num)+1):
                        if float(contigs_coverage[contig2][i]) != 0:
                            ratio=contigs_coverage[contig1][i]/contigs_coverage[contig2][i]
                            bin_similar_nor_sum[item][i]+=ratio
                            bin_similar_nor_item_num[item][i]+=1
    
    bin_similar_nor_avg={}
    for item in bin_similar_nor_sum.keys():
        bin_similar_nor_avg[item]={}
        for i in range(1, int(num)+1):
            try:
                bin_similar_nor_avg[item][i]=bin_similar_nor_sum[item][i]/bin_similar_nor_item_num[item][i]
            except:
                bin_similar_nor_avg[item][i]=0

    f=open('Similar_bin_in_'+str(final_iteration_folder)+'.txt','w')
    f2=open('Highly_possible_similar_bin_in_'+str(final_iteration_folder)+'.txt','w')
    f.write('Bin1'+'\t'+'Bin2'+'\t'+'Aligned(bp)'+'\t'+'Percetage_bin1(%)'+'\t'+'Percetage_bin2(%)'+'\t'+'Core_bin1_aligned(bp)'+'\t'+'Core_bin1_len(bp)'+'\t'+'Core_bin2_aligned(bp)'+'\t'+'Core_bin2_len(bp)'+'\t'+'Percetage_bin1_core(%)'+'\t'+'Percetage_bin2_core(%)'+'\n')
    f2.write('Bin1'+'\t'+'Bin2'+'\t'+'Aligned(bp)'+'\t'+'Percetage_bin1(%)'+'\t'+'Percetage_bin2(%)'+'\t'+'Core_bin1_aligned(bp)'+'\t'+'Core_bin1_len(bp)'+'\t'+'Core_bin2_aligned(bp)'+'\t'+'Core_bin2_len(bp)'+'\t'+'Percetage_bin1_core(%)'+'\t'+'Percetage_bin2_core(%)'+'\n')
    bin_sim_percentage, filtrated_bin, total_bins_4_mapping, mapping_group, confirmed_bin = {}, {}, {}, {}, {}
    for item in bin_similar.keys():
        bin1=item.split('\t')[0]
        bin2=item.split('\t')[1]
        if bin1 != bin2:
            length=float(bin_similar[item])
            sim1=float(round(50*length/int(bin_total_length[bin1]), 2)) ### 50 means
            sim2=float(round(50*length/int(bin_total_length[bin2]), 2)) ### 50 means
            bin1_core_aligned_len=bin_core_similar[item][bin1]
            bin2_core_aligned_len=bin_core_similar[item][bin2]
            temp1=float(round(bin_core_length[bin1], 2))
            if temp1 != 0:
                sim_core1=float(round(100*float(bin1_core_aligned_len)/temp1, 2))
            else:
                sim_core1=0
            temp2=float(round(bin_core_length[bin2], 2))
            if temp2 != 0:
                sim_core2=float(round(100*float(bin2_core_aligned_len)/temp2, 2))
            else:
                sim_core2=0
            bin_sim_percentage[item]=str(sim1)+'\t'+str(sim2)
            f.write(item+'\t'+str(bin_similar[item])+'\t'+str(sim1)+'\t'+str(sim2)+'\t'+str(bin1_core_aligned_len)+'\t'+str(bin_core_length[bin1])+'\t'+str(bin2_core_aligned_len)+'\t'+str(bin_core_length[bin2])+'\t'+str(sim_core1)+'\t'+str(sim_core2)+'\n')
            sum_sim=sim1+sim2
            sum_sim_core=sim_core1+sim_core2
            if sum_sim >= 100:
                # if sim1 >= 65 or sim2 >= 65:
                f2.write(item+'\t'+str(bin_similar[item])+'\t'+str(sim1)+'\t'+str(sim2)+'\t'+str(bin1_core_aligned_len)+'\t'+str(bin_core_length[bin1])+'\t'+str(bin2_core_aligned_len)+'\t'+str(bin_core_length[bin2])+'\t'+str(sim_core1)+'\t'+str(sim_core2)+'\n')
                filtrated_bin[item]=bin_similar[item]
                total_bins_4_mapping[bin1]=0
                total_bins_4_mapping[bin2]=0
                if sum_sim >= 180 or sum_sim_core >= 180 or sim1 >= 95 or sim2 >= 95 or sim_core1 >= 95 or sim_core2 >= 95:
                    confirmed_bin[item]=bin_similar[item]
                    del total_bins_4_mapping[bin1]
                    del total_bins_4_mapping[bin2]
    f.close()
    f2.close()

    ### Coverage filtration
    mapping_group['1'], mapping_group['2'], mapping_group['total'] = {}, {}, {}
    print(str(confirmed_bin))
    filtrated_bin2=copy.deepcopy(filtrated_bin)
    try:
        f_p_rep=open('Potential_replicate_but_cov_inconsistence_bins.txt','a')
    except:
        f_p_rep=open('Potential_replicate_but_cov_inconsistence_bins.txt','w')
    for item in filtrated_bin2.keys():
        if item not in confirmed_bin.keys():
            bin1=item.split('\t')[0]
            bin2=item.split('\t')[1]
            x=0
            if item in bin_similar_nor_avg.keys():
                for i in range(1, int(num)+1):
                    try:
                        bin2_normalized_cov=bin_contigs_avg_coverage[bin2][i]*bin_similar_nor_avg[item][i]
                    except:
                        bin2_normalized_cov=bin_contigs_avg_coverage[bin2][i]

                    delta=abs(bin_contigs_avg_coverage[bin1][i]-bin2_normalized_cov)
                    sum=bin_contigs_avg_coverage[bin1][i]+bin2_normalized_cov
                    if sum != 0:
                        qua=100*delta/sum
                    else:
                        qua=1000*delta
                        
                    if qua <= 10:
                        x+=1
        
            if x == num:
                try:
                    del total_bins_4_mapping[bin1]
                    del total_bins_4_mapping[bin2]
                except:
                    yzh=0
            else:
                f_p_rep.write(str(bin1)+'\t'+str(bin2)+'\n')
                # if step == 'final_drep':
                #     if bin1 not in mapping_group['total'].keys():
                #         mapping_group['1'][bin1]=0
                #         mapping_group['total'][bin1]=1
                #     if bin2 not in mapping_group['total'].keys():
                #         mapping_group['2'][bin2]=0
                #         mapping_group['total'][bin2]=2
                # else:
                del filtrated_bin[item]
    f_p_rep.close()

    bin_checkm, bin_checkm_o, qua_file = {}, {}, {}
    os.chdir(pwd+'/'+final_iteration_folder)
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                qua_file[file]=0
                for line in open(file, 'r'):
                    bin_id=str(line).strip().split('\t')[0]
                    bin_id_f=bin_id+'.fa'
                    bin_id_f2=bin_id+'.fasta'
                    bin_checkm_o[bin_id]=str(line).strip()
                    bin_checkm[bin_id_f]={}
                    bin_checkm[bin_id_f]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
                    bin_checkm[bin_id_f]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]={}
                    bin_checkm[bin_id_f2]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
                    bin_checkm[bin_id_f2]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())

    os.chdir(pwd)
    
    os.system('mv Similar_bin_in_final_iteration.txt Highly_possible_similar_bin_in_final_iteration.txt '+str(output_folder_name))
    os.chdir(str(output_folder_name))

    selected_bins, eliminated_bins={}, {}
    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}

    fy=open('Bins_wth_sameQua_difSize.txt','a')

    del_filtration={}
    for item in filtrated_bin.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]
        bin_score_1=marker_score[bin_checkm[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[bin_checkm[ID2]['marker lineage'].split('__')[0]]

        delta_1=bin_checkm[ID1]['Completeness']-bin_checkm[ID1]['Contamination']
        delta_2=bin_checkm[ID2]['Completeness']-bin_checkm[ID2]['Contamination']

        if bin_checkm[ID1]['Contamination'] >= 12 and float(delta_1) <= 80:
            del_filtration[item]=0
            eliminated_bins[ID1]=bin_checkm[ID1]
        
        if bin_checkm[ID2]['Contamination'] >= 12 and float(delta_2) <= 80:
            del_filtration[item]=0
            eliminated_bins[ID2]=bin_checkm[ID2]

    for item in del_filtration.keys():
        del filtrated_bin[item]

    for item in filtrated_bin.keys():
        ID1=str(item).split('\t')[0]
        ID2=str(item).split('\t')[1]
        bin_score_1=marker_score[bin_checkm[ID1]['marker lineage'].split('__')[0]]
        bin_score_2=marker_score[bin_checkm[ID2]['marker lineage'].split('__')[0]]

        if bin_score_1 == bin_score_2:
            CPN_CTN_1=bin_checkm[ID1]['Completeness']-5*bin_checkm[ID1]['Contamination']
            CPN_CTN_2=bin_checkm[ID2]['Completeness']-5*bin_checkm[ID2]['Contamination']
            gz1=bin_checkm[ID1]['Genome size']
            gz2=bin_checkm[ID2]['Genome size']
            msl1=bin_checkm[ID1]['Mean scaffold length']
            msl2=bin_checkm[ID2]['Mean scaffold length']

            if CPN_CTN_1 == CPN_CTN_2:
                if float(gz1) == float(gz2) and float(msl1) == float(msl2):
                    selected_bins[ID1]=bin_checkm[ID1]
                    eliminated_bins[ID2]=bin_checkm[ID2]
                    fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                else:
                    fy.write(ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            elif CPN_CTN_1 > CPN_CTN_2:
                selected_bins[ID1]=bin_checkm[ID1]
                eliminated_bins[ID2]=bin_checkm[ID2]
                fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            else:
                selected_bins[ID2]=bin_checkm[ID2]
                eliminated_bins[ID1]=bin_checkm[ID1]
                fx.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                f2x.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
        elif bin_score_1 > bin_score_2:
            selected_bins[ID1]=bin_checkm[ID1]
            eliminated_bins[ID2]=bin_checkm[ID2]
            fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
        else:
            selected_bins[ID2]=bin_checkm[ID2]   
            eliminated_bins[ID1]=bin_checkm[ID1]
            fx.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            f2x.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
    fx.close()
    f2x.close()
    fy.close()

    for item in eliminated_bins.keys():
       if item in selected_bins.keys():
           del selected_bins[item]

    eliminated_bins_checkm={}
    for item in eliminated_bins.keys():
        name_list=item.split('.')
        name_list.remove(name_list[-1])
        item_name='.'.join(name_list)
        eliminated_bins_checkm[item_name]=eliminated_bins[item]

    for item in qua_file.keys():
        file_name=item.split('.tsv')[0]+'_o.tsv'
        os.system('mv '+item+' '+str(file_name))

    f=open('Best_binset_bin_stats_ext.tsv', 'w')
    f2=open('Medium_quality_bins.txt','w')
    f3=open('Highly_quality_bins.txt','w')
    for item in bin_checkm_o.keys():
        if item not in eliminated_bins_checkm.keys():
            f.write(str(bin_checkm_o[item])+'\n')
            Completeness=str(bin_checkm_o[item]).split('Completeness\': ')[1].split(',')[0]
            Contamination=str(bin_checkm_o[item]).split('Contamination\': ')[1].split('}')[0]
            delta_value=float(Completeness)-5*float(Contamination)
            if delta_value >= 70:
                f3.write(str(bin_checkm_o[item])+'\n')
            elif delta_value >= 50:
                f2.write(str(bin_checkm_o[item])+'\n')
            else:
                continue
    f.close()
    f2.close()
    f3.close()

    f=open('Genome_group_all_list_'+str(output_folder_name)+'.txt','w')
    print('Parsing files in '+str(output_folder_name))
    os.chdir(pwd+'/'+final_iteration_folder)
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fa' in hz:
                    if file not in eliminated_bins.keys():
                        os.system('cp '+file+' '+pwd+'/'+str(output_folder_name))
            else:
                hz=str(file).split('.')[-1]
                if hz == 'fasta' or hz == 'fa':
                    if file not in eliminated_bins.keys():
                        os.system('cp '+file+' '+pwd+'/'+str(output_folder_name))
    
    os.chdir(pwd+'/'+str(output_folder_name))
    bin_folders, title={}, []
    for root, dirs, files in os.walk(pwd+'/'+str(output_folder_name)):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    bin_source_folder=str(file).split('_genomes.')[0]
                    bin_folders[bin_source_folder]=1


    for item in bin_folders.keys():
        try:
            os.chdir(pwd+'/'+item+'_genomes')
            for root, dirs, files in os.walk(pwd+'/'+item+'_genomes'):
                for file in files:
                    if 'Genome_group_all_list_' in file:
                        n=0
                        for line in open(file, 'r'): 
                            n+=1
                            if n == 1 and len(title) == 0:
                                f.write(str(line))
                                title.append(str(line))
                            else:
                                if str(line).strip().split('\t')[0] not in eliminated_bins.keys():
                                    f.write(str(line))
        except:
            print(str(item)+'_genomes no existed')
    f.close()

    os.chdir(pwd)
    os.system('mkdir '+str(output_folder_name)+'_comparison_files')
    os.system('rm temp.orfs.* Contigs_iteration_* *.nhi *.nhr *.nin *.nog *.nsd *.nsi *.nsq *.nhd temp_db.txt *_db.txt')
    os.system('rm *.bt2  *.bt2l')
    os.system('mv bins_depth_comparison.txt Contigs_sharing_bins.txt Contigs_iteration_* Coverage_matrix_for_binning_iteration_* Selected_bins_* Removed_bins_* Extract_bins_* Bins_similar_* Contigs_scoring_* Raw_bins_* Filtrated_* *_vs_* *_gc.txt *_bins_coverage.txt Test_Raw_Bins_Comparison_* Total_bins_similar_coverage_* Core_contigs_* core_contigs_* Single_mapping_depth.txt Abnor_normalization_record.txt Potential_similar_* Unfiltrated_potential_similar_* Potential_replicate_but_cov_inconsistence_bins.txt '+str(output_folder_name)+'_comparison_files')
    os.system('rm -rf CC_* Merge_binset*')

def final_binset_comparitor(final_iteration_folder, Coverage_maxtrix_list, datasets, num_threads, pwd, step):
    output_folder_name=final_iteration_folder

    present_bins, present_bins_org, bin_checkm, bin_status = {}, {}, {}, {}
    os.chdir(pwd+'/'+final_iteration_folder)
    #print(str(final_iteration_folder))
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                bin_checkm.update(parse_checkm(file))

    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if 'noclass' not in file and 'unbined' not in file:
                hz=file.split('.')[-1]
                if 'fasta' in hz or 'fa' in hz or 'fna' in hz:
                    present_bins[file]=0
                    if step == 'final_drep':
                        try:
                            bin_org_name=str(file).split('.')[0].split('_')[0]
                            present_bins_org[bin_org_name]=file
                        except:
                            bin_org_name=str(file).split('.')[0]
                            present_bins_org[bin_org_name]=file
                    else:
                        bin_org_name_list=str(file).split('.')
                        bin_org_name_list.remove(bin_org_name_list[-1])
                        bin_org_name='.'.join(bin_org_name_list)
                        present_bins_org[bin_org_name]=file

                    seq_num, seq_bp = 0, 0
                    for record in SeqIO.parse(file,'fasta'):
                        seq_num+=1
                        seq_bp+=len(record.seq)
                    bin_status[file]=seq_bp/seq_num
                    bin_checkm[file]['Genome size']=int(seq_bp)
                    bin_checkm[file]['Mean scaffold length']=float(seq_bp/seq_num)
    os.chdir(pwd)

    if step == 'final_drep':
        mod_bin1, mod_bin2, org_mod_bins = {}, {}, {}
        mod_folder=str(final_iteration_folder).split('_re-assembly')[0]+'_mod'
        for line in open(pwd+'/'+str(mod_folder)+'/Bin_name_mod.txt','r'):
            org_bin=str(line).strip().split('\t')[0]
            org_bin_name_list=org_bin.split('_genomes.')[1].split('.')
            # org_bin_name_list.remove(org_bin_name_list[-1])
            #org_bin_checkm_name='.'.join(org_bin_name_list)
            org_bin2=org_bin.split('_genomes.')[0]+'_genomes.'+org_bin_name_list[0]
            mod_bin=str(line).strip().split('\t')[1].strip()
            mod_bin_name2=str(mod_bin).split('.fa')[0]
            #mod_bin_checkm=str(line).strip().split('\t')[1].strip().split('.fa')[0]
            bin_checkm[org_bin]=bin_checkm[mod_bin]
            bin_checkm[org_bin+'sta']=bin_checkm[mod_bin]
            mod_bin1[org_bin]=mod_bin
            mod_bin2[org_bin+'sta']=mod_bin
            org_mod_bins[org_bin2]=mod_bin_name2

    mapping_group, filtrated_bin = {}, {}
    fx=open('Selected_bins_'+str(final_iteration_folder)+'.txt', 'w')
    f2x=open('Removed_bins_'+str(final_iteration_folder)+'.txt', 'w')
    fx.write('Selected_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    f2x.write('Removed_bin'+'\t'+'Similar_bins'+'\t'+'Aligned_length(bp)'+'\t'+'Aligned_A_percentage(%)'+'\t'+'Aligned_B_percentage(%)'+'\t'+'Average_Coverage_Variation'+'\t'+'Bin1_factors'+'\t'+'Bin2_factors'+'\n')
    try:
        print('Comparing priously found redundant bins')
        fy=open('Bins_wth_sameQua_difSize2.txt','w')
        for line in open('Bins_wth_sameQua_difSize.txt','r'):
            bin1=str(line).strip().split('\t')[0].split('---')[0]
            bin2=str(line).strip().split('\t')[0].split('---')[1]
            print('Comparing '+str(bin1)+' with '+str(bin2))
            if bin1 in present_bins.keys() and bin2 in present_bins.keys():
                bin1_completeness=bin_checkm[bin1]['Completeness']
                bin2_completeness=bin_checkm[bin2]['Completeness']
                bin1_contamination=bin_checkm[bin1]['Contamination']
                bin2_contamination=bin_checkm[bin2]['Contamination']
                bin1_qua=bin1_completeness-5*bin1_contamination
                bin2_qua=bin2_completeness-5*bin2_contamination

                if bin1_qua > bin2_qua:
                    if step != 'final_drep':
                        os.system('rm '+str(bin2))
                        print('Remove bin: '+str(bin2)+'. It is redundant bin with '+str(bin1))
                        fx.write(str(bin1)+'\t'+str(line))
                        f2x.write(str(bin2)+'\t'+str(line))
                    else:
                        os.system('rm '+str(mod_bin1[bin2])+' '+str(mod_bin2[bin2]))
                        print('Remove bin: '+str(bin2)+'. It is redundant bin with '+str(bin1))
                        fx.write(str(bin1)+'/'+str(mod_bin1[bin1])+'\t'+str(line))
                        f2x.write(str(bin2)+'/'+str(mod_bin1[bin2])+'\t'+str(line))
                elif bin2_qua > bin1_qua:
                    if step != 'final_drep':
                        os.system('rm '+str(bin1))
                        print('Remove bin: '+str(bin1)+'. It is redundant bin with '+str(bin2))
                        fx.write(str(bin2)+'\t'+str(line))
                        f2x.write(str(bin1)+'\t'+str(line))
                    else:
                        os.system('rm '+str(mod_bin1[bin1])+' '+str(mod_bin2[bin1]))
                        print('Remove bin: '+str(bin1)+'. It is redundant bin with '+str(bin2))
                        fx.write(str(bin2)+'/'+str(mod_bin1[bin2])+'\t'+str(line))
                        f2x.write(str(bin1)+'/'+str(mod_bin1[bin1])+'\t'+str(line))
                else:
                    fy.write(line)
                    filtrated_bin[bin1+'\t'+bin2]=0
        fy.close()
        os.system('mv Bins_wth_sameQua_difSize.txt Bins_wth_sameQua_difSize_'+str(final_iteration_folder)+'.txt')
        os.system('mv Bins_wth_sameQua_difSize2.txt Bins_wth_sameQua_difSize.txt')
    except:
        print('Bins_wth_sameQua_difSize.txt dose not exits')

    total_bins_4_mapping, mapping_group['1'], mapping_group['2'] = {}, {}, {}
    try:
        for line in open('Potential_replicate_but_cov_inconsistence_bins.txt','r'):
            bin1=str(line).strip().split('\t')[0]
            bin1_name_list=bin1.split('.')
            bin1_name_list.remove(bin1_name_list[-1])
            bin1_name='.'.join(bin1_name_list)
            bin2=str(line).strip().split('\t')[1]
            bin2_name_list=bin2.split('.')
            bin2_name_list.remove(bin2_name_list[-1])
            bin2_name='.'.join(bin2_name_list)
            axyzzz=0
            if step == 'final_drep':
                try:
                    mod_bin1_id=org_mod_bins[bin1_name]
                    axyzzz+=1
                except:
                    axyzzz+=0
            else:
                mod_bin1_id=bin1_name
                if bin1 in present_bins.keys():
                    axyzzz+=1
            
            if step == 'final_drep':
                try:
                    mod_bin2_id=org_mod_bins[bin2_name]
                    axyzzz+=1
                except:
                    axyzzz+=0
            else:
                mod_bin2_id=bin2_name
                if bin2 in present_bins.keys():
                    axyzzz+=1
                
            if axyzzz == 2:
                total_bins_4_mapping[mod_bin1_id]=1
                total_bins_4_mapping[mod_bin2_id]=1
                mapping_group['1'][mod_bin1_id]=1
                mapping_group['2'][mod_bin2_id]=1
                filtrated_bin[mod_bin1_id+'\t'+mod_bin2_id]=0
    except:
        print('There is no bin with potential_replicate_but_cov_inconsistence_bins')

    if step == 'final_drep':
        if len(total_bins_4_mapping) != 0:
            print(str(len(total_bins_4_mapping))+' bins are waiting for re-mapping')
            print(str(total_bins_4_mapping))
            f_g1=open('Merged_group1_contigs.fa','w')
            f_g2=open('Merged_group2_contigs.fa','w')
            os.chdir(pwd+'/'+final_iteration_folder)
            for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
                for file in files:
                    try:
                        bin_name=file.split('.')[0].split('_')[0]
                    except:
                        bin_name=file.split('.')[0]
                    if bin_name in mapping_group['1'].keys():
                        for record in SeqIO.parse(file, 'fasta'):
                            f_g1.write('>'+str(file)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
                    elif bin_name in mapping_group['2'].keys():
                        for record in SeqIO.parse(file, 'fasta'):
                            f_g2.write('>'+str(file)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
            os.chdir(pwd)
            f_g1.close()
            f_g2.close()

            bins_mo_depth = {}
            bins_mo_depth.update(final_iteration_mapping('Merged_group1_contigs.fa', datasets, num_threads, pwd))
            bins_mo_depth.update(final_iteration_mapping('Merged_group2_contigs.fa', datasets, num_threads, pwd))
    
            fmd=open('Single_mapping_depth.txt','w')
            for bins in bins_mo_depth.keys():
                fmd.write(str(bins)+'\t'+str(bins_mo_depth[bins])+'\n')
            fmd.close()

            fxx=open('Paired_bins_and_filtrateed_bins.txt','w')
            bin_num_org, xxx=len(filtrated_bin), 0
            filtrated_bin2=copy.deepcopy(filtrated_bin)
            for item in filtrated_bin2.keys():
                bin1=item.split('\t')[0]
                bin2=item.split('\t')[1]
                fxx.write(str(item)+'\n')
                try:
                    for i in range(1, len(bins_mo_depth[bin1])+1):
                        avg_cov1=bins_mo_depth[bin1][i]
                        avg_cov2=bins_mo_depth[bin2][i]
                        delta=abs(avg_cov1-avg_cov2)
                        avg=(avg_cov1+avg_cov2)/2
                        perc=100*delta/avg
                        if perc > 20:
                            del filtrated_bin[item]
                            fxx.write('\t'+'Not redundant'+'\n')
                except:
                    xxx+=1
            fxx.close()

            bin_num_after=len(filtrated_bin)
            delta_bin=bin_num_org-bin_num_after
            print('Removed '+str(delta_bin)+' potential redundant bin(s)')

    bin_checkm, bin_checkm_o, qua_file = {}, {}, {}
    os.chdir(pwd+'/'+final_iteration_folder)
    for root, dirs, files in os.walk(pwd+'/'+final_iteration_folder):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                qua_file[file]=0
                for line in open(file, 'r'):
                    bin_id=str(line).strip().split('\t')[0]
                    bin_id_f=bin_id+'.fa'
                    bin_id_f2=bin_id+'.fasta'
                    bin_checkm_o[bin_id]=str(line).strip()
                    bin_checkm[bin_id_f]={}
                    #bin_checkm[bin_id_f]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
                    bin_checkm[bin_id_f]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
                    try:
                        bin_checkm[bin_id_f]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip().replace('\'','').replace('\"',''))
                    except:
                        bin_checkm[bin_id_f]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]={}
                    #bin_checkm[bin_id_f2]['Connections']=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['marker lineage']=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\',')[0].strip()
                    bin_checkm[bin_id_f2]['Completeness']=float(str(line).strip().split('Completeness\': ')[1].split('}')[0].split(',')[0].strip())
                    try:
                        bin_checkm[bin_id_f2]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip().replace('\'','').replace('\"',''))
                    except:
                        bin_checkm[bin_id_f2]['Genome size']=float(str(line).strip().split('Genome size\': ')[1].split('}')[0].split(',')[0].strip())
                    bin_checkm[bin_id_f2]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
                    try:
                        bin_checkm[bin_id_f]['Mean scaffold length']=bin_status[bin_id_f]
                    except:
                        bin_checkm[bin_id_f2]['Mean scaffold length']=bin_status[bin_id_f2]
    os.chdir(pwd)

    #os.system('mv Similar_bin_in_final_iteration.txt Highly_possible_similar_bin_in_final_iteration.txt '+str(output_folder_name))
    os.chdir(str(output_folder_name))

    selected_bins, eliminated_bins={}, {}
    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}

###
    if step != 'final_drep':
        fy=open('Bins_wth_sameQua_difSize.txt','a')
    else:
        fxy=open('Final_potential_replicated_bins_wth_same_quality_value.txt','w')
        replicated_bins={}
###

    del_filtration={}
    for item in filtrated_bin.keys():
        # if step == 'final_drep':
        try:
            ID1=present_bins_org[str(item).split('\t')[0]]
            ID2=present_bins_org[str(item).split('\t')[1]]
        # else:
        #     ID1=str(item).split('\t')[0]
        #     ID2=str(item).split('\t')[1]

            bin_score_1=marker_score[bin_checkm[ID1]['marker lineage'].split('__')[0]]
            bin_score_2=marker_score[bin_checkm[ID2]['marker lineage'].split('__')[0]]

            delta_1=bin_checkm[ID1]['Completeness']-bin_checkm[ID1]['Contamination']
            delta_2=bin_checkm[ID2]['Completeness']-bin_checkm[ID2]['Contamination']

            if bin_checkm[ID1]['Contamination'] >= 12 and float(delta_1) <= 80:
                del_filtration[item]=0
                eliminated_bins[ID1]=bin_checkm[ID1]
            
            if bin_checkm[ID2]['Contamination'] >= 12 and float(delta_2) <= 80:
                del_filtration[item]=0
                eliminated_bins[ID2]=bin_checkm[ID2]
        except:
            xyxyxy=0

    for item in del_filtration.keys():
        del filtrated_bin[item]

    for item in filtrated_bin.keys():
        try:
            if step == 'final_drep':
                ID1=present_bins_org[str(item).split('\t')[0]]
                ID2=present_bins_org[str(item).split('\t')[1]]
            else:
                ID1=str(item).split('\t')[0]
                ID2=str(item).split('\t')[1]
            bin_score_1=marker_score[bin_checkm[ID1]['marker lineage'].split('__')[0]]
            bin_score_2=marker_score[bin_checkm[ID2]['marker lineage'].split('__')[0]]

            if bin_score_1 == bin_score_2:
                CPN_CTN_1=bin_checkm[ID1]['Completeness']-5*bin_checkm[ID1]['Contamination']
                CPN_CTN_2=bin_checkm[ID2]['Completeness']-5*bin_checkm[ID2]['Contamination']
                gz1=bin_checkm[ID1]['Genome size']
                gz2=bin_checkm[ID2]['Genome size']
                msl1=bin_checkm[ID1]['Mean scaffold length']
                msl2=bin_checkm[ID2]['Mean scaffold length']

                if CPN_CTN_1 == CPN_CTN_2:
                    if float(gz1) == float(gz2) and float(msl1) == float(msl2):
                        selected_bins[ID1]=bin_checkm[ID1]
                        eliminated_bins[ID2]=bin_checkm[ID2]
                        fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                        f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    else:
                        if step != 'final_drep':
                            fy.write(ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                        else:
                            fxy.write(ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                            replicated_bins[ID1]=ID2
                    # bin1=bin_checkm[ID1]['Connections']/bin_checkm[ID1]['Genome size']
                    # bin2=bin_checkm[ID2]['Connections']/bin_checkm[ID2]['Genome size']
                    # if bin1 <= bin2:
                    #     selected_bins[ID1]=bin_checkm[ID1]
                    #     eliminated_bins[ID2]=bin_checkm[ID2]
                    #     fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    #     f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    # else:
                    #     selected_bins[ID2]=bin_checkm[ID2]
                    #     eliminated_bins[ID1]=bin_checkm[ID1]
                    #     fx.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    #     f2x.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                elif CPN_CTN_1 > CPN_CTN_2:
                    selected_bins[ID1]=bin_checkm[ID1]
                    eliminated_bins[ID2]=bin_checkm[ID2]
                    fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                else:
                    selected_bins[ID2]=bin_checkm[ID2]
                    eliminated_bins[ID1]=bin_checkm[ID1]
                    fx.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                    f2x.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            elif bin_score_1 > bin_score_2:
                selected_bins[ID1]=bin_checkm[ID1]
                eliminated_bins[ID2]=bin_checkm[ID2]
                fx.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                f2x.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
            else:
                selected_bins[ID2]=bin_checkm[ID2]   
                eliminated_bins[ID1]=bin_checkm[ID1]
                fx.write(ID2+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
                f2x.write(ID1+'\t'+ID1+'---'+ID2+'\t'+str(filtrated_bin[item])+'\t'+str(bin_checkm[ID1])+'\t'+str(bin_checkm[ID2])+'\n')
        except:
            xyzzz=0
            
    fx.close()
    f2x.close()
    fy.close()
    #fxy.close()

    for item in eliminated_bins.keys():
       if item in selected_bins.keys():
           del selected_bins[item]

    eliminated_bins_checkm={}
    for item in eliminated_bins.keys():
        name_list=item.split('.')
        name_list.remove(name_list[-1])
        item_name='.'.join(name_list)
        eliminated_bins_checkm[item_name]=eliminated_bins[item]
        os.system('rm '+str(item))

    for item in qua_file.keys():
        file_name=item.split('.tsv')[0]+'_o.tsv'
        os.system('mv '+item+' '+str(file_name))

    f=open('Best_binset_bin_stats_ext.tsv', 'w')
    f2=open('Medium_quality_bins.txt','w')
    f3=open('Highly_quality_bins.txt','w')
    for item in bin_checkm_o.keys():
        if item not in eliminated_bins_checkm.keys():
            f.write(str(bin_checkm_o[item])+'\n')
            Completeness=str(bin_checkm_o[item]).split('Completeness\': ')[1].split(',')[0]
            Contamination=str(bin_checkm_o[item]).split('Contamination\': ')[1].split('}')[0].split(',')[0]
            delta_value=float(Completeness)-5*float(Contamination)
            if delta_value >= 70:
                f3.write(str(bin_checkm_o[item])+'\n')
            elif delta_value >= 50:
                f2.write(str(bin_checkm_o[item])+'\n')
            else:
                continue
    f.close()
    f2.close()
    f3.close()

    os.chdir(pwd)
    os.system('mkdir '+str(output_folder_name)+'_comparison_files')
    # os.system('rm temp.orfs.* Contigs_iteration_* *.nhi *.nhr *.nin *.nog *.nsd *.nsi *.nsq *.nhd temp_db.txt *_db.txt')
    os.system('rm *.bt2  *.bt2l')
    os.system('mv bins_depth_comparison.txt Contigs_sharing_bins.txt Contigs_iteration_* Coverage_matrix_for_binning_iteration_* Selected_bins_* Removed_bins_* Extract_bins_* Bins_similar_* Contigs_scoring_* Raw_bins_* Filtrated_* *_vs_* *_gc.txt *_bins_coverage.txt Test_Raw_Bins_Comparison_* Total_bins_similar_coverage_* Single_mapping_depth.txt Abnor_normalization_record.txt Potential_similar_* Unfiltrated_potential_similar_* Potential_replicate_but_cov_inconsistence_bins.txt '+str(output_folder_name)+'_comparison_files')
    # os.system('rm -rf CC_* Merge_binset*')

def record_bin_coverage(best_binset_from_multi_assemblies, coverage_file):
    pwd=os.getcwd()
    bin_contigs, bin_contigs_mock, total_bin_contigs, assembly_list, m={}, {}, {}, {}, 0
    print('Parsing bins '+str(best_binset_from_multi_assemblies))
    for root, dirs, files in os.walk(pwd+'/'+str(best_binset_from_multi_assemblies)):
        os.chdir(pwd+'/'+str(best_binset_from_multi_assemblies))
        for file in files:
            # try:
            #     hz=file.split('_genomes.')[-1]
            #     qz=file.split('_genomes.')[0]
            #     assembly_name_list=qz.split('_')
            #     assembly_name_list.remove(assembly_name_list[-1])
            #     assembly_name_list.remove(assembly_name_list[-1])
            #     assembly_name='_'.join(assembly_name_list)
            #     assembly_list[assembly_name]=1
            # except:
            #     hz=file.split('.')[-1]

            hz=file.split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                m+=1
                bin_contigs[file]={}
                bin_contigs_mock[file]={}
                for record in SeqIO.parse(file, 'fasta'):
                    total_bin_contigs[str(record.id)]=str(record.seq)
                    bin_contigs[file][str(record.id)]=str(record.seq)
                    bin_contigs_mock[file][str(record.id)]=0
    print('Parsed '+str(m)+' bins')

    os.chdir(pwd)
    
    # f=open(best_binset_from_multi_assemblies+'_refined_total_bins_contigs.fa', 'w')
    # for item in total_bin_contigs.keys():
    #     f.write('>'+str(item)+'\n'+str(total_bin_contigs[item])+'\n')
    # f.close()

    print('Recording the coverage of contigs from bins')
    n, contig_cov=0, {}
    for line in open(coverage_file,'r'):
        n+=1
        if n >= 2:
            ls=str(line).strip().split('\t')
            num=int((len(ls)-4)/3)
            ids=str(line).strip().split('\t')[0]
            contig_cov[ids]={}
            for i in range(1,num+1):
                contig_cov[ids][i]=float(str(line).strip().split('\t')[3*i+1])
        
    bin_contig_cov={}    
    for contig_id in contig_cov.keys():
        for bin_name in bin_contigs.keys():
            if bin_name not in bin_contig_cov.keys():
                bin_contig_cov[bin_name]={} 
            bin_contigs_mock[bin_name][contig_id]=1
            # if contig_id in bin_contigs[bin_name].keys():
            if len(bin_contigs_mock[bin_name]) == len(bin_contigs[bin_name]):
                bin_contig_cov[bin_name][contig_id]={} 
                for i in range(1,num+1):
                    bin_contig_cov[bin_name][contig_id][i]=contig_cov[contig_id][i]
            else:
                del bin_contigs_mock[bin_name][contig_id]
        
    print('Writing bin coverage matrix file')
    os.system('mkdir bin_coverage')
    os.chdir('bin_coverage')
    for bins in bin_contig_cov.keys():
        f=open(bins+'_coverage_matrix.txt', 'w')
        f.write('Bin'+'\t'+'Coverage'+'\n')
        for contigs in bin_contig_cov[bins].keys():
            f.write(str(contigs)+'\t'+str(bin_contig_cov[bins][contigs])+'\n')
        f.close()
    os.chdir(pwd)
    return bin_contig_cov, bin_contigs, contig_cov

def binset_comparitor(binset1, binset2, coverage1, coverage2, binset_coverage_total, pwd, num_threads):
    bins_binset1, bins_binset2 = {}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset1):
        os.chdir(pwd+'/'+binset1)
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                bins_binset1[file]=''
    
    for root, dirs, files in os.walk(pwd+'/'+binset2):
        os.chdir(pwd+'/'+binset2)
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                bins_binset2[file]=''
    os.chdir(pwd)

    try:
        os.mkdir('bin_comparison_folder')
    except:
        print('bin_comparison_folder existed')

    f_assem1=open('Merge_binset1_contigs.fa','w')
    f_assem1_red=open('Merge_binset1_contigs_redudant.fa','w')
    os.chdir(str(pwd)+'/'+binset1)
    bin_num, n, contig_len, change_id, seq_de_r = {}, 0, {}, {}, {}
    for bin1 in bins_binset1.keys():
        n+=1
        bin_num[n]=str(bin1)
        for record in SeqIO.parse(bin1,'fasta'):
            try:
                seq_de_r[str(bin1)+'|'+str(record.id)]+=1
                f_assem1_red.write('>'+str(bin1)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
            except:
                seq_de_r[str(bin1)+'|'+str(record.id)]=0
                f_assem1.write('>'+str(bin1)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
            # try:
            #     seq_id=str(record.id).split('|')[0]
            # except:
            #     seq_id=str(record.id)
            contig_len[record.id]=len(record.seq)
    f_assem1.close()
    os.chdir(pwd)

    f_assem2=open('Merge_binset2_contigs.fa','w')
    f_assem2_red=open('Merge_binset2_contigs_redudant.fa','w')
    os.chdir(str(pwd)+'/'+binset2)
    bin_num, n = {}, 0
    for bin2 in bins_binset2.keys():
        n+=1
        bin_num[n]=str(bin2)
        for record in SeqIO.parse(bin2,'fasta'):
            try:
                seq_de_r[str(bin2)+'|'+str(record.id)]+=1
                f_assem2_red.write('>'+str(bin2)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
            except:
                seq_de_r[str(bin2)+'|'+str(record.id)]=0
                f_assem2.write('>'+str(bin2)+'|'+str(record.id)+'\n'+str(record.seq)+'\n')
            # try:
            #     seq_id=str(record.id).split('|')[0]
            # except:
            #     seq_id=str(record.id)
            contig_len[record.id]=len(record.seq)
    f_assem2.close()
    os.chdir(pwd)

    contigs_coverage={}
    n1=0
    for line in open(coverage1, 'r'):
        n1+=1
        if n1 == 1:
            num=str(line).strip().count('drange')
        else:
            # try:
            #     ids=str(line).strip().split('\t')[0].split('|')[0] ### ids needs to be eltered 
            # except:
            ids=str(line).strip().split('\t')[0]
            contigs_coverage[str(ids)]={}
            i=1
            while i <= int(num):
                contigs_coverage[str(ids)][i]=str(line).strip().split('\t')[3*i+1]
                i+=1

    n1=0
    for line in open(coverage2, 'r'):
        n1+=1
        if n1 == 1:
            num=str(line).strip().count('drange')
        else:
            # try:
            #     ids=str(line).strip().split('\t')[0].split('|')[0]  ### ids needs to be eltered 
            # except:
            ids=str(line).strip().split('\t')[0]
            contigs_coverage[str(ids)]={}
            i=1
            while i <= int(num):
                contigs_coverage[str(ids)][i]=str(line).strip().split('\t')[3*i+1]
                i+=1

    try:
        os.system('makeblastdb -in Merge_binset1_contigs.fa -dbtype nucl -hash_index -parse_seqids -logfile temp_db.txt')
    except:
        os.system('makeblastdb -in Merge_binset1_contigs.fa -dbtype nucl -logfile temp_db.txt')
    
    os.system('blastn -db Merge_binset1_contigs.fa -query Merge_binset2_contigs.fa -outfmt 6 -out binset2_vs_binset1.txt -num_threads '+str(num_threads))

    blast_output=open('Filtrated_binset2_vs_binset1.txt','w')
    bin1_cov_ave, bin2_cov_ave, bins_blast_file = {}, {}, {}
    # contig_num = 0
    for line in open('binset2_vs_binset1.txt','r'):
        simi=str(line).strip().split('\t')[2]
        if float(simi) >= 99:
            length=eval(str(line).strip().split('\t')[3])
            contig2=str(line).strip().split('\t')[0]
            if '||' not in contig2:
                org_contig2_id=str(contig2).split('|')[1].strip()
            else:
                contig2a=contig2.replace('||','**')
                org_contig2_id=str(contig2a).split('|')[1].strip().replace('**','||')

            contig1=str(line).strip().split('\t')[1]
            if '||' not in contig1:
                org_contig1_id=str(contig1).split('|')[1].strip()
            else:
                contig1a=contig1.replace('||','**')
                org_contig1_id=str(contig1a).split('|')[1].strip().replace('**','||')

            contig2_length=contig_len[org_contig2_id]
            contig1_length=contig_len[org_contig1_id]
            iden2=float(length)/float(contig2_length)
            iden1=float(length)/float(contig1_length)
            iden_t=iden1+iden2
            if int(length) >= 1000 and iden_t >= 1:
                bin2=str(contig2).split('|')[0].strip()
                bin1=str(contig1).split('|')[0].strip()
                blast_output.write(str(line))
                re_line=str(line).replace(contig2, org_contig2_id).replace(contig1, org_contig1_id)
                bins_vs='Filtrated_'+str(bin2)+'_vs_'+str(bin1)+'.txt'
                if bins_vs not in bins_blast_file.keys():
                    f_bins_blast=open(bins_vs,'w')
                    bins_blast_file[bins_vs]=0
                else:
                    f_bins_blast=open(bins_vs,'a')
                f_bins_blast.write(str(re_line))
                f_bins_blast.close()
    blast_output.close()

    fabs=open('Abnor_normalization_record.txt','w')
    fabs.write('Bin1'+'\t'+'Bin2'+'\t'+'Coverage index'+'\t'+'Coverage1'+'\t'+'Coverage2'+'\n')
    fabs.close()
    try:
        os.mkdir('Bins_blast_output')
    except:
        print('Bins_blast_output existed')
    
    bin_group={}
    for bins_blast_output in bins_blast_file.keys():
        bin2=str(bins_blast_output).split('_vs_')[0].split('Filtrated_')[1]
        try:
            bin_group[bin2][bins_blast_output]=0
        except:
            bin_group[bin2]={}
            bin_group[bin2][bins_blast_output]=0

    ### Waiting for changing into multiple threads
    pool=Pool(processes=num_threads)
    results={}
    for target_bin in bin_group.keys():
        bin_dict=bin_group[target_bin]
        results[target_bin+'_bins_depth_comparison.txt']=0
        pool.apply_async(bin_depth_normalization, args=(target_bin, bin_dict, binset_coverage_total, contigs_coverage, pwd, 1))
    pool.close()
    pool.join()

    f=open('bins_depth_comparison.txt','w')
    for item in results.keys():
        try:
            for line in open(item,'r'):
                f.write(line)
            os.system('mv '+str(item)+' '+pwd+'/bin_comparison_folder')
        except:
            xxx=0
    f.close()

    potential_simialr_bin={}
    try:
        f=open('Potential_similar_bins.txt','a')
        f2=open('Unfiltrated_potential_similar_bins.txt','a')
    except:
        f=open('Potential_similar_bins.txt','w')
        f2=open('Unfiltrated_potential_similar_bins.txt','w')
    for line in open('bins_depth_comparison.txt','r'):
        xxx=0
        bin1=str(line).strip().split('\t')[0]
        bin1_nor_cov=str(line).strip().split('\t')[1].replace('{','').replace('}','')
        bin1_nor_cov_list=bin1_nor_cov.split(',')
        bin2=str(line).strip().split('\t')[2]
        bin2_nor_cov=str(line).strip().split('\t')[3].replace('{','').replace('}','')
        bin2_nor_cov_list=bin2_nor_cov.split(',')
        delta_dict={}
        for i in range(0, len(bin1_nor_cov_list)):
            cov1=float(bin1_nor_cov_list[i].split(':')[1].strip())
            cov2=float(bin2_nor_cov_list[i].split(':')[1].strip())
            delta=abs(cov1-cov2)
            sum=cov1+cov2
            if sum == 0:
                sum=1
            perc=100*delta/sum
            # if cov1 > cov2:
            #     perc=100*delta/cov2
            # else:
            #     perc=100*delta/cov1
            
            if perc <= 10:
                xxx+=1
                delta_dict[i+1]=perc
        
        if xxx == len(bin1_nor_cov_list):
            f.write(str(line).strip()+'\t'+str(delta_dict)+'\n')
            f2.write(str(line).strip()+'\t'+str(delta_dict)+'\n')
            potential_simialr_bin[bin1+'\t'+bin2]=str(line).strip()+'\t'+str(delta_dict)
        else:
            f2.write(str(line).strip()+'\t'+str(delta_dict)+'\n')
    f.close()
    f2.close()
    return potential_simialr_bin

def multiple_assembly_comparitor_main(Contig_list, BestBinSet_list, Coverage_list, datasets, step, num_threads):
    pwd=os.getcwd()
    print('Processing contigs: '+str(Contig_list))
    print('Processing bisnet: '+str(BestBinSet_list))
    print('Processing coverage: '+str(Coverage_list))
    print('Processing datasets: '+str(datasets))
    contig_num=len(Contig_list)

    ### Forming core-contigs of bins
    if step == 'initial_drep':
        fx=open('Core_contigs_total_bins.txt','w')
        CC_binset_list, CC_contig_list, binset_coverage_dict, binset_coverage_total, total_core_contigs =[], [], {}, {}, {}
        for i in range(0, len(BestBinSet_list)):
            binset_folder=BestBinSet_list[i]
            coverage_file=Coverage_list[i]
            contig_file=Contig_list[i]
            A=record_bin_coverage(binset_folder, coverage_file)
            bin_contig_cov=A[0]
            bin_contig=A[1]
            contig_cov=A[2]
            core_contigs, core_contigs_IQR_ave_coverage=core_contigs_filtration(bin_contig_cov, bin_contig, contig_cov, binset_folder, contig_file)
            binset_coverage_dict['CC_'+contig_file]=core_contigs_IQR_ave_coverage
            binset_coverage_total.update(core_contigs_IQR_ave_coverage)
            CC_binset_list.append('CC_'+binset_folder)
            CC_contig_list.append('CC_'+contig_file)
            f=open('Core_contigs_'+str(binset_folder)+'.txt','w')
            for bin_id in core_contigs.keys():
                f.write(str(bin_id)+'\t'+str(core_contigs[bin_id])+'\n')
                fx.write(str(bin_id)+'\t'+str(core_contigs[bin_id])+'\n')
            f.close()
            total_core_contigs.update(core_contigs)
        fx.close()

        if contig_num != 1:
            for i in range(1, contig_num):
                potential_simialr_bin=binset_comparitor(BestBinSet_list[0], BestBinSet_list[-1], Coverage_list[0], Coverage_list[-1], binset_coverage_total, pwd, num_threads)
                binset_checkm_connection_1=checkm_connections(BestBinSet_list[0])
                binset_checkm_connection_2=checkm_connections(BestBinSet_list[-1])

                selected_bins, eliminated_bins=bin_comparitor(potential_simialr_bin, binset_checkm_connection_1, binset_checkm_connection_2, i)
                binset_record, binset_genome_size={}, {}
                A_group=genome_contigs_recorder(BestBinSet_list[0], binset_record, binset_genome_size, Coverage_list[0])
                binset_record1=A_group[0]
                binset_genome_size1=A_group[1]
                binset_coverage_avg1=A_group[2]
                binset_GC_ratio1=A_group[3]
                all_bins_1=A_group[4]

                B_group=genome_contigs_recorder(BestBinSet_list[-1], binset_record, binset_genome_size, Coverage_list[-1])
                binset_record2=B_group[0]
                binset_genome_size2=B_group[1]
                binset_coverage_avg2=B_group[2]
                binset_GC_ratio2=B_group[3]
                all_bins_2=B_group[4]

                F=new_selected_bins_generator(selected_bins, eliminated_bins, all_bins_1, all_bins_2, binset_checkm_connection_1, binset_checkm_connection_2, i, BestBinSet_list[0], BestBinSet_list[-1], Contig_list[0], Contig_list[-1], Coverage_list[0], Coverage_list[-1], binset_coverage_avg1, binset_coverage_avg2)    
                new_Contigs=F[2]
                new_bin_folder=F[3]
                new_coverage=F[4]
                bin_avg_cov_new=F[5]

                Contig_list.remove(Contig_list[0])
                Contig_list.remove(Contig_list[-1])
                Contig_list.append(new_Contigs)
                print('Adding '+str(new_Contigs)+' to CC Contigs list')
                BestBinSet_list.remove(BestBinSet_list[0])
                BestBinSet_list.remove(BestBinSet_list[-1])
                BestBinSet_list.append(new_bin_folder)
                print('Adding '+str(new_bin_folder)+' to BestBinSet list')
                Coverage_list.remove(Coverage_list[0])
                Coverage_list.remove(Coverage_list[-1])
                Coverage_list.append(new_coverage)
                print('Adding '+str(new_coverage)+' to coverage list')
                binset_coverage_dict[new_Contigs]=bin_avg_cov_new
                print('Adding '+str(bin_avg_cov_new)+' of '+str(new_Contigs)+' to binset coverage dict')
        else:
            new_Contigs=Contig_list[-1]
            new_bin_folder=BestBinSet_list[-1]

        # new_bin_folder='Iteration_5_genomes'
        initial_drep_final_comparitor(new_bin_folder, Coverage_list, datasets, num_threads, pwd, step)
        os.system('rm -rf Iteration_*')

    elif step == 'second_drep' or step == 'final_drep':
        if len(BestBinSet_list) != 1:
            raise ValueError(
                "{} requires exactly one final binset folder".format(step)
            )
        final_binset_comparitor(BestBinSet_list[0], Coverage_list, datasets, num_threads, pwd, step)

    print(str(step)+' bin de-replacation done!')

if __name__ == '__main__': 
    step='final_drep' 
    ### 'initial_drep': remove redudant bins in S4; 'second_drep': remove redudant bins in the following step with regonized pair-redudant bins; 'final_drep': remove redudant bins after reassembly if paired-bins still present;
    ### Initial de-replication should put all the BestBinset folder of all assemblies into a list
    ### 'second_drep' and 'final_drep': put the single folder into the list
    BestBinSet_list=['1_RH_insert_cat_270_final.contigs.fa_BestBinsSet','2_RH_S001_insert_270_final.contigs.fa_BestBinsSet','3_RH_S002_insert_270_final.contigs.fa_BestBinsSet','4_RH_S003_insert_270_final.contigs.fa_BestBinsSet','5_RH_S004_insert_270_final.contigs.fa_BestBinsSet','6_RH_S005_insert_270_final.contigs.fa_BestBinsSet'] ### Beaware the ID of assembly should be different
    Contig_list=['1_RH_insert_cat_270_final.contigs.fa','2_RH_S001_insert_270_final.contigs.fa','3_RH_S002_insert_270_final.contigs.fa','4_RH_S003_insert_270_final.contigs.fa','5_RH_S004_insert_270_final.contigs.fa','6_RH_S005_insert_270_final.contigs.fa'] ##
    Coverage_list=['Coverage_matrix_for_binning_1_8_medium_cat_SPAdes_scaffolds.fasta.txt','Coverage_matrix_for_binning_2_9_medium_S001_SPAdes_scaffolds.fasta.txt','Coverage_matrix_for_binning_3_10_medium_S002_SPAdes_scaffolds.fasta.txt']
    datasets={'1':['PE_r1_RM2_S001__insert_270_mate1.fq','PE_r2_RM2_S001__insert_270_mate2.fq'],'2':['PE_r1_RM2_S001__insert_270_mate1.fq','PE_r2_RM2_S001__insert_270_mate2.fq']}
    
    num_threads=30
    final_folder_list='BestBinset_outlier_refined_filtrated_retrieved_retrieved_OLC_re-assembly_OLC'
    # multiple_assembly_comparitor_main(Contig_list, BestBinSet_list, Coverage_list, datasets, step, num_threads)
    pwd=os.getcwd()
    final_binset_comparitor(final_folder_list, Coverage_list, datasets, num_threads, pwd, step)
