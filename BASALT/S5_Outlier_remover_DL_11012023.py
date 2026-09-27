#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S5 (CheckM2 branch): Deep-learning-based outlier bin removal.

This module applies trained DL models to identify and remove
contaminated or low-quality bins, preparing cleaner binsets for
subsequent contig retrieval and reassembly steps.
"""

from Bio import SeqIO
import os, copy, math, glob, gc
import numpy as np
from multiprocessing import Pool
from basalt_runtime import coerce_count, depth_sample_count


def TNF_coverage_matrix(bin_contigs, bin_id, contigs_depth, ccc,
                        contigs_kmer, contigs_kmer2):
    """
    Write TNF and coverage features for all contigs within a bin.

    Parameters
    ----------
    bin_contigs : dict
        Mapping bin_id -> {contig_id -> contig_length}.
    bin_id : str
        Current bin identifier.
    contigs_depth : dict
        Mapping contig_id -> {cov_index -> depth}.
    ccc : dict
        Mapping bin_id -> contig_id -> coverage change ratios.
    contigs_kmer : dict
        Mapping contig_id -> k-mer feature vector.
    contigs_kmer2 : dict
        Helper dict used to mark contigs processed.
    """
    fout=open(str(bin_id)+'_contig_depth_TNF_matrix.txt','a')
    for contigs in bin_contigs[bin_id].keys():
        try:
            contigs_kmer2[contigs]=0
            try:
                fout.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+str(ccc[bin_id][contigs])+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
                # fx.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+str(ccc[bin_id][contigs])+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
            except: 
                fout.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+''+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
                # fx.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+''+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
        except:
            xyzzz=0
    fout.close()

def bin_kmer(kmer, contigs_bin, contigs_bin2, bin_contigs,
             contigs_depth, ccc):
    """
    Distribute k-mer features from a global k-mer file into per-bin matrices.

    Parameters
    ----------
    kmer : str
        Path to the k-mer feature file.
    contigs_bin : dict
        Mapping contig_id -> list of bins containing that contig.
    contigs_bin2 : dict
        Mapping contig_id -> number of bins containing it.
    bin_contigs : dict
        Mapping bin_id -> {contig_id -> contig_length}.
    contigs_depth : dict
        Mapping contig_id -> {cov_index -> depth}.
    ccc : dict
        Mapping bin_id -> contig_id -> coverage change ratios.
    """
    n=0
    for line in open(kmer,'r'):
        n+=1
        if n >= 2:
            contig_id=str(line).strip().split('\t')[0]
            kmers=str(line).strip().split('\t')[1:]
            try:
                contigs_bin2[contig_id]+=1
                for bin_id in contigs_bin[contig_id]:
                    fout=open(str(bin_id)+'_contig_depth_TNF_matrix.txt','a')
                    try:
                        fout.write(contig_id+'\t'+str(bin_contigs[bin_id][contig_id])+'\t'+str(contigs_depth[contig_id])+'\t'+str(ccc[bin_id][contig_id])+'\t'+str(str(kmers).replace('\'',''))+'\n')
                        # fx.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+str(ccc[bin_id][contigs])+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
                    except: 
                        fout.write(contig_id+'\t'+str(bin_contigs[bin_id][contig_id])+'\t'+str(contigs_depth[contig_id])+'\t'+''+'\t'+str(str(kmers).replace('\'',''))+'\n')
                        # fx.write(contigs+'\t'+str(bin_contigs[bin_id][contigs])+'\t'+str(contigs_depth[contigs])+'\t'+''+'\t'+str(str(contigs_kmer[contigs]).replace('\'',''))+'\n')
                    fout.close()
            except:
                xyzzz=0

def basic_information(bin_folder, depth_list, kmer_list, pwd, num_threads):
    """
    Collect basic information needed for DL outlier detection.

    Parameters
    ----------
    bin_folder : str
        Folder containing bin FASTA files.
    depth_list : list of str
        List of depth matrix filenames.
    kmer_list : list of str
        List of k-mer feature file paths.
    pwd : str
        Working directory.
    num_threads : int
        Number of processes to use when distributing k-mer features.

    Returns
    -------
    tuple
        (depth_TNF_matrix_file, contigs_depth, ccc, num_cov, bin_contigs_seq)
    """
    print('Reading basic information')
    # f=open('BestBinset.fasta','w')
    bin_contigs, bin_contigs_seq, total_contigs, bin_name, contigs_bin, contigs_bin2 = {}, {}, {}, {}, {}, {}
    os.chdir(pwd+'/'+bin_folder)
    for root, dirs, files in os.walk(pwd+'/'+bin_folder):
        for file in files:
            hz=file.split('.')[-1]
            if 'fa' in hz:
                bin_contigs[file]={}
                bin_contigs_seq[file]={}
                bn=('.').join(file.split('.')[:-1])
                # print(bn)
                bin_name[bn]=file
                for record in SeqIO.parse(file,'fasta'):
                    bin_contigs[file][record.id]=str(len(record.seq))
                    bin_contigs_seq[file][record.id]=str(record.seq)
                    
                    try:
                        contigs_bin[record.id].append(file)
                        contigs_bin2[record.id]+=1
                    except:
                        contigs_bin[record.id]=[file]
                        contigs_bin2[record.id]=1
                    if record.id not in total_contigs.keys():
                        # f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                        total_contigs[record.id]=0
    os.chdir(pwd)
    # f.close()

    # os.system('perl calc.kmerfreq.pl -i BestBinset.fasta -o BestBinset.kmer.txt')

    contigs_depth = {}
    for depth_file in depth_list:
        print('Processing '+str(depth_file))
        n=0
        for line in open(depth_file,'r'):
            n+=1
            if n == 1:
                num=str(line).strip().count('drange')
                # nx=int(int(num)/2)
            else:
                contig_id=str(line).strip().split('\t')[0]
                contigs_depth[contig_id]={} 
                for i in range(1,num+1):
                    # contigs_depth[contig_id][i+1]=float(str(line).strip().split('\t')[3+2*i])
                    contigs_depth[contig_id][i]=float(str(line).strip().split('\t')[3*i+1])

    ccc={}
    if num >= 2:
        for bin_id in bin_contigs.keys():
            # print('Processing '+str(bin_id))
            ccc[bin_id]={}
            for c_id in bin_contigs[bin_id].keys():
                cx=0
                ccc[bin_id][c_id]={}
                # print('Processing '+str(c_id))
                for i in range(2,num+1):
                    cx+=1
                    try:
                        ccc[bin_id][c_id][cx]=contigs_depth[c_id][i]/contigs_depth[c_id][1]
                    except:
                        ccc[bin_id][c_id][cx]=contigs_depth[c_id][i]
    print(str(kmer_list))
    
    for bin_id in bin_contigs.keys():
        fout=open(str(bin_id)+'_contig_depth_TNF_matrix.txt','w')
        fout.write(bin_id+'\n')
        fout.close()

    pool=Pool(processes=num_threads)
    for kmer in kmer_list:
        # bin_kmer(kmer, contigs_bin, contigs_bin2, bin_contigs, contigs_depth, ccc)
        pool.apply_async(bin_kmer, args=(kmer, contigs_bin, contigs_bin2, bin_contigs, contigs_depth, ccc,))
    pool.close()
    pool.join()

    fx=open('Bin_contig_depth_TNF_matrix.txt','w')
    for bin_id in bin_contigs.keys():
        for line in open(str(bin_id)+'_contig_depth_TNF_matrix.txt','r'):
            if len(line) != 0:
                fx.write(line.strip()+'\n')
        os.system('rm '+str(bin_id)+'_contig_depth_TNF_matrix.txt')
    fx.close() 
    os.system('mkdir BestBinset_outlier_test')
    os.system('mv Bin_contig_depth_TNF_matrix.txt BestBinset_outlier_test')
    # os.system('rm BestBinset.fasta BestBinset.kmer.txt')
    # os.chdir(pwd+'/BestBinset_outlier_test')
    # os.system('cat *_contig_depth_TNF_matrix.txt > Bin_contig_depth_TNF_matrix2.txt')
    # os.system('rm *_contig_depth_TNF_matrix.txt')
    # os.system('mv Bin_contig_depth_TNF_matrix2.txt Bin_contig_depth_TNF_matrix.txt')
    # os.chdir(pwd)
    return 'Bin_contig_depth_TNF_matrix.txt', contigs_depth, ccc, num, bin_contigs_seq

def outlier_predictor(depth_TNF_matrix, contigs_depth, bin_contigs,
                      datasets, lr, hifi_list, num_threads, nx):
    """
    Run DL-based outlier prediction and summarise contig/bins statistics.

    Parameters
    ----------
    depth_TNF_matrix : str
        Path to the combined depth+TNF feature matrix.
    contigs_depth : dict
        Mapping contig_id -> {cov_index -> depth}.
    bin_contigs : dict
        Mapping bin_id -> {contig_id -> sequence}.
    datasets : dict
        Mapping dataset_id -> [R1, R2] short-read files.
    lr : list of str
        Long-read datasets (ONT/PB).
    hifi_list : list of str
        HiFi read datasets.
    num_threads : int
        Number of threads for remapping and depth testing.
    nx : int
        Unused parameter kept for compatibility.

    Returns
    -------
    tuple
        (confirmed_outlier, contig_clean_status, bin_avg_coverage,
         bin_core_contig_num)
    """
    print('Predicting potential outliers')
    os.system('ensemble.py')
    # outlier_predictor.py
    contig_status, bin_positive_contigs, bin_avg_coverage, bin_core_contig_num, contig_clean_status = {}, {}, {}, {}, {}
    for line in open('Predicted_potential_outlier.txt','r'):
        num=str(line).strip().count('\t')
        if num == 0:
            bin_id=str(line).strip()
            bin_core_contig_num[bin_id]=0
            # contig_status[bin_id]={}
            # contig_status[bin_id]['Real']=[]
            contig_status[bin_id]=[]
            contig_clean_status[bin_id]={}
            bin_positive_contigs[bin_id]={}
        else:
            status=str(line).strip().split('\t')[0]
            contig_id=str(line).strip().split('\t')[1]
            # coverage_matrix=str(line).strip().split('\t')[3]
            if status == 'Real':
                bin_core_contig_num[bin_id]+=1
                for cov in contigs_depth[contig_id].keys():
                    try:
                        bin_positive_contigs[bin_id][cov].append(contigs_depth[contig_id][cov])
                    except:
                        bin_positive_contigs[bin_id][cov]=[]
                        bin_positive_contigs[bin_id][cov].append(contigs_depth[contig_id][cov])
                contig_clean_status[bin_id][contig_id]=contigs_depth[contig_id]
            elif status == 'Contaminated':
                contig_status[bin_id].append(contig_id)
    
    ft=open('Outlier_test_file.txt','w')
    ft.write('Coverage Index'+'\t'+'Coverage list'+'\t'+'Q1'+'\t'+'AVG'+'\t'+'Q3'+'\t'+'Upper filtration value'+'\t'+'Lower filtration value'+'\n')
    bin_positive_contigs_IQR = {}
    # bin_positive_contigs_IQR, bin_avg_coverage_media = {}, {}
    for bin_id in bin_positive_contigs.keys():
        bin_positive_contigs_IQR[bin_id]={}
        bin_avg_coverage[bin_id]={}
        # bin_avg_coverage_media[bin_id]={}
        for cov in bin_positive_contigs[bin_id].keys():
            if len(bin_positive_contigs[bin_id][cov]) >= 1: ### Contig number
                bin_positive_contigs_IQR[bin_id][cov]={}
                a=np.array(bin_positive_contigs[bin_id][cov])
                Q1=np.percentile(a,25)
                Q3=np.percentile(a,75)
                ### consider the number of contigs if the number of contigs < = 5
                # Q2=np.percentile(a,50)
                avg=np.mean(a)
                IQR=Q3-Q1
                upper = Q3 + IQR
                lower = Q1 - IQR
                bin_positive_contigs_IQR[bin_id][cov]['lower']=lower
                bin_positive_contigs_IQR[bin_id][cov]['upper']=upper
                bin_avg_coverage[bin_id][cov]=avg
                # bin_avg_coverage_media[bin_id][cov]=Q2
                ft.write(str(bin_id)+'\t'+str(cov)+'\t'+str(a)+'\t'+str(Q1)+'\t'+str(avg)+'\t'+str(Q3)+'\t'+str(upper)+'\t'+str(lower)+'\n')
    ft.close()
    
    # print(str(bin_positive_contigs_IQR))

    f=open('Bin_core_seq_avg_depth.txt','w')
    for bin_id in bin_avg_coverage.keys():
        f.write(str(bin_id))
        for cov in bin_avg_coverage[bin_id].keys():
            bin_avg_coverage[bin_id][cov]=bin_avg_coverage[bin_id][cov]/bin_core_contig_num[bin_id]
            f.write('\t'+str(bin_avg_coverage[bin_id][cov]))
        f.write('\n')
    f.close()

    # print(str(contig_status))

    f=open('Potential_contaminted_seq_vari.txt','w')
    f2=open('Remapping.fasta','w')
    fx=open('Remapped_depth_test.txt','w')
    # f2=open('Confirmed_potential_contaminted_seq_vari.txt','w')
    f.write('Judge_status'+'\t'+'Bin_id'+'\t'+'Contig_id'+'\t'+'Contig Cov.'+'\t'+'Cov Vari.(%)'+'\n')
    confirmed_outlier, t_contigs = {}, {}
    for bin_id in contig_status.keys():
        confirmed_outlier[bin_id]=[]
        if len(contig_status[bin_id]) != 0:
            for contig_id in contig_clean_status[bin_id].keys():
                if contig_id not in contig_status[bin_id]:
                    f.write('Clean'+'\t'+str(bin_id)+'\t'+str(contig_id)+'\t'+str(contig_clean_status[bin_id][contig_id])+'\n')

            print('Judging '+str(bin_id))
            for contig_id in contig_status[bin_id]:
                try:
                    judge=0
                    tt=str(bin_id)+'\t'+str(contig_id)+'\t'+str(contigs_depth[contig_id])
                    for cov in contigs_depth[contig_id].keys():
                        d=contigs_depth[contig_id][cov]
                        ad=bin_avg_coverage[bin_id][cov]
                        # am=bin_avg_coverage_media[bin_id][cov]
                        if ad != 0:
                            vari=100*abs(d-ad)/ad
                            # vari2=100*abs(d-am)/am
                            if vari < 30:
                            # if vari < 30 or vari2 < 30:
                                judge+=1
                            else:
                                try:
                                    if d >= bin_positive_contigs_IQR[bin_id][cov]['lower'] and d <= bin_positive_contigs_IQR[bin_id][cov]['upper']:
                                        judge+=1
                                except:
                                    xyz=0
                        else:
                            judge+=1
                        tt+='\t'+str(cov)+':'+str(vari)
                    tt+='\n'
                    
                    if judge == nx:
                        f.write('Clean'+'\t'+tt)
                        print('Clean '+str(contig_id)+':'+str(judge))
                    else:
                        # f.write('Putative contamination'+'\t'+tt)
                        try:
                            t_contigs[contig_id]+=1
                        except:
                            t_contigs[contig_id]=0
                            f2.write('>'+str(contig_id)+'\n'+str(bin_contigs[bin_id][contig_id])+'\n')
                        confirmed_outlier[bin_id].append(contig_id)
                        fx.write(str(contig_id)+'\t'+str(contigs_depth[contig_id])+'\n')
                        # print('Contaminated '+str(contig_id)+':'+str(judge))
                except:
                    xyz=0

    f2.close()
    fx.close()
    # f.close()

    ### Re-mapping
    print('Re-mapping')
    total_fa='Remapping.fasta'
    if len(datasets) != 0 or len(hifi_list) != 0:
        os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))

    ###
    bam_sorted1=''
    if len(datasets) != 0:        
        n = 0
        for item in datasets.keys():
            print('Mapping '+str(datasets[item]))
            n+=1
            os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets[item][0])+' -2 '+str(datasets[item][1])+' -S '+str(item)+'.sam -q --no-unal')
            os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(item)+'.sam -o '+str(item)+'.bam')
            # py2
            os.system('samtools sort -@ '+str(num_threads)+' '+str(item)+'.bam '+str(item)+'_sorted') 

            try:
                with open(str(item)+'_sorted.bam', 'r') as fh:
                    pass
            except FileNotFoundError:
                print('samtools sorting '+str(item)+'.bam failed. Redoing')
                # py3
                os.system('samtools sort -@ '+str(num_threads)+' -o '+str(item)+'_sorted.bam '+str(item)+'.bam')
            os.system('rm '+str(item)+'.sam')

            if item == '1':
                bam_sorted1='1_sorted.bam'
            else:
                bam_sorted1+=' '+str(item)+'_sorted.bam'
    
    bam_sorted2=''
    if len(lr) != 0 or len(hifi_list) != 0:
        print('Mapping '+str(lr)+' to contigs/scaffolds')
        if len(hifi_list) != 0:
            hifi_s=[]
            for item in hifi_list:
                name_l=str(item).split('.')
                name_l.remove(name_l[-1])
                name='.'.join(name_l)
                hifi_s.append(name+'_split_R1.fa')
                hifi_s.append(name+'_split_R2.fa')

            for i in range(1,len(hifi_list)+1):
                print(str(hifi_list[i-1])+' mapped against Remapped.fasta')
                os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(hifi_s[2*i-2])+' -2 '+str(hifi_s[2*i-1])+' -S hifi'+str(i)+'.sam -f --no-unal')
                os.system('samtools view -@ '+str(num_threads)+' -b -S hifi'+str(i)+'.sam -o hifi'+str(i)+'.bam')
                os.system('rm hifi'+str(i)+'.sam')

                print('Sorting bam file')
                ### py2
                os.system('samtools sort -@ '+str(num_threads)+' hifi'+str(i)+'.bam hifi'+str(i)+'_sorted') 

                try:
                    with open('hifi'+str(i)+'_sorted.bam', 'r') as fh:
                        pass
                except FileNotFoundError:
                    print('samtools sorting hifi'+str(i)+'.bam failed. Redoing')
                    ### py3
                    os.system('samtools sort -@ '+str(num_threads)+' -o hifi'+str(i)+'_sorted.bam hifi'+str(i)+'.bam' )
            
                if i == 1:
                    bam_sorted2='hifi1_sorted.bam'
                else:
                    bam_sorted2+=' hifi'+str(i)+'_sorted.bam'


        if len(lr) != 0:
            for i in range(1, len(lr)+1):
                # if lr_type == 'ont': 
                os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Remapping.fasta '+str(lr[i-1])+' > lr'+str(i)+'.sam')
                # elif lr_type == 'pb': 
                #     os.system('minimap2 -t '+str(num_threads)+' -ax map-pb Remapping.fasta '+str(lr[i-1])+' > lr'+str(i)+'.sam')
                # elif lr_type == 'hifi': 
                #     os.system('minimap2 -t '+str(num_threads)+' -ax map-hifi Remapping.fasta '+str(lr[i-1])+' > lr'+str(i)+'.sam')
                print(str(lr[i-1])+' mapped against Remapped.fasta')
                os.system('samtools view -@ '+str(num_threads)+' -b -S lr'+str(i)+'.sam -o lr'+str(i)+'.bam')
                os.system('rm lr'+str(i)+'.sam')

                print('Sorting bam file')
                ### py2
                os.system('samtools sort -@ '+str(num_threads)+' lr'+str(i)+'.bam lr'+str(i)+'_sorted') 

                try:
                    with open('lr'+str(i)+'_sorted.bam', 'r') as fh:
                        pass
                except FileNotFoundError:
                    print('samtools sorting lr'+str(i)+'.bam failed. Redoing')
                    ### py3
                    os.system('samtools sort -@ '+str(num_threads)+' -o lr'+str(i)+'_sorted.bam lr'+str(i)+'.bam' )
            
                if i == 1:
                    bam_sorted2='lr1_sorted.bam'
                else:
                    bam_sorted2+=' lr'+str(i)+'_sorted.bam'

    if len(bam_sorted1) !=0 and len(bam_sorted2) !=0:
        bam_sorted=bam_sorted1+' '+bam_sorted2
    elif len(bam_sorted1) !=0 and len(bam_sorted2) ==0:
        bam_sorted=bam_sorted1
    elif len(bam_sorted1) ==0 and len(bam_sorted2) !=0:
        bam_sorted=bam_sorted2

    try:
        os.system('jgi_summarize_bam_contig_depths --outputDepth Re-mapped_depth.txt '+str(bam_sorted))
        nxxyy=0
        for line in open('Re-mapped_depth.txt','r'):
            nxxyy+=1
            if nxxyy == 2:
                break
    except:
        os.system('jgi_summarize_bam_contig_depths --outputDepth Re-mapped_depth.txt '+str(bam_sorted))

    os.system('rm '+str(bam_sorted)+' *.bt2')
    
    for item in datasets.keys():
        os.system('rm '+str(item)+'.bam')

    if len(lr) != 0:
        for i in range(1, len(lr)+1):
            os.system('rm lr'+str(i)+'_sorted.bam lr'+str(i)+'.bam')

    contigs_depth_remapped, n={}, 0
    for line in open('Re-mapped_depth.txt','r'):
        n+=1
        if n == 1:
            num=depth_sample_count(line)
            # nx=int(int(num)/2)
        else:
            contig_id=str(line).strip().split('\t')[0]
            contigs_depth_remapped[contig_id]={} 
            for i in range(0,num):
                contigs_depth_remapped[contig_id][i+1]=float(str(line).strip().split('\t')[3+2*i])

    ### Re-judge
    fxy=open('Rejudge_clean.txt','w')
    confirmed_outlier2={}
    for bin_id in confirmed_outlier.keys():
        confirmed_outlier2[bin_id]=[]
        print('Re-judging '+str(bin_id))
        for contig_id in confirmed_outlier[bin_id]:
            try:
                judge=0
                tt=str(bin_id)+'\t'+str(contig_id)+'\t'+str(contigs_depth_remapped[contig_id])
                for cov in contigs_depth_remapped[contig_id].keys():
                    d=contigs_depth_remapped[contig_id][cov]
                    ad=bin_avg_coverage[bin_id][cov]
                    if ad != 0:
                        vari=100*abs(d-ad)/ad
                        if vari < 30:
                            judge+=1
                        else:
                            try:
                                if d >= bin_positive_contigs_IQR[bin_id][cov]['lower'] and d <= bin_positive_contigs_IQR[bin_id][cov]['upper']:
                                    judge+=1
                            except:
                                xyz=0
                    else:
                        judge+=1
                    tt+='\t'+str(cov)+':'+str(vari)
                tt+='\n'
                    
                if judge == nx:
                    f.write('Clean'+'\t'+tt)
                    fxy.write('Clean'+'\t'+tt)
                    print('Clean '+str(contig_id)+':'+str(judge))
                else:
                    f.write('Contaminated'+'\t'+tt)
                    confirmed_outlier2[bin_id].append(contig_id)
                    print('Contaminated '+str(contig_id)+':'+str(judge))
            except:
                    tt=str(bin_id)+'\t'+str(contig_id)+'\t'+str(contigs_depth[contig_id])
                    for cov in contigs_depth[contig_id].keys():
                        tt+='\t'+str(cov)+':'+str(vari)
                    tt+='\n'
                    f.write('Contaminated'+'\t'+tt)
    f.close()
    fxy.close()

    os.system('mkdir BestBinset_outlier_basic_information')
    os.system('mv Predicted_potential_outlier.txt Bin_contig_depth_TNF_matrix.txt Bin_core_seq_avg_depth.txt BestBinset_outlier_basic_information')
    return contig_status, bin_avg_coverage, confirmed_outlier

def checkm_eval(bin_contigs, bin_folder, confirmed_outlier, pwd, num_threads):
    """
    Evaluate outlier removal using CheckM2 for each bin.

    Parameters
    ----------
    bin_contigs : dict
        Mapping bin_id -> {contig_id -> sequence}.
    bin_folder : str
        Original BestBinset folder name.
    confirmed_outlier : dict
        Mapping bin_id -> list of contig_ids marked as outliers.
    pwd : str
        Working directory.
    num_threads : int
        Number of threads for CheckM2.
    """
    outlier_binset=bin_folder+'_outlier_refined'
    os.system('mkdir '+outlier_binset)
    os.chdir(outlier_binset) 
    for bin_id in confirmed_outlier.keys():
        bin_name_list=bin_id.split('.')
        bin_name_list.remove(bin_name_list[-1])
        bin_name='.'.join(bin_name_list)
        f=open(bin_name+'.fa','w')
        for seq_id in bin_contigs[bin_id].keys():
            if seq_id not in confirmed_outlier[bin_id]:
                f.write('>'+str(seq_id)+'\n'+str(bin_contigs[bin_id][seq_id])+'\n')
        f.close()
    os.chdir(pwd)

    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(outlier_binset)+' -x fa -o '+str(outlier_binset)+'_checkm')
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(outlier_binset)+' '+str(outlier_binset)+'_checkm')
    os.chdir(str(outlier_binset)+'_checkm/')
    print('Parsing '+outlier_binset+' checkm output')
    
    refined_checkm, n = {}, 0
    for line in open('quality_report.tsv','r'):
        n+=1
        if n >= 2:
            binID=str(line).strip().split('\t')[0].strip()
            genome_size=str(line).strip().split('\t')[8].strip()
            completeness=str(line).strip().split('\t')[1].strip()
            contamination=str(line).strip().split('\t')[2].strip()
            N50=str(line).strip().split('\t')[6].strip()

            refined_checkm[str(binID)]={}
            refined_checkm[str(binID)]['N50']=int(N50)
            refined_checkm[str(binID)]['Completeness']=float(completeness)
            refined_checkm[str(binID)]['Genome size']=coerce_count(genome_size)
            refined_checkm[str(binID)]['Contamination']=float(contamination)
    
    os.chdir(pwd+'/'+str(bin_folder))
    orig_checkm={}
    for root, dirs, files in os.walk(pwd+'/'+str(bin_folder)):
        for file in files:
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:
                        binID=str(line).strip().split('\t')[0].strip()
                        if binID in refined_checkm.keys():
                            completeness=float(str(line).strip().split('\t')[2].strip())
                            contamination=float(str(line).strip().split('\t')[3].strip())
                            N50=float(str(line).strip().split('\t')[4].strip())
                            genome_size=coerce_count(str(line).strip().split('\t')[1].strip())
                            qua=float(completeness)-5*float(contamination)

                            refined_qua=refined_checkm[str(binID)]['Completeness']-5*refined_checkm[str(binID)]['Contamination']
                            if qua > refined_qua:
                                os.system('rm '+pwd+'/'+str(outlier_binset)+'/'+binID+'.fa')
                                os.system('cp '+binID+'.fa '+pwd+'/'+str(outlier_binset))
                                refined_checkm[str(binID)]['N50']=int(N50)
                                refined_checkm[str(binID)]['Completeness']=float(completeness)
                                refined_checkm[str(binID)]['Genome size']=float(genome_size)
                                refined_checkm[str(binID)]['Contamination']=float(contamination)
                        else:
                            os.system('cp '+binID+'.fa '+pwd+'/'+str(outlier_binset))
                            refined_checkm[str(binID)]={}
                            refined_checkm[str(binID)]['N50']=int(N50)
                            refined_checkm[str(binID)]['Completeness']=float(completeness)
                            refined_checkm[str(binID)]['Genome size']=float(genome_size)
                            refined_checkm[str(binID)]['Contamination']=float(contamination)
    
    os.chdir(pwd+'/'+str(outlier_binset))
    f=open('quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')

    for binID in refined_checkm.keys():
        f.write(str(binID)+'\t'+str(refined_checkm[str(binID)]['Genome size'])+'\t'+str(refined_checkm[str(binID)]['Completeness'])+'\t'+str(refined_checkm[str(binID)]['Contamination'])+'\t'+str(refined_checkm[str(binID)]['N50'])+'\n')
    f.close()
    os.chdir(pwd)

def kmer_cal(input_file, output_file):
    """
    Calculate k-mer features for the given FASTA and write to output file.
    """
    os.system('calc.kmerfreq.pl -i '+str(input_file)+' -o '+str(output_file))

def outlier_remover_main(bin_folder, depth_list, datasets, lr, hifi_list,
                         assemblies_list, pwd, num_threads):
    """
    Entry point for S5 (CheckM2 branch): DL-based outlier removal.

    Parameters
    ----------
    bin_folder : str
        Folder containing original BestBinset bins.
    depth_list : list of str
        List of coverage matrix files used for feature extraction.
    datasets : dict
        Mapping dataset_id -> [R1, R2] short-read files.
    lr : list of str
        Long-read datasets (ONT/PB).
    hifi_list : list of str
        HiFi read datasets.
    assemblies_list : list of str
        List of modified assembly filenames.
    pwd : str
        Working directory.
    num_threads : int
        Number of threads to use.

    Returns
    -------
    None
        Writes refined binsets and intermediate logs to disk.
    """
    # kmer_list=glob.glob(r'*.kmer.txt')
    # if len(kmer_list)==0:
    kmer_list=[]
    pool=Pool(processes=num_threads)
    for i in range(0, len(assemblies_list)):
        print('Parsing kmer', str(assemblies_list[i]))
        kmer_list.append(str(assemblies_list[i])+'.kmer.txt')
        pool.apply_async(kmer_cal, args=(str(assemblies_list[i]), str(assemblies_list[i])+'.kmer.txt',))
    pool.close()
    pool.join()

    A=basic_information(bin_folder, depth_list, kmer_list, pwd, num_threads)
    depth_TNF_matrix, contigs_depth, ccc, nx, bin_contigs=A[0], A[1], A[2], A[3], A[4]
    A=outlier_predictor(depth_TNF_matrix, contigs_depth, bin_contigs, datasets, lr, hifi_list, num_threads, nx)
    contig_status, bin_avg_coverage, confirmed_outlier= A[0], A[1], A[2]
    checkm_eval(bin_contigs, bin_folder, confirmed_outlier, pwd, num_threads)
    del contigs_depth, ccc, bin_contigs, contig_status, bin_avg_coverage, confirmed_outlier ### Release ram
    gc.collect()
    print('Outlier removal done!')

if __name__ == '__main__': 
    bin_folder='BestBinset'
    depth_list=['Coverage_matrix_for_binning_1_HumanGut_cat.fasta.txt']
    # depth_list=['Coverage_matrix_for_binning_1_RH_insert_cat_270_final.contigs.fa.txt','Coverage_matrix_for_binning_2_RH_S001_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_3_RH_S002_insert_270_final.contigs.fa.txt',
    # 'Coverage_matrix_for_binning_4_RH_S003_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_5_RH_S004_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_6_RH_S005_insert_270_final.contigs.fa.txt']
    # datasets={'1':['RH_S001_insert_270_mate1.fq','RH_S001_insert_270_mate2.fq'], '2':['RH_S002_insert_270_mate1.fq','RH_S002_insert_270_mate2.fq'], '3':['RH_S003_insert_270_mate1.fq','RH_S003_insert_270_mate2.fq'],
    # '4':['RH_S004_insert_270_mate1.fq','RH_S004_insert_270_mate2.fq'],'5':['RH_S005_insert_270_mate1.fq','RH_S005_insert_270_mate2.fq']}
    datasets={} ### let it empty if you do not have these sequences, e.g. datasets={}
    lr=[] ### ONT and Pb dataset. Let it empty if you do not have these sequences, e.g. lr=[]
    hifi_list=['SRR15275210.fastq','SRR15275211.fastq','SRR15275212.fastq','SRR15275213.fastq','SRR17687125.fastq']  #### Hifi dataset
    # datasets={'1':['RH_S001_insert_270_1.fastq','RH_S001_insert_270_2.fastq'], '2':['RH_S002_insert_270_1.fastq','RH_S002_insert_270_2.fastq'], '3':['RH_S003_insert_270_1.fastq','RH_S003_insert_270_2.fastq'],
    # '4':['RH_S004_insert_270_1.fastq','RH_S004_insert_270_2.fastq'],'5':['RH_S005_insert_270_1.fastq','RH_S005_insert_270_2.fastq']}
    # assembly_list=['1_RH_insert_cat_270_final.contigs.fa','2_RH_S001_insert_270_final.contigs.fa','3_RH_S002_insert_270_final.contigs.fa','4_RH_S003_insert_270_final.contigs.fa','5_RH_S004_insert_270_final.contigs.fa','6_RH_S005_insert_270_final.contigs.fa']
    assembly_list=['1_HumanGut_cat.fasta']
    pwd=os.getcwd()
    num_threads=60
    outlier_remover_main(bin_folder, depth_list, datasets, lr, hifi_list, assembly_list, pwd, num_threads)
