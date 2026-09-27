#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S5 (CheckM branch): Deep-learning-based outlier bin removal.

CheckM-specific implementation of DL-based outlier removal using
CheckM quality reports in combination with learned models.
"""

from Bio import SeqIO
import os, copy, math, glob, gc
import numpy as np
from multiprocessing import Pool
from basalt_runtime import depth_sample_count

def TNF_coverage_matrix(bin_contigs, bin_id, contigs_depth, ccc, contigs_kmer, contigs_kmer2):
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

def bin_kmer(kmer, contigs_bin, contigs_bin2, bin_contigs, contigs_depth, ccc):
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
    # print(str(bin_contigs))

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
    # print(str(kmer_list))
    
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

def outlier_predictor(depth_TNF_matrix, contigs_depth, bin_contigs, datasets, num_threads, nx):
    print('Predicting potential outliers')
    os.system('ensemble.py')
    # outlier_predictor.py
    contig_status, bin_positive_contigs, bin_avg_coverage, bin_core_contig_num = {}, {}, {}, {}
    for line in open('Predicted_potential_outlier.txt','r'):
        num=str(line).strip().count('\t')
        if num == 0:
            bin_id=str(line).strip()
            bin_core_contig_num[bin_id]=0
            # contig_status[bin_id]={}
            # contig_status[bin_id]['Real']=[]
            contig_status[bin_id]=[]
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
            elif status == 'Contaminated':
                contig_status[bin_id].append(contig_id)
    
    bin_positive_contigs_IQR = {}
    # bin_positive_contigs_IQR, bin_avg_coverage_media = {}, {}
    for bin_id in bin_positive_contigs.keys():
        bin_positive_contigs_IQR[bin_id]={}
        bin_avg_coverage[bin_id]={}
        # bin_avg_coverage_media[bin_id]={}
        for cov in bin_positive_contigs[bin_id].keys():
            if len(bin_positive_contigs[bin_id][cov]) >= 1:
                bin_positive_contigs_IQR[bin_id][cov]={}
                a=np.array(bin_positive_contigs[bin_id][cov])
                Q1=np.percentile(a,25)
                Q3=np.percentile(a,75)
                # Q2=np.percentile(a,50)
                avg=np.mean(a)
                IQR=Q3-Q1
                upper = Q3 + IQR
                lower = Q1 - IQR
                bin_positive_contigs_IQR[bin_id][cov]['lower']=lower
                bin_positive_contigs_IQR[bin_id][cov]['upper']=upper
                bin_avg_coverage[bin_id][cov]=avg
                # bin_avg_coverage_media[bin_id][cov]=Q2
    
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
                        # f.write('Contaminated'+'\t'+tt)
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
    # os.system('bwa index '+str(total_fa))
    os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))

    n = 0
    for item in datasets.keys():
        print('Mapping '+str(datasets[item]))
        n+=1
        # os.system('bwa mem -t '+str(num_threads)+' '+str(total_fa)+' '+str(datasets[item][0])+' '+str(datasets[item][1])+' > '+str(item)+'.sam') ### bwa
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets[item][0])+' -2 '+str(datasets[item][1])+' -S '+str(item)+'.sam -q --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(item)+'.sam -o '+str(item)+'.bam')
        ## py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(item)+'.bam '+str(item)+'_sorted') 

        try:
            with open(str(item)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(item)+'.bam failed. Redoing')
            ## py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(item)+'_sorted.bam '+str(item)+'.bam')
        os.system('rm '+str(item)+'.sam')

        if item == '1':
            bam_sorted='1_sorted.bam'
        else:
            bam_sorted+=' '+str(item)+'_sorted.bam'
    
    try:
        os.system('jgi_summarize_bam_contig_depths --outputDepth Re-mapped_depth.txt '+str(bam_sorted))
        nxxyy=0
        for line in open('Re-mapped_depth.txt','r'):
            nxxyy+=1
            if nxxyy == 2:
                break
    except:
        os.system('jgi_summarize_bam_contig_depths --outputDepth Re-mapped_depth.txt '+str(bam_sorted))
    # os.system('rm  *.bt2')
    os.system('rm '+str(bam_sorted)+' *.bt2')
    
    for item in datasets.keys():
        os.system('rm '+str(item)+'.bam')

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
    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(outlier_binset)+' '+str(outlier_binset)+'_checkm')
    os.chdir(str(outlier_binset)+'_checkm/storage/')
    print('Parsing '+outlier_binset+' checkm output')
    refined_checkm={}
    for line in open('bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].split('}')[0].strip()
        # GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
        refined_checkm[str(binID)]={}
        refined_checkm[str(binID)]['marker lineage']=str(taxon)
        refined_checkm[str(binID)]['Completeness']=float(completeness)
        refined_checkm[str(binID)]['Genome size']=float(eval(genome_size))
        refined_checkm[str(binID)]['Contamination']=float(contamination)
    
    os.chdir(pwd+'/'+str(bin_folder))
    orig_checkm={}
    for root, dirs, files in os.walk(pwd+'/'+str(bin_folder)):
        for file in files:
            if '.tsv' in file:
                for line in open(file,'r'):
                    binID=str(line).strip().split('{\'')[0].strip()
                    if binID in refined_checkm.keys():
                        genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
                        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
                        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].split('}')[0].strip()
                        qua=float(completeness)-5*float(contamination)

                        refined_qua=refined_checkm[str(binID)]['Completeness']-5*refined_checkm[str(binID)]['Contamination']
                        if qua > refined_qua:
                            os.system('rm '+pwd+'/'+str(outlier_binset)+'/'+binID+'.fa')
                            os.system('rm '+pwd+'/'+str(outlier_binset)+'/'+binID+'.fasta')
                            os.system('cp '+binID+'.fa '+pwd+'/'+str(outlier_binset))
                            os.system('cp '+binID+'.fasta '+pwd+'/'+str(outlier_binset))
                            refined_checkm[str(binID)]['marker lineage']=str(taxon)
                            refined_checkm[str(binID)]['Completeness']=float(completeness)
                            refined_checkm[str(binID)]['Genome size']=float(eval(genome_size))
                            refined_checkm[str(binID)]['Contamination']=float(contamination)
                    else:
                        os.system('cp '+binID+'.fa '+pwd+'/'+str(outlier_binset))
                        os.system('cp '+binID+'.fasta '+pwd+'/'+str(outlier_binset))
                        refined_checkm[str(binID)]={}
                        refined_checkm[str(binID)]['marker lineage']=str(taxon)
                        refined_checkm[str(binID)]['Completeness']=float(completeness)
                        refined_checkm[str(binID)]['Genome size']=float(eval(genome_size))
                        refined_checkm[str(binID)]['Contamination']=float(contamination)
    
    os.chdir(pwd+'/'+str(outlier_binset))
    f=open('bin_stats_ext.tsv','w')
    for binID in refined_checkm.keys():
        f.write(str(binID)+'\t'+str(refined_checkm[str(binID)])+'\n')
    f.close()
    os.chdir(pwd)

def kmer_cal(input_file, output_file):
    """
    Calculate k-mer features for the given FASTA and write to output file.
    """
    os.system('calc.kmerfreq.pl -i '+str(input_file)+' -o '+str(output_file))

def outlier_remover_main(bin_folder, depth_list, datasets, assemblies_list,
                         pwd, num_threads):
    """
    Entry point for S5 (CheckM branch): DL-based outlier removal using CheckM.

    Parameters
    ----------
    bin_folder : str
        Folder containing original BestBinset bins.
    depth_list : list of str
        List of coverage matrix files used for feature extraction.
    datasets : dict
        Mapping dataset_id -> [R1, R2] short-read files.
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
    A=outlier_predictor(depth_TNF_matrix, contigs_depth, bin_contigs, datasets, num_threads, nx)
    contig_status, bin_avg_coverage, confirmed_outlier= A[0], A[1], A[2]
    checkm_eval(bin_contigs, bin_folder, confirmed_outlier, pwd, num_threads)
    del contigs_depth, ccc, bin_contigs, contig_status, bin_avg_coverage, confirmed_outlier ### Release ram
    gc.collect()
    print('Outlier removal done!')

if __name__ == '__main__': 
    bin_folder='BestBinset'
    # depth_list=['Coverage_matrix_for_binning_1_Opera_unpolished_cat_contigs.fasta.txt']
    depth_list=['Coverage_matrix_for_binning_1_assembly.fa.txt','Coverage_matrix_for_binning_4_assembly.fa.txt']
    # datasets={'1':['RH_S001_insert_270_mate1.fq','RH_S001_insert_270_mate2.fq'], '2':['RH_S002_insert_270_mate1.fq','RH_S002_insert_270_mate2.fq'], '3':['RH_S003_insert_270_mate1.fq','RH_S003_insert_270_mate2.fq'],
    # '4':['RH_S004_insert_270_mate1.fq','RH_S004_insert_270_mate2.fq'],'5':['RH_S005_insert_270_mate1.fq','RH_S005_insert_270_mate2.fq']}
    # datasets={'1':['RH_S001_insert_270_1.fastq','RH_S001_insert_270_2.fastq'], '2':['RH_S002_insert_270_1.fastq','RH_S002_insert_270_2.fastq'], '3':['RH_S003_insert_270_1.fastq','RH_S003_insert_270_2.fastq'],
    # '4':['RH_S004_insert_270_1.fastq','RH_S004_insert_270_2.fastq'],'5':['RH_S005_insert_270_1.fastq','RH_S005_insert_270_2.fastq']}
    datasets={'1':['PE_r1_SRR12358675_sub_1.fastq','PE_r2_SRR12358675_sub_2.fastq'],'2':['PE_r1_SRR12358676_sub_1.fastq','PE_r2_SRR12358676_sub_2.fastq']}
    assembly_list=['1_assembly.fa','4_assembly.fa']
    pwd=os.getcwd()
    num_threads=60
    outlier_remover_main(bin_folder, depth_list, datasets, assembly_list, pwd, num_threads)
