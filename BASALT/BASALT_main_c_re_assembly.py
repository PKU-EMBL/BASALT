#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Re-assembly entry point for the BASALT pipeline (CheckM-based branch).

This module orchestrates OLC-based reassembly, hybrid reassembly, and
final polishing of bins when long reads or HiFi datasets are provided.
"""

import time
import sys
import os
from Bio import SeqIO
from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import *
from S8_OLC_new_checkm import *
from S9_Reassembly_checkm import *
from S9p_Hybrid_Reassembly_checkm import *
from S10_OLC_new_checkm import *
from glob import glob
from Cleanup import *
from basalt_runtime import read_checkpoint_step, require_model_directory


def BASALT_main_c_re_assembly(assembly_list, datasets, num_threads, lr_list, hifi_list,
                              hic_list, eb_list, ram, continue_mode, functional_module,
                              sensitive, refinement_paramter, max_ctn, min_cpn,
                              pwd, QC_software, output_folder):
    """
    Run BASALT reassembly and polishing workflows.

    Parameters mirror those of the main pipeline and include information
    about assemblies, short reads, long reads/HiFi, quality thresholds,
    and runtime configuration.
    """
    ### Record the last accomplished step
    pwd=os.getcwd()

    require_model_directory()

    #### Program start
    last_step=0
    if continue_mode == 'last':
        last_step=read_checkpoint_step('Basalt_checkpoint.txt')
        if not os.path.exists('Basalt_checkpoint.txt'):
            f_cp_m=open('Basalt_checkpoint.txt', 'w')
            f_cp_m.close()
    else:
        print('Start a new project')
        cleanup(assembly_list)
        f_cp_m=open('Basalt_checkpoint.txt', 'w')
        f_cp_m.close()

    print('BASALT started from step: '+str(last_step))
    try:
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    except:
        fx=open('Basalt_log.txt','w')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    fx.close()
    
    try:
        f=open('BASALT_command.txt','a')
        f.write('BASALT restarted from step: '+str(last_step)+'\n')
        f.close()
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT restarted from step: '+str(last_step)+'\n')
        fx.close()
    except:
        f=open('BASALT_command.txt','w')
        f.write('BASALT started from 1st step'+'\n')
        f.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        f.write('Datasets: ')
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from 1st step'+'\n')
        fx.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        fx.write('Datasets: ')
        datasets_list=datasets.split(',')
        x=0
        for ds in datasets_list:
            x+=1
            ds1=str(ds).split('[')[1].replace('\'','').replace(']','').replace(' ','')
            if x == 1:
                f.write(str(ds1))
                fx.write(str(ds1))
            else:
                f.write('/'+str(ds1))
                fx.write('/'+str(ds1))
        f.write('\n')
        fx.write('\n')
        try:
            f.write('Long reads: '+str(str(lr_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
            fx.write('Long reads: '+str(str(lr_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        except:
            f.write('Long reads: []'+'\n')
            fx.write('Long reads: []'+'\n')
        f.write('Num threads: '+str(num_threads)+'\n'+'Ram: '+str(ram)+'\n'+'Refinement mode: '+str(refinement_paramter)+'\n')
        f.write('Autobinning mode: '+str(sensitive)+'\n'+'Functional mode: '+str(functional_module)+'\n'+'Continue mode: '+str(continue_mode)+'\n')
        f.write('Max contamination: '+str(max_ctn)+'\n'+'Min completeness: '+str(min_cpn)+'\n')
        f.close()
        fx.write('Num threads: '+str(num_threads)+'\n'+'Ram: '+str(ram)+'\n'+'Refinement mode: '+str(refinement_paramter)+'\n')
        fx.write('Autobinning mode: '+str(sensitive)+'\n'+'Functional mode: '+str(functional_module)+'\n'+'Continue mode: '+str(continue_mode)+'\n')
        fx.write('Max contamination: '+str(max_ctn)+'\n'+'Min completeness: '+str(min_cpn)+'\n')
        fx.close()

    x = 0
    datasets2=copy.deepcopy(datasets)
    for item in datasets.keys():
        hz_list=datasets[item][0].split('.')
        if len(hz_list) >= 2:
            if hz_list[-1] == 'fq' or hz_list[-1] == 'fastq':
                x=1
            elif hz_list[-1] == 'zip':
                x=2
            elif hz_list[-1] == 'gz':
                if hz_list[-2] == 'tar':
                    x=3
                else:
                    x=4
        else:
            print('Input format error! Please check the input file.')
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.'+'\n')
            fx.close()

    if x > 1:
        datasets={}
        for item in datasets2.keys():
            datasets[item]=[]
            if x == 2:
                f1_d=str(datasets2[item][0]).split('.zip')[0]
                f2_d=str(datasets2[item][1]).split('.zip')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('unzip '+str(datasets2[item][0]))

                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:
                    os.system('unzip '+str(datasets2[item][1]))

            elif x == 3:
                f1=str(datasets2[item][0])
                f2=str(datasets2[item][1])
                f1_d=str(f1).split('.tar.gz')[0]
                f2_d=str(f2).split('.tar.gz')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(f1))
                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(f2))

            elif x == 4:
                f1=str(datasets2[item][0])
                f2=str(datasets2[item][1])
                f1_d=f1.split('.gz')[0]
                f2_d=f2.split('.gz')[0]
                if os.path.exists(pwd+'/PE_r1_'+str(f1_d)):
                    z=0
                else:
                    os.system('gunzip -c '+f1+' > '+f1_d)
                if os.path.exists(pwd+'/PE_r2_'+str(f2_d)):
                    z=0
                else:                    
                    os.system('gunzip -c '+f2+' > '+f2_d)

            if '.fq' not in f1_d and '.fastq' not in f1_d:
                f1_d=f1_d+'.fq'
            if '.fq' not in f2_d and '.fastq' not in f2_d:
                f2_d=f2_d+'.fq'
            datasets[item].append(f1_d)
            datasets[item].append(f2_d)

    x=0
    assembly_list2=copy.deepcopy(assembly_list)
    for item in assembly_list:
        hz_list=assembly_list[0].split('.')
        if len(hz_list) >= 2:
            if hz_list[-1] == 'fa' or hz_list[-1] == 'fasta' or hz_list[-1] == 'fna':
                x=1
            elif hz_list[-1] == 'zip':
                x=2
            elif hz_list[-1] == 'gz':
                if hz_list[-2] == 'tar':
                    x=3
                else:
                    x=4
        else:
            print('Input format error! Please check the input file.')
            print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.')
            fx=open('Basalt_log.txt','a')
            fx.write('Input format error! Please check the input file.'+'\n')
            fx.write('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.'+'\n')
            fx.close()
    
    if x > 1:
        assembly_list=[]
        for item in assembly_list2:
            if x == 2:
                f_d=str(item).split('.zip')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:   
                    os.system('unzip '+str(item))
            
            elif x == 3:
                f_d=str(item).split('.tar.gz')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:
                    os.system('tar -zxf '+str(item))

            elif x == 4:
                f_d=str(item).split('.gz')[0]
                if os.path.exists(pwd+'/'+str(f_d)):
                    z=0
                else:
                    os.system('gunzip -c '+str(item)+' > '+str(f_d))

            if '.fa' not in f_d and '.fna' not in f_d and '.fasta' not in f_d:
                f_d=f_d+'.fa'
            assembly_list.append(f_d)

    ### Reassembly module
    if functional_module == 'reassembly' or functional_module == 'all':
        if last_step < 8:
            if len(assembly_list) > 1:
                print('Starting OLC process')
                if len(lr_list) == 0:
                    try:
                        bin_n=0
                        os.chdir(pwd+'/BestBinset_outlier_refined_filtrated_retrieved_retrieved')
                        for root, dirs, files in os.walk(pwd+'/BestBinset_outlier_refined_filtrated_retrieved_retrieved'):
                            for file in files:
                                hz=file.split('.')[-1]
                                if 'fa' in hz or 'fna' in hz:
                                    bin_n+=1
                        os.chdir(pwd)
                        if bin_n != 0:
                            target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_retrieved'
                        else:
                            target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'
                    except:
                        print('There is not bin in BestBinset_outlier_refined_filtrated_retrieved_retrieved. Try BestBinset_outlier_refined_filtrated_retrieved')
                        target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'

                    print('Processing with bins in '+str(target_bin_folder))
                    bin_comparison_folder='BestBinset_comparison_files'
                    step='assemblies_OLC'
                    aligned_len_cutoff=500
                    similarity_cutoff=99
                    coverage_extension=95
                    mod_bin_folder='' ### sometimes []
                    # OLC_main(target_bin_folder, step, bin_comparison_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, mod_bin_folder) ### OLC 07232023
                    OLC_main_checkm(target_bin_folder, step, bin_comparison_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, ram, mod_bin_folder)
                    # OLC_main(target_bin_folder, step, bin_comparison_or_reassembly_binset_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, ram, orig_binset)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'8th contig OLC done!'+'\t'+target_bin_folder+'_OLC')
                    f_cp_m.close()
                else:
                    for line in open('Basalt_checkpoint.txt', 'r'):
                        if '7th' == str(line[0:3]):
                            target_bin_folder=str(line).strip().split('\t')[1].strip()
                    # target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'
                    print('Long-read presented. Skip OLC process.')
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'8th contig OLC did not perform!'+'\t'+target_bin_folder)
                    f_cp_m.close()
            else:
                for line in open('Basalt_checkpoint.txt', 'r'):
                    if '7th' == str(line[0:3]):
                        target_bin_folder=str(line).strip().split('\t')[1].strip()

                print('Long-read presented. Skip OLC process.')
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'8th contig OLC did not perform!'+'\t'+target_bin_folder)
                f_cp_m.close()

        if last_step < 9:
            print('Starting reassembly')
            for line in open('Basalt_checkpoint.txt', 'r'):
                if '8th contig OLC' in line:
                    binset_folder=str(line).strip().split('\t')[1].strip()
                    print('Processing bins in '+str(binset_folder))

            datasets_list={}
            for ds in datasets.keys():
                datasets_list[ds]=[]
                datasets_list[ds].append('PE_r1_'+str(datasets[ds][0]))
                datasets_list[ds].append('PE_r2_'+str(datasets[ds][1]))
            
            if len(lr_list) == 0:
                hybri_reassembly='n'
                re_assembly_main(binset_folder, datasets_list, lr_list, hybri_reassembly, ram, num_threads)
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'9th reassembly done.'+'\t'+str(binset_folder)+'_re-assembly')
                f_cp_m.close()
            else:
                x=0
                lr_list2=copy.deepcopy(lr_list)
                for item in lr_list:
                    hz_list=lr_list[0].split('.')
                    if len(hz_list) >= 2:
                        if hz_list[-1] == 'fq' or hz_list[-1] == 'fastq':
                            x=1
                        elif hz_list[-1] == 'zip':
                            x=2
                        elif hz_list[-1] == 'gz':
                            if hz_list[-2] == 'tar':
                                x=3
                            else:
                                x=4
                    else:
                        print('Input format error! Please check the input file.')
                        print('BASALT supports the input (1) sequence files in  .gz, .zip, and .tar.gz; (2) and assemlies in .fa, .fna, .fasta.')
                
                if x > 1:
                    lr_list=[]
                    for item in lr_list2:
                        if x == 2:
                            f_d=str(item).split('.zip')[0]
                            if os.path.exists(pwd+'/'+str(f_d)):
                                z=0
                            else:   
                                os.system('unzip '+str(item))
                        
                        elif x == 3:
                            f_d=str(item).split('.tar.gz')[0]
                            if os.path.exists(pwd+'/'+str(f_d)):
                                z=0
                            else:
                                os.system('tar -zxf '+str(item))

                        elif x == 4:
                            f_d=str(item).split('.gz')[0]
                            if os.path.exists(pwd+'/'+str(f_d)):
                                z=0
                            else:
                                os.system('gunzip -c '+str(item)+' > '+str(f_d))
                        
                        if '.fq' not in f_d and '.fastq' not in f_d:
                            f_d=f_d+'.fq'
                        lr_list.append(f_d)
                polished_binset=binset_folder
                sr_folder='BestBinset_sr_bins_seq'
                lr_folder='BestBinset_long_read'
                # polished_binset='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished'
                # sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
                # lr_folder='BestBinset_outlier_refined_filtrated_retrieved_long_read'
                hybrid_re_assembly_main(polished_binset, sr_folder, lr_folder, ram, num_threads)
                os.system('rm -rf SPAdes_corrected_reads')
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'9th reassembly done.'+'\t'+str(polished_binset)+'_re-assembly')
                f_cp_m.close()
        
        if last_step < 10:
            print('Starting reassemblies OLC')
            for line in open('Basalt_checkpoint.txt', 'r'):
                if '9th reassembly done' in line:
                    target_bin_folder=str(line).strip().split('\t')[1].strip()
                    print('Processing bins in '+str(target_bin_folder))

            bin_comparison_or_reassembly_binset_folder=target_bin_folder+'_binset'
            orig_binset=target_bin_folder.split('_re-assembly')[0]+'_mod'
            step='OLC_after_reassembly'
            aligned_len_cutoff=500
            similarity_cutoff=98
            coverage_extension=90

            reassembly_OLC_main(target_bin_folder, step, bin_comparison_or_reassembly_binset_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, ram, orig_binset)
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'10th reassemblies OLC done.'+'\t'+target_bin_folder+'_OLC')
            f_cp_m.close()
        
        try:
            if last_step < 11:
                if len(lr_list) != 0:
                    ### Long read present. Start polishing
                    print('Starting 2nd run polishing')
                    target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC'
                    remained_seq_list=glob(r'Remained_seq1.fq_*')
                    i=0
                    for item in remained_seq_list:
                        num=int(item.split('.fq_')[1])
                        if num > i:
                            i = num
                    os.system('mv Remained_seq1.fq_'+str(num)+' 2nd_Remained_seq1.fq')
                    os.system('mv Remained_seq2.fq_'+str(num)+' 2nd_Remained_seq2.fq')
                    unmapped_datasets={'1':['2nd_Remained_seq1.fq','2nd_Remained_seq2.fq']}
                    # unmapped_datasets={'1':['Remained_seq1.fq','Remained_seq2.fq']}
                    polishing_main(target_bin_folder, unmapped_datasets, assembly_list, lr_list, 2, 'bwa', 3, ram, num_threads)
                    remained_seq_list=glob(r'Remained_seq1.fq_*')
                    i=0
                    for item in remained_seq_list:
                        num=int(item.split('.fq_')[1])
                        if num > i:
                            i = num
                    
                    for i in range(1,num+1):
                        os.system('rm Remained_seq1.fq_'+str(i)+' '+'Remained_seq2.fq_'+str(i))
                    os.system('mkdir Remained_seq')
                    os.system('mv 2nd_Remained_seq1.fq '+pwd+'/Remained_seq/Remained_seq1.fq')
                    os.system('mv 2nd_Remained_seq2.fq '+pwd+'/Remained_seq/Remained_seq2.fq')
                    os.system('tar -zcvf Remained_seq.tar.gz Remained_seq')
                    os.system('rm -rf Remained_seq')
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'11th 2nd Polishing done!')
                    f_cp_m.close()
                else:
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'11th did not perform polish')
                    f_cp_m.close()
        except:
            print('Did not perform the 2nd polish. Skit the step.')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'11th 2nd polishing did not perform')
            f_cp_m.close()

        if last_step < 12:
            ### Final de-replication: S4 function
            for line in open('Basalt_checkpoint.txt', 'r'):
                try:
                    if '10th reassemblies OLC done' in line:
                        final_folder=str(line).strip().split('\t')[1].strip()
                except:
                    xyz=0

                try:
                    if '11th 2nd Polishing done' in line:
                        final_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_MAGs_polished'
                except:
                    xyz=0
                
                try:
                    if '11th 2nd did not perform' in line:
                        final_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC'
                except:
                    xyz=0

            if len(assembly_list) != 1:                
                coverage_matrix_list, bestbinset_list = [], []
                # for line in open('Bestbinset_list.txt','r'):
                #     bestbinset_list.append(str(line).strip())
            
                # for line in open('Coverage_matrix_list.txt','r'):
                #     coverage_matrix_list.append(str(line).strip())

                # datasets_fq={}
                # for item in datasets.keys():
                #     datasets_fq[item]=[]
                #     datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
                #     datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))

                # multiple_assembly_comparitor_main(drep_list, bestbinset_list, coverage_matrix_list, datasets_fq, 'final_drep', num_threads)
                try:
                    final_binset_comparitor(final_folder, coverage_matrix_list, datasets, num_threads, pwd, 'final_drep')
                except:
                    print('Final drep did not perform')
            #     # final_binset_comparitor(final_folder, coverage_matrix_list, datasets_fq, num_threads, pwd, 'final_drep')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'12th final de-replication done!')
            f_cp_m.close()

        if last_step < 13:
            ### Final binset
            os.system('mv '+str(final_folder)+' '+str(output_folder))
            rename={}
            os.chdir(pwd+'/'+str(output_folder))
            rename_file=open('Rename.txt', 'w')
            for root, dirs, files in os.walk(pwd+'/'+str(output_folder)):
                for file in files:
                    hz=file.split('.')[-1]
                    if 'fa' in hz or 'fna' in hz or 'fasta' in hz:
                        if '_' in file:
                            nl=file.split('.')
                            nl.remove(nl[-1])
                            full_name='.'.join(nl)
                            bin_name=file.split('_')[0].split('.')[0]
                            os.system('mv '+file+' '+bin_name+'.fa')
                            rename[full_name]=bin_name
                            rename_file.write(file+' '+bin_name+'.fa'+'\n')
            rename_file.close()
            
            for root, dirs, files in os.walk(pwd+'/'+str(output_folder)):
                for file in files:
                    if 'bin_stats_ext.tsv' in file:
                        f=open(file+'2','w')
                        for line in open(file):
                            n=str(line).strip().split('\t')[0]
                            c=str(line).strip().split('\t')[1]
                            if n in rename.keys():
                                f.write(str(rename[n])+'\t'+str(c)+'\n')
                            else:
                                f.write(line)
                        f.close()
                        
                        os.system('mv '+file+'2 '+file)
            os.system('rm *bin_stats_ext_o.tsv')
            os.chdir(pwd)

            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'13th BASALT binning done!')
            f_cp_m.close()
    print('BASALT main program accomplished!')
    print('BASALT will continue to cleanup or compress all the temp files. Please wait for a little bit longer')

    ### Cleanup
    # if functional_module == 'all':
    #     os.system('rm *.njs *.ndb *.nto *.ntf *.not *.nos')
    #     os.mkdir('Coverage_depth_connection_SimilarBin_files_backup')
    #     os.system('mv *.depth.txt Coverage_matrix_* Combat_* condense_connections_* Connections_* Similar_bins.txt Coverage_depth_connection_SimilarBin_files_backup')
    #     os.system('tar -zcvf Coverage_depth_connection_SimilarBin_files_backup.tar.gz Coverage_depth_connection_SimilarBin_files_backup')
    #     os.system('rm -rf Coverage_depth_connection_SimilarBin_files_backup')
    #     os.system('rm -rf *_kmer bin_coverage Bin_coverage_after_contamination_removal bin_comparison_folder bin_extract-eleminated-selected_contig Bins_blast_output')
    #     os.system('tar -zcvf Group_comparison_files.tar.gz *_comparison_files')
    #     os.system('tar -zcvf Group_Bestbinset.tar.gz *_BestBinsSet')
    #     # os.system('tar -zcvf Group_checkm.tar.gz *_checkm')
    #     os.system('tar -zcvf Group_genomes.tar.gz *_genomes')
    #     os.system('rm -rf *_sr_bins_seq')
    #     os.system('tar -zcvf Binsets_backup.tar.gz BestBinse*') ###
    #     os.system('rm -rf *_comparison_files *_checkm *_genomes *_BestBinset *_BestBinsSet BestBinse* Deep_retrieved_bins coverage_deep_refined_bins S6_coverage_filtration_matrix S6_TNF_filtration_matrix split_blast_output TNFs_deep_refined_bins')
    #     # os.system('rm *_checkpoint.txt')
    #     os.system('rm -rf Merged_seqs_*')
    #     os.system('rm -rf *.bt2 Outlier_in_threshold* Summary_threshold* Refined_total_bins_contigs.fa Total_bins.fa') 
    #     for i in range(1,20):
    #         os.system('rm -rf *_deep_retrieval_'+str(i))
    #     os.system('rm -rf *_MP_1 *_MP_2 *_gf_lr_polished *_gf_lr *_gf_lr_mod *_gf_lr_checkm *_long_read')
    #     os.system('rm Bin_reads_summary.txt Depth_total.txt Basalt_log.txt Assembly_mo_list.txt Assembly_MoDict.txt *_gf_lr_blasted.txt Bestbinset_list.txt Bin_extract_contigs_after_coverage_filtration.txt Bin_lw.txt Bin_record_error.txt')
    #     os.system('rm Bins_folder.txt BLAST_output_error.txt Concoct_* condensed.cytoscape* cytoscape.*')
    #     os.system('rm Hybrid_re-assembly_status.txt Mapping_log_* OLC_merged_error_blast_results.txt Potential_contaminted_seq_vari.txt Reassembled_bins_comparison.txt Rejudge_clean.txt')
    #     os.system('rm Remained_seq* Remapped_depth_test.txt Re-mapped_depth.txt Remapping.fasta TNFs_exceptional_contigs.txt Total_contigs_after_OLC_reassembly.fa')
    #     # os.system('rm PE_r1_* PE_r2_*')
    #     for i in range(1, len(assembly_list)+1):
    #         os.system('rm '+str(i)+'_'+str(assembly_list[i-1]))
    # print('Done!')

if __name__ == '__main__': 
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    datasets={'1':['RM2_S001_insert_270_mate1.fq','RM2_S001_insert_270_mate2.fq'], '2':['RM2_S002_insert_270_mate1.fq','RM2_S002_insert_270_mate2.fq']}
    lr_list=[]
    hic_list=[]
    hifi_list=[]
    eb_list=[] ### if you want to use extra binner, you may use eb_list=['m', 'v']; 'm': metabinner, 'v': vamb
    num_threads=30
    ram=120
    pwd=os.getcwd()
    refinement_paramter='deep' ### 1st refinement mode: 1. quick; 2. deep
    # autobining_parameters='more-sensitive' ### 'quick', 'sensitive', 'more-sensitive'
    functional_module='all' ### 1. all; 2. autobinning; 3. refinement; 4. reassembly
    continue_mode='last' ### parameter type: 1. continue(last): start from the latest step finished from the last time; 2. new: start from beginning
    # hybri_reassembly='n' ### Use Unicycler to re-assembly. e.g. --hybrid y / --hybrid n; defalt no
    max_ctn, min_cpn=20, 35 
    QC_software='checkm2'
    output_folder='Final_binset'
    sensitive='sensitive' ### defaul sensitive; option: quick, more-sensitive
    BASALT_main_c_re_assembly(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitive, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
    # BASALT_main(assembly_list, datasets, num_threads, lr_list, ram, continue_mode, functional_module, autobining_parameters, refinement_paramter, hybri_reassembly, max_ctn, min_cpn, pwd)
