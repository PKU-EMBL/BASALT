#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Autobinning entry point for the BASALT pipeline (CheckM-based branch).

This module runs the initial binning workflows, including multiple
autobinners and optional extra binners, and records progress in
``Basalt_checkpoint.txt`` to support resumable execution.
"""

import time
import sys
import os
from Bio import SeqIO
from S1_Autobinners_2qc_11152023 import *
from S1e_extra_binners import *
from S2_BinsAbundance_PE_connections_multiple_processes_pool_checkm import *
from S3_Bins_comparator_within_group_checkm import *
from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import *
from glob import glob
from Cleanup import *
from basalt_runtime import read_checkpoint_step, require_model_directory


def BASALT_main_c_autobinning(assembly_list, datasets, num_threads, lr_list, hifi_list,
                              hic_list, eb_list, ram, continue_mode, functional_module,
                              sensitive, refinement_paramter, max_ctn, min_cpn,
                              pwd, QC_software, output_folder):
    """
    Run the BASALT autobinning module using CheckM-based quality control.

    Parameters
    ----------
    assembly_list : list
        List of assembly fasta file paths.
    datasets : dict
        Mapping from dataset id to paired-end read paths.
    num_threads : int
        Number of CPU threads to use.
    lr_list : list
        Long read datasets (ONT/PB).
    hifi_list : list
        HiFi read datasets.
    hic_list : list
        Hi-C datasets (if available).
    eb_list : list
        Extra binner identifiers to use in addition to the default set.
    ram : int
        Available RAM in gigabytes.
    continue_mode : str
        Either ``'last'`` to resume from checkpoints or ``'new'``.
    functional_module : str
        Functional mode flag from the top-level CLI (e.g. 'all', 'autobinning').
    sensitive : str
        Autobinning sensitivity level ('quick', 'sensitive', 'more-sensitive').
    refinement_paramter : str
        Refinement mode ('quick' or 'deep').
    max_ctn : int
        Maximum allowed contamination for bins.
    min_cpn : int
        Minimum required completeness for bins.
    pwd : str
        Working directory.
    QC_software : str
        Quality check software identifier, expected 'checkm' or 'checkm2'.
    output_folder : str
        Output folder name for final binsets.
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

    ### Autobinner
    if functional_module == 'autobinning' or functional_module == 'all':
        if last_step == 0:
            print('----------------------------------')
            print('Running autobinner')
            print('Processing files in '+str(pwd))
            fx=open('Basalt_log.txt','a')
            fx.write('Running autobinner'+'\n')
            fx.write('Processing files in '+str(pwd)+'\n')
            fx.close()
    
            insert_size=100
            A=autobinner_main(assembly_list, datasets, lr_list, hifi_list, insert_size, num_threads, ram, sensitive, QC_software, pwd)
            bins_folders_dic=A[0]
            connections_total_dict=A[1]
            depth_total=A[2]
            assembly_MoDict=A[3]

            if len(eb_list) != 0:
                for binner in eb_list:
                    for i in range(1,len(assembly_list)+1):
                        assembly_file=str(i)+'_'+str(assembly_list[i-1])
                        depth_file=str(i)+'_assembly.depth.txt'

                        try:
                            extra_bin_folder=extra_binner(binner, datasets, assembly_file, depth_file, num_threads, ram, pwd, QC_software)
                            for item in extra_bin_folder:
                                bins_folders_dic[str(assembly_list[i-1])].append(item)
                            print(f"[SUCCESS] Extra binner '{binner}' completed for {assembly_file}")
                            fx=open('Basalt_log.txt','a')
                            fx.write(f"[SUCCESS] Extra binner '{binner}' completed for {assembly_file}\n")
                            fx.close()
                        except Exception as e:
                            print(f"[WARNING] Extra binner '{binner}' failed for {assembly_file}. Skipping...")
                            print(f"[ERROR DETAILS] {str(e)}")
                            try:
                                fx=open('Basalt_log.txt','a')
                                fx.write(f"[WARNING] Extra binner '{binner}' failed for {assembly_file}: {str(e)}\n")
                                fx.close()
                            except:
                                pass
                            continue

            f1=open('Connections_total_dict.txt','w')
            for item in connections_total_dict.keys():
                f1.write(str(item)+'\t'+str(connections_total_dict[item])+'\n')
            f1.close()

            f1=open('Depth_total.txt','w')
            for item in depth_total.keys():
                f1.write(str(item)+'\t'+str(depth_total[item])+'\n')
            f1.close()

            f1=open('Assembly_MoDict.txt','w')
            for item in assembly_MoDict.keys():
                f1.write(str(item)+'\t'+str(assembly_MoDict[item])+'\n')
            f1.close()

            # bins_folders, n={}, 0
            f1=open('Bins_folder.txt','w')
            for item in bins_folders_dic.keys():
                # n+=1
                # bins_folders[str(n)]=[]
                f1.write(str(item)+'\t'+str(bins_folders_dic[item])+'\n')
                # genomes_list=str(bins_folders_dic[item]).strip().replace('[','').replace(']','').replace('\'','').replace(' ','').split(',')
                # for genomes_folder in genomes_list:
                #     genomes_folder_name_list=genomes_folder.split('_')
                #     genomes_folder_name_list.remove(genomes_folder_name_list[-1])
                #     genomes_folder_name='_'.join(genomes_folder_name_list)
                #     bins_folders[str(n)].append(genomes_folder_name)
            f1.close()    
            os.system('rm *.bam')
            
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('1st autobinner done!')
            f_cp_m.close()
            fx=open('Basalt_log.txt','a')
            fx.write('1st autobinner done!'+'\n')
            fx.close()
            #connections_total_dict={'9groups_assembly.fa':'condense_connections_9groups_assembly.fa.txt', 'AN1-32_assembly.fa':'condense_connections_AN1-32_assembly.fa.txt'}
            #depth_total={'9groups_assembly.fa':'1_assembly.depth.txt', 'AN1-32_assembly.fa':'2_assembly.depth.txt'}
            #assembly_MoDict={'9groups_assembly.fa':'1_9groups_assembly.fa', 'AN1-32_assembly.fa':'2_AN1-32_assembly.fa'}
            #bins_folders={'9groups_assembly.fa':['1_9groups_assembly.fa_0.3_genomes','1_9groups_assembly.fa_0.5_genomes','1_9groups_assembly.fa_0.7_genomes','1_9groups_assembly.fa_0.9_genomes','1_9groups_assembly.fa_1000_concoct_genomes','1_9groups_assembly.fa_400_concoct_genomes','1_9groups_assembly.fa_200_concoct_genomes','1_9groups_assembly.fa_200_metabat_genomes','1_9groups_assembly.fa_300_metabat_genomes','1_9groups_assembly.fa_400_metabat_genomes','1_9groups_assembly.fa_500_metabat_genomes'], 'AN1-32_assembly.fa':['2_AN1-32_assembly.fa_0.3_maxbin2_genomes','2_AN1-32_assembly.fa_0.9_maxbin2_genomes','2_AN1-32_assembly.fa_0.7_maxbin2_genomes','2_AN1-32_assembly.fa_0.5_maxbin2_genomes','2_AN1-32_assembly.fa_200_metabat_genomes','2_AN1-32_assembly.fa_300_metabat_genomes','2_AN1-32_assembly.fa_400_metabat_genomes','2_AN1-32_assembly.fa_500_metabat_genomes','2_AN1-32_assembly.fa_200_concoct_genomes','2_AN1-32_assembly.fa_400_concoct_genomes','2_AN1-32_assembly.fa_1000_concoct_genomes']}

        ### Bin selecting process
        if last_step < 2:
            connections_total_dict, depth_total, assembly_MoDict, bins_folders = {}, {}, {}, {}
            n=0
            for line in open('Connections_total_dict.txt','r'):
                n+=1
                # connections_total_dict[str(line).strip().split('\t')[1].strip().split('_')[0]]=str(line).strip().split('\t')[1].strip()
                connections_total_dict[str(n)]=str(line).strip().split('\t')[1].strip()
            n=0
            for line in open('Depth_total.txt','r'):
                n+=1
                depth_total[str(n)]=str(line).strip().split('\t')[1].strip()
                # depth_total[str(line).strip().split('\t')[1].strip().split('_')[0]]=str(line).strip().split('\t')[1].strip()
            n=0
            for line in open('Assembly_MoDict.txt','r'):
                n+=1
                assembly_MoDict[str(n)]=str(line).strip().split('\t')[1].strip()
                # assembly_MoDict[str(line).strip().split('\t')[1].strip().split('_')[0]]=str(line).strip().split('\t')[1].strip()

            n=0
            for line in open('Bins_folder.txt','r'):
                n+=1
                # assembly_name=str(line).strip().split('[')[1].replace('\'','').split('_')[0]
                bins_folders[str(n)]=[]
                genomes_list=str(line).strip().split('\t')[1].strip().replace('[','').replace(']','').replace('\'','').replace(' ','').split(',')
                for genomes_folder in genomes_list:
                    genomes_folder_name_list=genomes_folder.split('_')
                    genomes_folder_name_list.remove(genomes_folder_name_list[-1])
                    genomes_folder_name='_'.join(genomes_folder_name_list)
                    bins_folders[str(n)].append(genomes_folder_name)
                    
            Ax=binsabundance_pe_connections(bins_folders, depth_total, connections_total_dict, assembly_MoDict, num_threads)

            fx=open('Basalt_log.txt','a')
            fx.write('Starting best bin selection: within each group'+'\n')
            print('Starting best bin selection: within each group')
            coverage_matrix_list, bestbinset_list, assembly_mo_list=[], [], []
            for item in bins_folders.keys():
                print('---------------------------')
                print('Processing group '+str(item)+' '+str(bins_folders[item]))
                fx.write('---------------------------'+'\n')
                fx.write('Processing group '+str(item)+' '+str(bins_folders[item])+'\n')
                coverage_matrix_list.append(Ax[item])
                genome_folder_list=bins_folders[item]
                genomes_folder_name_list=[]
                for item2 in genome_folder_list:
                    genomes_folder_name_list.append(item2+'_genomes')
                # Bx=bins_comparator_multiple_groups(genome_folder_list, assembly_MoDict[item])
                Bx=bins_comparator_multiple_groups(genomes_folder_name_list, assembly_MoDict[item])
                print('Adding '+str(Bx)+' to the bestbinset list')
                fx.write('Adding '+str(Bx)+' to the bestbinset list'+'\n')
                bestbinset_list.append(Bx)
                print('Adding '+str(assembly_MoDict[item])+' to the assembly modified list')
                fx.write('Adding '+str(assembly_MoDict[item])+' to the assembly modified list'+'\n')
                assembly_mo_list.append(assembly_MoDict[item])

            f1=open('Coverage_matrix_list.txt','w')
            for item in coverage_matrix_list:
                f1.write(str(item)+'\n')
            f1.close()

            f1=open('Bestbinset_list.txt','w')
            for item in bestbinset_list:
                f1.write(str(item)+'\n')
            f1.close()

            f1=open('Assembly_mo_list.txt','w')
            for item in assembly_mo_list:
                f1.write(str(item)+'\n')
            f1.close()

            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'2nd bin selection within group done!')
            f_cp_m.close()

            #coverage_matrix_list=['Coverage_matrix_for_binning_1_I4_assembly.fa.txt','Coverage_matrix_for_binning_2_I6_assembly.fa.txt','Coverage_matrix_for_binning_3_I12_assembly.fa.txt','Coverage_matrix_for_binning_4_I_assembly.fa.txt']
            #bestbinset_list=['1_I4_assembly.fa_BestBinsSet','2_I6_assembly.fa_BestBinsSet','3_I12_assembly.fa_BestBinsSet','4_I_assembly.fa_BestBinsSet']
            #assembly_mo_list=['1_I4_assembly.fa', '2_I6_assembly.fa', '3_I12_assembly.fa', '4_I_assembly.fa']

        if last_step < 3:
            if len(assembly_list) != 1:
                print('Starting best bin selection: within multiple groups')
                if last_step == 2:
                    coverage_matrix_list, bestbinset_list, assembly_mo_list=[], [], []
                    for line in open('Coverage_matrix_list.txt','r'):
                        coverage_matrix_list.append(str(line).strip())

                    for line in open('Bestbinset_list.txt','r'):
                        bestbinset_list.append(str(line).strip())
                    
                    for line in open('Assembly_mo_list.txt','r'):
                        assembly_mo_list.append(str(line).strip())

                # multiple_assembly_comparitor_main(assembly_mo_list, bestbinset_list, coverage_matrix_list, num_threads)
                datasets_fq={}
                for item in datasets.keys():
                    datasets_fq[item]=[]
                    datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
                    datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))
                multiple_assembly_comparitor_main(assembly_mo_list, bestbinset_list, coverage_matrix_list, datasets_fq, 'initial_drep', num_threads)

                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'3rd bin selection within multiple groups done!')
                f_cp_m.close()
            else:
                os.system('cp -r '+str(assembly_list[0])+'_comparison_files BestBinset_comparison_files')
                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'3rd bin selection did not perform, because there is only one assembly!')
                f_cp_m.close()
        print('Autobinning accomplished')

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
    BASALT_main_c_autobinning(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitive, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
    # BASALT_main(assembly_list, datasets, num_threads, lr_list, ram, continue_mode, functional_module, autobining_parameters, refinement_paramter, hybri_reassembly, max_ctn, min_cpn, pwd)
