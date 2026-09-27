#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Refinement entry point for the BASALT pipeline (CheckM-based branch).

This module performs dereplication, DL-based outlier removal, contig
retrieval, and optional polishing to refine the quality of MAGs.
"""

import time
import sys
import os
from Bio import SeqIO
from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import *
from S5_Outlier_remover_DL_checkm import *
from S6_retrieve_contigs_from_PE_contigs_checkm import *
from S7_Contigs_retrieve_within_group_checkm import *
from S7lr_finding_sr_contigs_basing_lr_and_polishing_checkm import *
from glob import glob
from Cleanup import *
from basalt_runtime import pe_read_names, read_checkpoint_step, require_model_directory


def BASALT_main_c_refinement(assembly_list, datasets, num_threads, lr_list, hifi_list,
                             hic_list, eb_list, ram, continue_mode, functional_module,
                             sensitive, refinement_paramter, max_ctn, min_cpn,
                             pwd, QC_software, output_folder):
    """
    Run BASALT refinement workflows for an existing autobinning result.

    The function reads and updates checkpoint files, allowing interrupted
    refinement runs to be resumed.
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

    ### Refinement
    if functional_module == 'refinement' or functional_module == 'all':
        if last_step < 4:
            print('Starting outlier removal process')
            coverage_matrix_list, connections_list, assembly_mo_list, bestbinset_list = [], [], [], []
            for line in open('Coverage_matrix_list.txt','r'):
                coverage_matrix_list.append(str(line).strip())
            for line in open('Assembly_mo_list.txt','r'):
                assembly_mo_list.append(str(line).strip())
            for line in open('Bestbinset_list.txt','r'):
                bestbinset_list.append(str(line).strip())
            
            for item in assembly_list:
                connections_list.append('condense_connections_'+item+'.txt')
            if len(bestbinset_list) == 1:
                print('Copying '+str(bestbinset_list[0])+' to BestBinset')
                os.system('cp -r '+str(bestbinset_list[0])+' BestBinset')
            # contig_outlier_remover_main('BestBinset', coverage_matrix_list, connections_list, num_threads, ram)
            outlier_remover_main('BestBinset', coverage_matrix_list, datasets, assembly_mo_list, pwd, num_threads)
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'4th outlier removal done!')
            f_cp_m.close()

        if last_step < 5:
            print('Starting contig retrival process')
            best_binset_from_multi_assemblies='BestBinset_outlier_refined'
            outlier_remover_folder='BestBinset_outlier_refined'
            
            coverage_matrix_list, connections_list, assembly_mo_list = [], [], []
            for line in open('Coverage_matrix_list.txt','r'):
                coverage_matrix_list.append(str(line).strip())
            for line in open('Assembly_mo_list.txt','r'):
                assembly_mo_list.append(str(line).strip())
            for line in open('Assembly_MoDict.txt','r'):
                item=str(line).strip().split('\t')[0].strip()
                connections_list.append('condense_connections_'+item+'.txt')
                
            if sensitive == 'more-sensitive':
                refinement_paramter = 'deep'
            # elif refinement_paramter != 'deep':
            #     refinement_paramter = 'quick'
            Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, continue_mode, min_cpn, max_ctn, assembly_mo_list, connections_list, coverage_matrix_list, refinement_paramter, pwd)
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'5th contig retrieve done!')
            f_cp_m.close()
            os.system('cp -r BestBinset_outlier_refined_filtrated_retrieved BestBinset_outlier_refined_filtrated_retrieved_backup')

        if last_step < 6:
            ### Second de-replication: S4 function
            if len(assembly_list) != 1:
                print('Secondary de-repplication start')
                #drep_list='BestBinset_outlier_refined_filtrated_retrieved'

                coverage_matrix_list, bestbinset_list = [], []
                for line in open('Bestbinset_list.txt','r'):
                    bestbinset_list.append(str(line).strip())
            
                for line in open('Coverage_matrix_list.txt','r'):
                    coverage_matrix_list.append(str(line).strip())

                datasets_fq={}
                for item in datasets.keys():
                    datasets_fq[item]=pe_read_names(datasets[item])

                # multiple_assembly_comparitor_main(drep_list, bestbinset_list, coverage_matrix_list, datasets_fq, 'second_drep', num_threads)
                final_binset_comparitor('BestBinset_outlier_refined_filtrated_retrieved', coverage_matrix_list, datasets_fq, num_threads, pwd, 'second_drep')
            f_cp_m=open('Basalt_checkpoint.txt', 'a')
            f_cp_m.write('\n'+'6th secondary de-repplication done!')
            f_cp_m.close()

        if last_step < 7:
            if len(lr_list) == 0:
                print('Starting contig retrival within group')
                # Contig_recruiter_main('BestBinset_outlier_refined', num_threads, parameter, 35, 20)
                best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                outlier_remover_folder='BestBinset_outlier_refined'
                cpn_cutoff=35 ### The minimun completeness to keep
                ctn_cutoff=5 ### The maximal contaminaition to keep

                coverage_matrix_list, connections_list, assembly_mo_list = [], [], []
                for line in open('Coverage_matrix_list.txt','r'):
                    coverage_matrix_list.append(str(line).strip())
                for line in open('Assembly_mo_list.txt','r'):
                    assembly_mo_list.append(str(line).strip())
                for line in open('Assembly_MoDict.txt','r'):
                    item=str(line).strip().split('\t')[0].strip()
                    connections_list.append('condense_connections_'+item+'.txt')

                # Contig_retrieve_within_group_main(best_binset_after_contig_retrieve, outlier_remover_folder, num_threads, continue_mode, cpn_cutoff, ctn_cutoff, assembly_mo_list, connections_list, coverage_matrix_list)

                f_cp_m=open('Basalt_checkpoint.txt', 'a')
                f_cp_m.write('\n'+'7th Skip contig retrieve within group!'+'\t'+best_binset_after_contig_retrieve)
                # f_cp_m.write('\n'+'7th contig retrieve within group done!'+'\t'+best_binset_after_contig_retrieve+'_retrieved')
                f_cp_m.close()
            else:
                ### Long read present. Start polishing
                print('Starting 1st run polishing')
                if functional_module == 'all':
                    best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                    polishing_main(best_binset_after_contig_retrieve, datasets, assembly_list, lr_list, 1, 'bw2', 3, ram, num_threads)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'7th Polishing done!'+'\t'+best_binset_after_contig_retrieve+'_MAGs_polished')
                    f_cp_m.close()
                elif functional_module == 'refinement':
                    best_binset_after_contig_retrieve='BestBinset_outlier_refined_filtrated_retrieved'
                    polishing_main(best_binset_after_contig_retrieve, datasets, assembly_list, lr_list, 1, 'bwa', 3, ram, num_threads)
                    f_cp_m=open('Basalt_checkpoint.txt', 'a')
                    f_cp_m.write('\n'+'7th Polishing done!'+'\t'+best_binset_after_contig_retrieve+'_MAGs_polished')
                    f_cp_m.close()
    print('Refinemnt done!')

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
    BASALT_main_c_refinement(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitive, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
    # BASALT_main(assembly_list, datasets, num_threads, lr_list, ram, continue_mode, functional_module, autobining_parameters, refinement_paramter, hybri_reassembly, max_ctn, min_cpn, pwd)
