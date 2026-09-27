#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S1p: Merge closely related bins within the same assembly group.

This module analyses paired-end connection files, depth matrices, and
CheckM statistics to identify candidate bin pairs that likely belong
to the same genome, and records them for downstream merging.
"""

from Bio import SeqIO
import sys, os, threading, copy
from multiprocessing import Pool
from collections import Counter
from time import ctime,sleep
from basalt_runtime import depth_sample_count


def seq_recorder(bin_folder, pwd):
    """
    Index contig IDs for each bin in the given folder.

    Returns
    -------
    seqs_record : dict
        Mapping contig_id -> pipe-separated list of bin filenames.
    file_seqs_record : dict
        Mapping bin filename -> list of contig IDs it contains.
    """
    seqs_record, file_seqs_record = {}, {}
    for root,dirs,files in os.walk(pwd+'/'+bin_folder):
        for file in files:
            # print('Reading', file
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                file_seqs_record[file]=[]
                for record in SeqIO.parse(pwd+'/'+bin_folder+'/'+str(file),'fasta'):
                    file_seqs_record[file].append(record.id)
                    try:
                        seqs_record[str(record.id)]+='|'+file
                    except:
                        seqs_record[str(record.id)]=file
    return seqs_record, file_seqs_record

def bin_group(PE_connection_file, seqs_record, pwd):
    """
    Build bin–bin connection scores based on paired-end connections.

    Parameters
    ----------
    PE_connection_file : str
        Path to the PE connection summary file.
    seqs_record : dict
        Mapping from contig IDs to bin filenames as produced by
        ``seq_recorder``.
    pwd : str
        Working directory.

    Returns
    -------
    dict
        Mapping bin1 -> bin2 for candidate bin pairs that pass
        connection score thresholds.
    """
    genome_connection, bin_group1, m={}, {}, 0
    for line in open(pwd+'/'+str(PE_connection_file), 'r'):
        m+=1
        if m >= 2:
            node1=str(line).strip().split('\t')[0]
            node2=str(line).strip().split('\t')[2]
            num_connections=float(str(line).strip().split('\t')[3])
            try:
                node1_bin_list=seqs_record[str(node1)].split('|')
                node2_bin_list=seqs_record[str(node2)].split('|')
                for bin1 in node1_bin_list:
                    for bin2 in node2_bin_list:
                        if bin1 != bin2:
                            x=[bin1,bin2]
                            x.sort()
                            y=x[0]+'\t'+x[1]

                            try:
                                genome_connection[y]+=num_connections
                            except:
                                genome_connection[y]=num_connections
                            
                            try:
                                bin_group1[x[0]][y]=genome_connection[y]
                            except:
                                bin_group1[x[0]]={}
                                bin_group1[x[0]][y]=genome_connection[y]
            except:
                xyzzz=0

    f=open('Bin_connections.txt','a')
    for i,j in bin_group1.items():
        for y in j.keys():
            f.write(str(i)+'\t'+str(y)+'\t'+str(genome_connection[y])+'\n')
    f.close()

    #### Ranking the score
    bin_group2={}
    for i in bin_group1.keys():
        bin_group2[i]={}
        xyz=0
        maxxyz=xyz
        for y in bin_group1[i].keys():
            if bin_group1[i][y] >= 50 :
                if bin_group1[i][y] > maxxyz:
                    bin_group2[i][y]=bin_group1[i][y]
                    maxxyz=bin_group1[i][y]
    
    bin_group3=copy.deepcopy(bin_group2)
    for i in bin_group3.keys():
        for y in bin_group3[i].keys():
            if bin_group3[i][y] == 0:
                del bin_group2[i][y]

    pair_bins={}
    if len(bin_group2) != 0:
        f=open('Bin_connections_selected.txt','a')
        for i,j in bin_group2.items():
            for y in j.keys():
                bin1=str(y).split('\t')[0]
                bin2=str(y).split('\t')[1]
                pair_bins[bin1]=bin2
                f.write(str(i)+'\t'+str(y)+'\t'+str(genome_connection[y])+'\n')
        f.close()
    return pair_bins

def parse_checkm(bin_folder_checkm):
    """
    Parse CheckM bin statistics from ``bin_stats_ext.tsv``.

    Returns a dictionary keyed by bin filename ('.fa' / '.fasta').
    """
    checkm={}
    for line in open(str(bin_folder_checkm)+'/storage/bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        # genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].strip()
        # GC=round(float(str(line).strip().split('\'GC\':')[1].split(',')[0].strip())*100, 1)
        # Mean_scaffold_length=str(line).strip().split('Mean scaffold length\': ')[1].split(',')[0].split('}')[0].strip()
        checkm[str(binID)+'.fa']['completeness']=float(completeness)
        checkm[str(binID)+'.fa']['contamination']=float(contamination)
        checkm[str(binID)+'.fa']['lineage']=taxon
        checkm[str(binID)+'.fasta']['completeness']=float(completeness)
        checkm[str(binID)+'.fasta']['contamination']=float(contamination)
        checkm[str(binID)+'.fasta']['lineage']=taxon
    return checkm

def depth_eval(file_seqs_record, depth_file, pair_bins, bin_folder):
    """
    Evaluate candidate bin pairs using per-sample depth profiles.

    For each bin pair, this function calculates total and average
    coverage across samples and records pairs with consistent depth
    patterns and compatible CheckM taxonomy/completeness.
    """
    n, depth_matrix, num = 0, {}, 0
    for line in open(depth_file,'r'):
        n+=1
        if n == 1:
            num=depth_sample_count(line)
        else:
            contig=str(line).strip().split('\t')[0]
            depth_matrix[contig]=[]
            for i in range(0, num):
                cov=float(str(line).strip().split('\t')[3+i*2])
                depth_matrix[contig].append(cov)
    
    f=open('Pair_bins_'+bin_folder+'.txt','w')
    f.write('Bin1'+'\t'+'Bin2'+'\t'+'Bin1 contig NO.'+'\t'+'Bin2 contig NO.')
    for i in range(1, num+1):
        f.write('\t'+'Bin1 total cov-'+str(i)+'\t'+'Bin1 avg. cov-'+str(i)+'\t'+'Bin2 total cov-'+str(i)+'\t'+'Bin2 avg. cov-'+str(i)+'\t'+'Var(%)')
    f.write('\n')

    bin1_depth, bin2_depth = {}, {}
    for bin1 in pair_bins.keys():
        bin2=pair_bins[bin1]
        f.write(str(bin1)+'\t'+str(bin2))

        x=0
        for contig in file_seqs_record[bin1]:
            x+=1
            for i in range(1, num+1):
                try:
                    bin1_depth[i]+=depth_matrix[contig][i-1]
                except:
                    bin1_depth[i]=depth_matrix[contig][i-1]

        y=0
        for contig in file_seqs_record[bin2]:
            y+=1
            for i in range(1, num+1):
                try:
                    bin2_depth[i]+=depth_matrix[contig][i-1]
                except:
                    bin2_depth[i]=depth_matrix[contig][i-1]
        f.write('\t'+str(x)+'\t'+str(y))

        judge=0
        for i in bin1_depth.keys():
            bin1_depth_avg=bin1_depth[i]/x
            bin2_depth_avg=bin2_depth[i]/y
            f.write('\t'+str(bin1_depth[i])+'\t'+str(bin1_depth_avg)+'\t'+str(bin2_depth[i])+'\t'+str(bin2_depth_avg))
            delta=abs(bin1_depth_avg-bin2_depth_avg)
            mean=(bin1_depth_avg+bin2_depth_avg)/2
            perc=100*delta/mean
            f.write('\t'+str(perc))
            if perc <= 20:
                judge+=1
        f.write('\n')

        pair_bin_passed={}
        if judge == num:
            f2=open('Potential_paired_bins_selected.txt','a')
            name_list=bin_folder.split('_')
            name_list.remove(name_list[-1])
            bin_folder_checkm='_'.join(name_list)+'_checkm'
            checkm=parse_checkm(bin_folder_checkm)
            bin1_tax=bin_folder_checkm[bin1]['lineage']
            bin2_tax=bin_folder_checkm[bin2]['lineage']
            bin1_cmp=bin_folder_checkm[bin1]['completeness']
            bin1_ctn=bin_folder_checkm[bin1]['contamination']
            bin2_cmp=bin_folder_checkm[bin2]['completeness']
            bin2_ctn=bin_folder_checkm[bin1]['contamination']
            f2.write(bin1+'\t'+bin1_tax+'\t'+bin1_cmp+'\t'+bin1_ctn+'\t'+bin2+'\t'+bin2_tax+'\t'+bin2_cmp+'\t'+bin2_ctn+'\n')
            cmp_sum=bin1_cmp+bin2_cmp
            if cmp_sum <= 105 and bin1_tax == bin2_tax:
                pair_bin_passed[bin1]=bin2
        
    f.close()
    f2.close()

def merge_bin_within_the_same_assembly(assembly_binning_group, depth_files,
                                       PE_connections_files, assembly_names,
                                       num_threads):
    """
    Top-level driver to identify mergeable bin pairs within one assembly.
    """
    print('-------------------------------')
    print(str(assembly_binning_group))
    print(str(depth_files))
    print(str(PE_connections_files))
    print(str(assembly_names))
    print('-------------------------------')

    f=open('Bin_connections.txt','w')
    f.write('Group'+'\t'+'Bin1'+'\t'+'Bin2'+'\t'+'Connections'+'\n')
    f.close()

    f1=open('Bin_connections_selected.txt','w')
    f1.write('Group'+'\t'+'Bin1'+'\t'+'Bin2'+'\t'+'Connections'+'\n')
    f1.close()

    f=open('Potential_paired_bins_selected.txt','w')
    f.write('Bin1'+'\t'+'Lineage'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'Bin2'+'\t'+'Lineage'+'\t'+'Completeness'+'\t'+'Contamination'+'\n')
    f.close()
    
    pwd=os.getcwd()
    for item in assembly_binning_group.keys():
        for bin_folder in assembly_binning_group[item]:
            seqs_record={}
            A=seq_recorder(bin_folder, pwd)
            seqs_record.update(A[0])
            file_seqs_record=A[1]
        
            pair_bins={}
            PE_connections_file=PE_connections_files[item]
            pair_bins.update(bin_group(PE_connections_file, seqs_record, pwd))

            if len(pair_bins) != 0:
                depth_file=depth_files[item]
                depth_eval(file_seqs_record, depth_file, pair_bins, bin_folder)

if __name__ == '__main__': 
    assembly_binning_group={'1':['1_8_medium_cat_SPAdes_scaffolds.fasta_0.3_maxbin2_genomes', '1_8_medium_cat_SPAdes_scaffolds.fasta_300_metabat_genomes','1_8_medium_cat_SPAdes_scaffolds.fasta_100_concoct_genomes'],'2':['2_9_medium_S001_SPAdes_scaffolds.fasta_0.3_maxbin2_genomes', '2_9_medium_S001_SPAdes_scaffolds.fasta_200_metabat_genomes','2_9_medium_S001_SPAdes_scaffolds.fasta_146_concoct_genomes'],'3':['3_10_medium_S002_SPAdes_scaffolds.fasta_0.3_maxbin2_genomes', '3_10_medium_S002_SPAdes_scaffolds.fasta_200_metabat_genomes','3_10_medium_S002_SPAdes_scaffolds.fasta_152_concoct_genomes']
    }

    depth_files={'1':'1_assembly.depth.txt','2':'2_assembly.depth.txt','3':'3_assembly.depth.txt'}
    PE_connections_files={'1':'condense_connections_8_medium_cat_SPAdes_scaffolds.fasta.txt','2':'condense_connections_9_medium_S001_SPAdes_scaffolds.fasta.txt','3':'condense_connections_10_medium_S002_SPAdes_scaffolds.fasta.txt'}
    assembly_names={'1':'1_8_medium_cat_SPAdes_scaffolds.fasta','2':'2_9_medium_S001_SPAdes_scaffolds.fasta','3':'3_10_medium_S002_SPAdes_scaffolds.fasta'}
    num_threads=2

    merge_bin_within_the_same_assembly(assembly_binning_group, depth_files, PE_connections_files, assembly_names, num_threads)
