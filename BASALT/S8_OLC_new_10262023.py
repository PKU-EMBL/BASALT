#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S8 (CheckM2 branch): OLC-based reassembly of selected bins.

This module performs overlap-layout-consensus (OLC) reassembly of
bins using selected contigs to improve continuity and assembly
quality.
"""

from Bio import SeqIO
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna
    generic_dna_t=1
except ImportError:
    generic_dna_t=0
import os, copy, glob
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from multiprocessing import Pool
from basalt_runtime import (
    clean_stale_olc_intermediates,
    evaluate_bins_checkm2,
    expand_alignment_group,
    fasta_has_sequence,
    positive_thread_count,
    worker_thread_budget,
)


def elongate_contig_selector(eliminated_bin, threshold, pwd,
                             eliminated_bin_containing_folder):
    """
    Select contigs within an eliminated bin that are suitable for OLC elongation.

    This function loads an eliminated bin, computes per-contig coverage,
    applies PCA-based outlier detection, and writes new bin candidates
    where contigs marked as outliers are removed.

    Parameters
    ----------
    eliminated_bin : str
        File name of the bin FASTA to be evaluated (within the
        ``eliminated_bin_containing_folder``).
    threshold : float
        Outlier threshold used in the PCA / IQR based filter.
    pwd : str
        Working directory path.
    eliminated_bin_containing_folder : str
        Folder in which the eliminated bin FASTA and coverage files reside.

    Returns
    -------
    list of str
        List of file names for newly generated candidate bins after
        removing outlier contigs. The list may be empty if no candidates
        are produced.
    """
    # def elongate_contig_selector(eliminated_bin, threshold, pwd, eliminated_bin_containing_folder, all_kmer):
    try:
        vs_contig_seq, num_contig, i = {}, 0, 1
        while i <= 3 and len(vs_contig_seq) == 0:
            i+=1
            os.system('cp '+pwd+'/'+eliminated_bin_containing_folder+'/'+str(eliminated_bin)+' '+pwd)
            for record in SeqIO.parse(eliminated_bin,'fasta'):
                num_contig+=1
                vs_contig_seq[record.id]=record.seq
            # os.system('calc.kmerfreq.pl -i '+str(eliminated_bin)+' -o '+str(eliminated_bin)+'.kmer.txt')
            if i == 4 and len(vs_contig_seq) == 0:
                fbin_record_error=open('Bin_record_error.txt','a')
                fbin_record_error.write('Recorded '+str(eliminated_bin)+' error'+'\n')
                fbin_record_error.close()

        depth_file=str(eliminated_bin).split('_genomes.')[0].split('_assembly.fa_')[0].split('_')[0]+'_assembly.depth.txt'

        n, contig_cov = 0, {}
        for line in open(depth_file,'r'):
            n+=1
            if n == 1:
                ls=str(line).strip().split('\t')
                num=(len(ls)-4)/3
            else:
                ids=str(line).strip().split('\t')[0]
                if ids in vs_contig_seq.keys():
                    contig_cov[ids]={}
                    for i in range(1,int(num)+1):
                        contig_cov[ids][i]=float(str(line).strip().split('\t')[3*i+1])

        pwd=os.getcwd()
        print('Transfroming coverage matrix')
        os.system('mkdir '+eliminated_bin+'_outlier')
        bin_outlier, coverage_data, contigs_ids, coverage_list={}, {}, [], []
        for contig in contig_cov.keys():
            contigs_ids.append(contig)
            num_coverage=len(contig_cov[contig])
            for i in range(1, num_coverage+1):
                if i not in coverage_data.keys():
                    coverage_data[i]=[]
                coverage_data[i].append(contig_cov[contig][i])
                coverage_list.append(contig_cov[contig][i])

        coverage_array=np.array(coverage_list).reshape((num_contig,num_coverage))

        A=PCA_slector(coverage_array, num_contig)
        newData=A[0]
        explained_variance_ratio=A[1]
        # bin_outlier=outliner_remover(eliminated_bin, contigs_ids, threshold, newData, explained_variance_ratio, bin_outlier)
        seperate_outlier=outliner_remover(eliminated_bin, contigs_ids, threshold, newData, explained_variance_ratio, pwd)
        bin_outlier.update(seperate_outlier)

        # print('Calculating TNFs of', eliminated_bin)
        # bin_TNFs_outlier={}
        # Bins_TNFs, bin_contig_list = [], []
        # for contig in vs_contig_seq.keys():
        #     lis=str(all_kmer[contig]).split('\t')
        #     bin_contig_list.append(contig)
        #     for i in range(1, len(lis)):
        #         Bins_TNFs.append(lis[i])

        # TNF_array=np.array(Bins_TNFs).reshape((num_contig, 256))
        # A=PCA_slector(TNF_array, num_contig)
        # newData=A[0]
        # explained_variance_ratio=A[1]
        # seperate_outlier=outliner_remover(eliminated_bin, bin_contig_list, threshold, newData, explained_variance_ratio, pwd)
        # bin_TNFs_outlier.update(seperate_outlier)

        f=open('Record_coverage_outlier_'+eliminated_bin+'.txt', 'w')
        # f1=open('Record_TNFs_outlier_'+eliminated_bin+'.txt', 'w')
        selected_bins=[]
        for item in bin_outlier.keys():
            f.write(str(item)+'\t'+str(bin_outlier[item])+'\n')
            # f1.write(str(item)+'\t'+str(bin_TNFs_outlier[item])+'\n')
            f2=open(str(item)+'_'+eliminated_bin,'w')
            selected_bins.append(str(item)+'_'+eliminated_bin)
            for item2 in vs_contig_seq.keys():
                # if item2 not in bin_TNFs_outlier.keys():
                # if item2 not in bin_outlier[item].keys() and item2 not in bin_TNFs_outlier.keys():
                if item2 not in bin_outlier[item].keys():
                    f2.write('>'+str(item2)+'\n'+str(vs_contig_seq[item2])+'\n')
            f2.close()
        f.close()  
        # f1.close()
    except:
        selected_bins=[]
    os.system('mv Record_coverage_outlier_'+eliminated_bin+'.txt '+pwd+'/'+eliminated_bin+'_outlier')
    # os.system('mv Record_coverage_outlier_'+eliminated_bin+'.txt Record_TNFs_outlier_'+eliminated_bin+'.txt '+pwd+'/'+eliminated_bin+'_outlier')
    return selected_bins

def outliner_remover(bin_id, contigs_ids, threshold, item_data, explained_variance_ratio, pwd):
    """
    Identify outlier contigs using an IQR-based rule on projected features.

    Parameters
    ----------
    bin_id : str
        Identifier of the bin being evaluated.
    contigs_ids : list of str
        List of contig identifiers corresponding to ``item_data`` entries.
    threshold : list of float
        List of multiplicative factors applied to IQR when computing
        upper/lower bounds (one value per run).
    item_data : list of float
        One-dimensional feature vector (e.g. first PCA component values)
        for each contig in ``contigs_ids``.
    explained_variance_ratio : numpy.ndarray
        Explained variance ratios of the PCA transformation used to
        generate ``item_data``, written into the summary files.
    pwd : str
        Working directory path where per-bin outlier folders are created.

    Returns
    -------
    dict
        Mapping ``threshold_value -> dict(contig_id -> feature_value)``
        describing which contigs were flagged as outliers for each
        threshold.
    """
    print('Finding outliner from', bin_id)
    four = pd.Series(item_data).describe()
    bin_outlier={}
    for item in threshold:
        bin_outlier[item]={}
        f1=open('Outlier_in_threshold'+str(item)+'_'+str(bin_id)+'.txt', 'w')
        f2=open('Summary_threshold'+str(item)+'_'+str(bin_id)+'.txt', 'w')
        #print(four)
        #print('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))
        Q1 = four['25%']
        Q3 = four['75%']
        IQR = Q3 - Q1
        upper1 = Q3 + float(item) * IQR
        lower1 = Q1 - float(item) * IQR
        #print(upper1, lower1)
        n1, outliner_record=0, {}
        for i in range(0, len(item_data)):
            if item_data[i] > float(upper1) or item_data[i] < float(lower1):
                n1+=1
                f1.write(str(contigs_ids[i])+'\t'+str(item_data[i])+'\n')
                bin_outlier[item][str(contigs_ids[i])]=str(item_data[i])
                # if str(item) not in bin_outlier.keys():
                #     bin_outlier[str(item)]=[]
                #     bin_outlier[str(item)].append(str(contigs_ids[i]))
                # else:
                #     bin_outlier[str(item)].append(str(contigs_ids[i]))
        f1.close()
        f2.write(str(four)+'\n'+str('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))+'\n'+'Upper:'+str(upper1)+'\t'+'Lower:'+str(lower1)+'\n'+str(n1)+' outliers in '+str(bin_id)+' under the threshold of '+str(item)+'\n'+'Explained variance ratio:'+str(explained_variance_ratio)+'\n')
        f2.close()
        os.system('mv Outlier_in_threshold'+str(item)+'_'+str(bin_id)+'.txt Summary_threshold'+str(item)+'_'+str(bin_id)+'.txt '+pwd+'/'+str(bin_id)+'_outlier')
        print(n1, 'outliers in', str(bin_id), 'with threshold of', item)
        print('-------------------------')
    return bin_outlier

def PCA_slector(data_array, num_contig):
    """
    Project feature matrix to a single principal component.

    Parameters
    ----------
    data_array : numpy.ndarray
        Two-dimensional array of shape ``(num_contig, num_features)``
        containing per-contig features (coverage, TNF, etc.).
    num_contig : int
        Number of contigs represented in ``data_array``.

    Returns
    -------
    list of float
        Single principal-component score per contig.
    numpy.ndarray
        Array of explained variance ratios for the PCA model.
    """
    pca = PCA(n_components=1)
    pca.fit(data_array)
    explained_variance_ratio=pca.explained_variance_ratio_
    print(explained_variance_ratio)
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

def record_seq(target_bin, eliminated_bin):
    """
    Concatenate sequences from a target bin and an eliminated bin.

    Parameters
    ----------
    target_bin : str
        File name of the target bin FASTA.
    eliminated_bin : str
        File name of the eliminated bin FASTA whose contigs may be
        recruited or merged.

    Returns
    -------
    dict
        Mapping ``contig_id -> sequence_str`` for the target bin.
    dict
        Mapping ``contig_id -> length`` for the target bin.
    dict
        Mapping ``contig_id -> sequence`` for the eliminated bin.
    dict
        Mapping ``contig_id -> length`` for the eliminated bin.
    str
        Name of the merged FASTA file written to the working directory.
    dict
        Mapping ``contig_id -> sequence`` for all contigs from both bins.
    dict
        Dictionary used to record which merged bin files have been
        generated; keys are merged file names.
    """
    merged_bin_recorded={}
    f=open('Merged_'+target_bin+'_'+eliminated_bin,'w')
    merged_bin_recorded['Merged_'+target_bin+'_'+eliminated_bin]=0
    target_contig_len, target_contig_seq, total_seq={}, {}, {}
    for record in SeqIO.parse(target_bin,'fasta'):
        target_contig_len[record.id]=len(record.seq)
        target_contig_seq[record.id]=str(record.seq)
        total_seq[record.id]=record.seq
        f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')

    vs_contig_len, vs_contig_seq={}, {}
    for record in SeqIO.parse(eliminated_bin,'fasta'):
        vs_contig_len[record.id]=len(record.seq)
        vs_contig_seq[record.id]=record.seq
        total_seq[record.id]=record.seq
        f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    f.close()

    # print(len(target_contig_seq['PB_1-15800'])
    return target_contig_seq, target_contig_len, vs_contig_seq, vs_contig_len, 'Merged_'+target_bin+'_'+eliminated_bin, total_seq, merged_bin_recorded

def blast_1(target_bin, eliminated_bin, target_contig_seq, target_contig_len,
            vs_contig_seq, vs_contig_len, aligned_len_cutoff,
            similarity_cutoff, num_threads, blast_name, folder_name):
    """
    Run BLAST between target and eliminated bins and group similar contigs.

    Parameters
    ----------
    target_bin : str
        Path to the FASTA file of the target bin being elongated.
    eliminated_bin : str
        Path to the FASTA file of the eliminated bin to search against.
    target_contig_seq : dict
        Mapping ``contig_id -> sequence_str`` for target contigs.
    target_contig_len : dict
        Mapping ``contig_id -> length`` for target contigs.
    vs_contig_seq : dict
        Mapping ``contig_id -> sequence`` for contigs from the eliminated bin.
    vs_contig_len : dict
        Mapping ``contig_id -> length`` for contigs from the eliminated bin.
    aligned_len_cutoff : int
        Minimum aligned length to consider a BLAST hit valid.
    similarity_cutoff : float
        Minimum percentage identity for a BLAST hit to be kept.
    num_threads : int
        Number of BLAST threads to use (the implementation currently
        sets it to 1 internally).
    blast_name : str
        Base name for BLAST output files.
    folder_name : str
        Folder where filtered BLAST outputs and helper files will be moved.

    Returns
    -------
    dict
        Mapping ``group_id -> list(stripped_alignment_lines)`` describing
        groups of contigs that are similar enough to be considered in
        the same OLC elongation group.
    """
    pwd = os.getcwd()
    os.system('makeblastdb -in '+eliminated_bin+' -dbtype nucl -hash_index -parse_seqids -logfile '+eliminated_bin+'_db.txt')
    # os.system('blastn -query '+target_bin+' -db '+eliminated_bin+' -evalue 1e-20 -num_threads '+str(num_threads)+' -outfmt 6 -out '+str(blast_name))
    os.system('blastn -query '+target_bin+' -db '+eliminated_bin+' -evalue 1e-20 -num_threads '+str(positive_thread_count(num_threads))+' -outfmt 6 -out '+str(blast_name))
    os.system('rm *_db.txt')

    f=open(str(blast_name)+'_filtration_1.txt','w')
    f2=open(str(blast_name)+'_filtration_start_end_1.txt','w')
    f3=open('BLAST_output_error.txt','a')
    query_id, subject_id = {}, {}
    for line in open(blast_name,'r'):
        query=str(line).strip().split('\t')[0]
        subject=str(line).strip().split('\t')[1]
        try:
            similarity=float(str(line).strip().split('\t')[2])
            aligned_length=int(str(line).strip().split('\t')[3])
        except:
            similarity=1
            aligned_length=1

        try:
            query_start=int(str(line).strip().split('\t')[6])
            query_end=int(str(line).strip().split('\t')[7])
            subject_start=int(str(line).strip().split('\t')[8])
            subject_end=int(str(line).strip().split('\t')[9])
        except:
            query_start=1
            query_end=0
            subject_start=1
            subject_end=0
            f3.write(str(line))

        if similarity >= similarity_cutoff and aligned_length >= aligned_len_cutoff:
            f.write(line)
            try:
                if query_start == 1 or query_start == target_contig_len[query] or query_end == 1 or query_end == target_contig_len[query]:
                    if subject_start == 1 or subject_end == vs_contig_len[subject] or subject_start == vs_contig_len[subject] or subject_end == 1:
                        f2.write(line)
                        if query not in query_id.keys():
                            query_id[query]=1
                        else:
                            query_id[query]+=1
            
                        if subject not in subject_id.keys():
                            subject_id[subject]=1
                        else:
                            subject_id[subject]+=1
            except:
                ferror=open('OLC_merged_error_blast_results.txt','a')
                ferror.write(str(target_bin)+' extending with '+str(eliminated_bin)+' error blast result: '+'\n'+str(line).strip()+'\n')
                ferror.close()
    f2.close()
    f.close()
    f3.close()
    
    redance_id, r_query, alignment_dict={}, {}, {}
    f=open(str(blast_name)+'_filtration_merge_1.txt','w')
    try:
        for line in open(str(blast_name)+'_filtration_start_end_1.txt','r'):
            query=str(line).strip().split('\t')[0]
            subject=str(line).strip().split('\t')[1]
            if query_id[query] > 1 or subject_id[subject] > 1:
                f.write(str(line))
                r_query[query]=1
                alignment_dict[query+' '+subject]=str(line).strip()
    except:
        xyzzzz=0
    f.close()

    f=open(str(blast_name)+'_nr_seq.txt','w')
    for contigs in target_contig_seq.keys():
        if contigs not in r_query.keys():
            f.write('>'+str(contigs)+'\n'+str(target_contig_seq[contigs])+'\n')
    f.close()

    print(str(blast_name)+' splitting blast output')
    blast_group={}
    for item in r_query.keys():
        num2, num1=0, 1
        blast_group[item]={}
        blast_group[item][item]=1
        blast_group[item], grouping_converged=expand_alignment_group(blast_group[item], alignment_dict)
        if not grouping_converged:
            print('Warning: BLAST grouping stopped at the iteration limit for '+str(item))

    nr_blast_group={}
    for item in blast_group.keys():
        test_list=[]
        for item2 in blast_group[item].keys():
            test_list.append(item2)
        test_list.sort()
        test_dict={}
        for item3 in test_list:
            test_dict[item3]=1
        test=str(test_dict)
        if test not in nr_blast_group.keys():
            nr_blast_group[test]=1
        else:
            continue
        
    # print(len(nr_blast_group)
    f=open(str(blast_name)+'_possible_similar_contig_group.txt','w')
    try:
        os.mkdir(target_bin+'_split_blast_output')
    except:
        print(target_bin+'_split_blast_output exists')
    n, blast_group, mv_file=1, {}, {}
    for item in nr_blast_group.keys():
        blast_group[n]=[]
        filtration_query={}
        f2=open(str(blast_name)+'_group_'+str(n)+'_split_blast_output.txt','w')
        mv_file[str(blast_name)+'_group_'+str(n)+'_split_blast_output.txt']=0
        for line in open(str(blast_name)+'_filtration_merge_1.txt','r'):
            query=str(line).strip().split('\t')[0]
            subject=str(line).strip().split('\t')[1]
            if '\''+query+'\'' in item or '\''+subject+'\'' in item:
                f2.write(line)
                blast_group[n].append(line.strip())
        f2.close()

        for line in open(str(blast_name)+'_group_'+str(n)+'_split_blast_output.txt','r'):
            query=str(line).strip().split('\t')[0]
            if query in target_contig_seq.keys():
                filtration_query[query]=''
        
        if len(filtration_query) < 2: ### Query from blast output shall have at least two contigs from target_contig_seq
            os.system('rm '+str(blast_name)+'_group_'+str(n)+'_split_blast_output.txt')
            del blast_group[n]
            del mv_file[str(blast_name)+'_group_'+str(n)+'_split_blast_output.txt']
        else:
            f.write(str(n)+'\t'+str(item)+'\n')
            n+=1
    f.close()

    for item in mv_file.keys():
        os.system('mv '+str(item)+' '+target_bin+'_split_blast_output')
    
    os.system('mv '+str(blast_name)+'_nr_seq.txt '+str(blast_name)+'_filtration_1.txt '+str(blast_name)+'_filtration_start_end_1.txt '+str(blast_name)+'_possible_similar_contig_group.txt '+str(blast_name)+'_filtration_merge_1.txt '+folder_name)
    return blast_group

def elongation_sub_contig(merged_seq, query_seq, iteration_num,
                          aligned_len_cutoff, similarity_cutoff,
                          num_threads):
    """
    Refine an already merged contig by aligning it back to the query contig.

    Parameters
    ----------
    merged_seq : str
        Path to the FASTA file containing the merged sequence from the
        previous elongation iteration.
    query_seq : str
        Path to the FASTA file of the original query contig.
    iteration_num : int
        Current iteration index of the elongation process.
    aligned_len_cutoff : int
        Minimum alignment length for BLAST hits.
    similarity_cutoff : float
        Minimum percentage identity for BLAST hits.
    num_threads : int
        Number of threads for BLAST.

    Returns
    -------
    int
        Status flag indicating whether a further elongation iteration
        should be attempted (1) or the process should stop (0).
    """
    pwd = os.getcwd()

    query_len, merged_seq_len=0, 0
    for record in SeqIO.parse(query_seq,'fasta'):
        query_len=len(record.seq)

    for record in SeqIO.parse(merged_seq,'fasta'):
        merged_seq_len=len(record.seq)

    ### positive direction blast
    os.system('makeblastdb -in '+str(query_seq)+' -dbtype nucl -hash_index -parse_seqids -logfile '+str(query_seq)+'_db.txt')
    os.system('blastn -query '+str(merged_seq)+' -db '+str(query_seq)+' -evalue 1e-20 -num_threads '+str(num_threads)+' -outfmt 6 -out '+str(merged_seq)+'_blast_'+str(iteration_num)+'_'+str(iteration_num)+'.txt')

    ### BLAST filtration
    positive_blast, negative_blast, n1, n2, aligned_len1, aligned_len2 = 0, {}, 0, 0, 0, 0
    for line in open(str(merged_seq)+'_blast_'+str(iteration_num)+'_'+str(iteration_num)+'.txt','r'):
        n1+=1
        if n1 == 1:
            query=str(line).strip().split('\t')[0]
            subject=str(line).strip().split('\t')[1]
            similarity=float(str(line).strip().split('\t')[2])
            aligned_length=int(str(line).strip().split('\t')[3])
            query_start=int(str(line).strip().split('\t')[6])
            query_end=int(str(line).strip().split('\t')[7])
            subject_start=int(str(line).strip().split('\t')[8])
            subject_end=int(str(line).strip().split('\t')[9])
            if similarity >= similarity_cutoff and aligned_length >= aligned_len_cutoff:
                if query_start == 1 or query_start == merged_seq_len or query_end == 1 or query_end == merged_seq_len:
                    if subject_start == 1 or subject_end == query_len or subject_start == query_len or subject_end == 1:
                        positive_blast=str(line).strip()
                    
    os.system('rm '+str(query_seq)+' '+str(merged_seq)+'_blast_'+str(iteration_num)+'_'+str(iteration_num)+'.txt')
    # os.system('rm *.nhd *.nin *.nsq *.nhr *.nog *.nsd *.nsi *.nhi *.ndb *.nos *.not *.ntf *.nto *.perf '+str(query_seq)+' blast_'+str(iteration_num)+'_'+str(iteration_num)+'.txt')
    return positive_blast

def blast_2(target_bin,target_contig_seq, merged_bin, total_seq, threshold_item, iteration_num, num_threads, folder_name):
    os.system('makeblastdb -in '+str(merged_bin)+' -dbtype nucl -hash_index -parse_seqids -logfile '+str(merged_bin)+'_db.txt')
    os.system('blastn -query '+str(target_bin)+' -db '+str(merged_bin)+' -evalue 1e-20 -num_threads '+str(positive_thread_count(num_threads))+' -outfmt 6 -out blast_'+str(target_bin)+'_self_merged_'+str(threshold_item)+'.txt')
    # os.system('blastn -query '+str(target_bin)+' -db '+str(merged_bin)+' -evalue 1e-20 -num_threads '+str(num_threads)+' -outfmt 6 -out blast_'+str(target_bin)+'_self_merged_'+str(threshold_item)+'.txt')
    os.system('rm *.perf')

    merged_bin_seq={}
    for record in SeqIO.parse(merged_bin,'fasta'):
        merged_bin_seq[record.id]=record.seq

    being_merged_contigs, being_merged_contigs2, being_merged_contigs3={}, {}, {}
    f=open('Filtrated_blast_'+str(target_bin)+'_self_merged_'+str(threshold_item)+'.txt','w')
    for line in open('blast_'+str(target_bin)+'_self_merged_'+str(threshold_item)+'.txt','r'):
        query=str(line).strip().split('\t')[0]
        subject=str(line).strip().split('\t')[1]
        similarity=float(str(line).strip().split('\t')[2])
        query_start=int(str(line).strip().split('\t')[6])
        query_end=int(str(line).strip().split('\t')[7])
        subject_start=int(str(line).strip().split('\t')[8])
        subject_end=int(str(line).strip().split('\t')[9])
        query_aligned_length=int(query_end)-int(query_start)+1
        identity=similarity*query_aligned_length/len(target_contig_seq[query])
        # if query == '9-6065':
        #     print(identity
        if identity >= 99:
            f.write(line)
            being_merged_contigs[query]=target_contig_seq[query]
            if subject not in being_merged_contigs2.keys():
                being_merged_contigs2[subject]=[]
                being_merged_contigs3[subject]={}
                being_merged_contigs2[subject].append(subject_start)
                being_merged_contigs2[subject].append(subject_end)
                being_merged_contigs3[subject][query]=1
            else:
                being_merged_contigs2[subject].append(subject_start)
                being_merged_contigs2[subject].append(subject_end)
                being_merged_contigs3[subject][query]=1
    f.close()

    revised_merged_seq, removed_merged_seq, elongated_seq_status, being_merged_contigs_copy={}, {}, {}, {}
    being_merged_contigs_copy=copy.deepcopy(being_merged_contigs)
    for item in being_merged_contigs2.keys():
        start=min(being_merged_contigs2[item])-1
        # print(start
        end=max(being_merged_contigs2[item])-1
        # print(end
        total_length=0
        for item2 in being_merged_contigs3[item].keys():
            total_length+=len(total_seq[item2])
        revised_length=int(end) - int(start) + 1

        if 100*total_length/revised_length >= 95:
            revised_merged_seq[item]=str(merged_bin_seq[item])[start:end]   
        else:
            removed_merged_seq[item]=str(merged_bin_seq[item])
            for item2 in being_merged_contigs3[item].keys():
                if item2 in being_merged_contigs_copy.keys():
                    del being_merged_contigs[item2]

        if 100*total_length/revised_length >= 90:
            if '90-100' not in elongated_seq_status.keys():
                elongated_seq_status['90-100']={}
            for item2 in being_merged_contigs3[item].keys():
                elongated_seq_status['90-100'][item2]=1
        elif 100*total_length/revised_length >= 80 and 100*total_length/revised_length < 90:
            if '80-90' not in elongated_seq_status.keys():
                elongated_seq_status['80-90']={}
            for item2 in being_merged_contigs3[item].keys():
                elongated_seq_status['80-90'][item2]=1
        elif 100*total_length/revised_length >= 70 and 100*total_length/revised_length < 80:
            if '70-80' not in elongated_seq_status.keys():
                elongated_seq_status['70-80']={}
            for item2 in being_merged_contigs3[item].keys():
                elongated_seq_status['70-80'][item2]=1
        elif 100*total_length/revised_length >= 60 and 100*total_length/revised_length < 70:
            if '60-70' not in elongated_seq_status.keys():
                elongated_seq_status['60-70']={}
            for item2 in being_merged_contigs3[item].keys():
                elongated_seq_status['60-70'][item2]=1
        else:
            if '<60' not in elongated_seq_status.keys():
                elongated_seq_status['<60']={}
            for item2 in being_merged_contigs3[item].keys():
                elongated_seq_status['<60'][item2]=1
        
    delete_orginal_contigs=[]
    f=open('Elongated_seq_status_'+str(threshold_item)+'_'+str(target_bin)+'.txt','w')
    for item in elongated_seq_status.keys():
        f.write(str(item)+'\t'+str(elongated_seq_status[item])+'\n')
        if item == '<60' or item == '60-70':
            for item2 in elongated_seq_status[item].keys():
                delete_orginal_contigs.append(item2)
    f.close()
    os.system('mv Elongated_seq_status_'+str(threshold_item)+'_'+str(target_bin)+'.txt '+folder_name)

    f=open('Revised_merged_bin_'+str(threshold_item)+'_'+str(target_bin),'w')
    for item in revised_merged_seq.keys():
        f.write('>'+item+'\n'+str(revised_merged_seq[item])+'\n')
    f.close()

    f=open('Removed_merged_seq_'+str(threshold_item)+'_'+str(target_bin),'w')
    for item in removed_merged_seq.keys():
        f.write('>'+item+'\n'+str(removed_merged_seq[item])+'\n')
    f.close()
    os.system('mv Removed_merged_seq_'+str(threshold_item)+'_'+str(target_bin)+' '+folder_name)

    name_lis=target_bin.split('.')
    name_lis.remove(name_lis[-1])
    new_name='.'.join(name_lis)+'.'+str(iteration_num)+'_'+str(threshold_item)+'.fa'
    f=open(new_name,'w')
    for record in SeqIO.parse(target_bin,'fasta'):
        if record.id not in being_merged_contigs.keys() and 'PB_' not in record.id:
            f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            
    for line in open('Revised_merged_bin_'+str(threshold_item)+'_'+str(target_bin),'r'):
        f.write(str(line))
    f.close()

    f=open('Deleted_potential_contaminated_contig_bin_'+str(threshold_item)+'_'+target_bin,'w')
    for record in SeqIO.parse(target_bin,'fasta'):
        if record.id not in being_merged_contigs.keys() and 'PB_' not in record.id and record.id not in delete_orginal_contigs:
            f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            
    for line in open('Revised_merged_bin_'+str(threshold_item)+'_'+str(target_bin),'r'):
        f.write(str(line))
    f.close()

    os.system('mv Revised_merged_bin_'+str(threshold_item)+'_'+str(target_bin)+' '+folder_name)
    os.system('mv Deleted_potential_contaminated_contig_bin_'+str(threshold_item)+'_'+target_bin+' '+folder_name)
    os.system('mv Filtrated_blast_'+str(target_bin)+'_self_merged_* '+folder_name)
    # os.system('rm *.nhd *.nin *.nsq *.nhr *.nog *.nsd *.nsi *.nhi')

def seq_merge(total_seq, query, subject, query_start, query_end,
              subject_start, subject_end, num_seq, target_bin):
    """
    Merge two contigs based on BLAST alignment coordinates.

    Parameters
    ----------
    total_seq : dict
        Mapping ``contig_id -> sequence`` for all contigs involved in
        the current OLC operation.
    query : str
        Contig ID of the query sequence in the BLAST alignment.
    subject : str
        Contig ID of the subject sequence in the BLAST alignment.
    query_start : int
        Start position of the alignment on the query contig (1-based).
    query_end : int
        End position of the alignment on the query contig (1-based).
    subject_start : int
        Start position of the alignment on the subject contig (1-based).
    subject_end : int
        End position of the alignment on the subject contig (1-based).
    num_seq : int
        Incremental index used when naming merged sequence files.
    target_bin : str
        Base name of the target bin currently being elongated.

    Returns
    -------
    dict
        Updated ``total_seq`` mapping including the new merged sequence.
    dict
        Mapping ``contig_id -> 1`` for contigs that have been merged.
    str
        File name of the merged sequence FASTA written to disk.
    """
    print('Merging sequences')
    processed_contigs={}
    delta_target_seq_alignment=abs(query_end-query_start)+1
    delta_vs_seq_alignment=abs(subject_end-subject_start)+1
    if delta_vs_seq_alignment == len(total_seq[subject]): ### In this case, subject sequence will totally covered by query seq
        merged_seq=total_seq[query]
    else:
        target_seq_start_nucl=total_seq[query][query_start-1]
        elong_seq_start_nucl=total_seq[subject][subject_start-1]
        if delta_target_seq_alignment == len(total_seq[query]):
            if subject_end > subject_start: ### Positive direction
                if subject_start == 1:
                    if target_seq_start_nucl == elong_seq_start_nucl:
                        merged_seq=total_seq[query]+total_seq[subject][subject_end:]
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=total_seq[query]+c_seq[subject_end:]
                else: ### Alignment start from Suject End
                    if target_seq_start_nucl == elong_seq_start_nucl:
                        merged_seq=total_seq[subject][:subject_start-1]+total_seq[query] ### Be careful about the subject start
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=c_seq[:subject_start-1]+total_seq[query]
            else: ### Reverse direction
                # r_seq=str(total_seq[subject])[::-1]
                # subject_reverse_start_position=len(subject)-subject_start
                if subject_end == 1:
                    # if total_seq[query][-1] == r_seq[-1]:
                    if total_seq[query][-1] == total_seq[subject][0]:
                        # merged_seq=r_seq[:subject_reverse_start_position]+total_seq[query]
                        merged_seq=total_seq[subject][subject_start:][::-1]+total_seq[query]
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1])).complement()
                        merged_seq=c_seq+total_seq[query]
                else: ### Alignment start from Suject End
                    if total_seq[query][0] == total_seq[subject][-1]:
                        # merged_seq=total_seq[query]+r_seq[:subject_start-1] ### ?
                        # merged_seq=total_seq[query]+r_seq[delta_vs_seq_alignment:]
                        merged_seq=total_seq[query]+total_seq[subject][:subject_end][::-1]
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][:subject_end][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][:subject_end][::-1])).complement()
                        # merged_seq=total_seq[query]+c_seq[:subject_start-1]
                        merged_seq=total_seq[query]+c_seq
        else: ### part of query and subject seqs merge together
            if subject_end > subject_start: 
                if subject_start == 1: ### Positive direction
                    if total_seq[query][query_start-1] == total_seq[subject][0]:
                        # merged_seq=total_seq[query][:query_end-1]+total_seq[subject][subject_end:]
                        merged_seq=total_seq[query]+total_seq[subject][subject_end:] #check
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        # merged_seq=total_seq[query][:query_end-1]+c_seq[subject_end:]
                        merged_seq=total_seq[query]+c_seq[subject_end:]
                else: ### Subject End
                    if total_seq[query][0] == total_seq[subject][subject_start-1]:
                        merged_seq=total_seq[subject][:subject_start-1]+total_seq[query] #check
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=c_seq[:subject_start-1]+total_seq[query]
            else: ### Reverse seq
                if subject_end == 1:
                    if total_seq[query][query_end-1] == total_seq[subject][0]:
                        merged_seq=total_seq[subject][subject_start:][::-1]+total_seq[query] ### Checked
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1])).complement()
                        merged_seq=c_seq+total_seq[query]
                else:
                    if total_seq[query][query_end-1] == total_seq[subject][subject_end-1]:
                        merged_seq=total_seq[query]+total_seq[subject][:subject_end-1][::-1] ### Checked
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][:subject_end-1][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][:subject_end-1][::-1])).complement()
                        merged_seq=total_seq[query]+c_seq

    if '--' in subject:
        subject1=subject.split('--')[0]
        subject2=subject.split('--')[-1]
        merged_subject=subject1+'--'+subject2
    else:
        merged_subject=subject

    if '--' in str(query):
        query_ori=str(query).split('--')[0]
        merged_file_name=target_bin+'--'+query_ori+'--'+merged_subject+'_'+str(num_seq)+'_merged_seq.txt'
        merged_seq_file=open(merged_file_name, 'w')
    else:
        merged_file_name=target_bin+'--'+query+'--'+merged_subject+'_'+str(num_seq)+'_merged_seq.txt'
        merged_seq_file=open(merged_file_name, 'w')

    merged_seq_file.write('>'+str(query)+'--'+str(subject)+'\n'+str(merged_seq)+'\n')
    merged_seq_file.close()
    total_seq[query+'--'+subject]=merged_seq
    processed_contigs[query]=1
    processed_contigs[subject]=1
    return total_seq, processed_contigs, merged_seq, merged_file_name

def elongation_main(blast_filtration_list, total_seq, target_contig_seq, aligned_len_cutoff, similarity_cutoff, num_threads, target_bin):
    pwd=os.getcwd()
    merge={}
    ### 1st run of merging
    processed_contigs, contigs_pool, contigs_query_pool={}, {}, {}
    contigs_pool['temp'], test, num_seq=1, 0, 0
    if len(blast_filtration_list) >= 2:
        while len(contigs_pool) != 0 and test != contigs_pool:
            test=contigs_pool
            n=0
            for item in blast_filtration_list:
                n+=1
                query=str(item).strip().split('\t')[0]
                subject=str(item).strip().split('\t')[1]
                contigs_pool[query]=1
                contigs_query_pool[query]=1
                contigs_pool[subject]=1
                if 'temp' in contigs_pool.keys():
                    del contigs_pool['temp']

            n=0
            for item in blast_filtration_list:
                n+=1
                query=str(item).strip().split('\t')[0]
                subject=str(item).strip().split('\t')[1]
                similarity=float(str(item).strip().split('\t')[2])
                aligned_length=int(str(item).strip().split('\t')[3])
                query_start=int(str(item).strip().split('\t')[6])
                query_end=int(str(item).strip().split('\t')[7])
                subject_start=int(str(item).strip().split('\t')[8])
                subject_end=int(str(item).strip().split('\t')[9])

                if n == 1:
                    num_seq+=1
                    A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, target_bin)
                    total_seq.update(A[0])
                    processed_contigs.update(A[1])
                # merge[str(query)+'--'+str(subject)]=str(A[2])
                    merge_seq_name=A[3]
                    del contigs_pool[query]
                    del contigs_pool[subject]
                    blast_filtration_list.remove(item)
            
            # print(contigs_pool
            # print(merge_seq_name
            num1, num2 = 0, 1

            while len(contigs_pool) != 0 and num1 != num2:
                num1=len(contigs_pool)
                for item in blast_filtration_list:
                    query=str(item).strip().split('\t')[0]
                    subject=str(item).strip().split('\t')[1]
                    similarity=float(str(item).strip().split('\t')[2])
                    aligned_length=int(str(item).strip().split('\t')[3])
                    query_start=int(str(item).strip().split('\t')[6])
                    query_end=int(str(item).strip().split('\t')[7])
                    subject_start=int(str(item).strip().split('\t')[8])
                    subject_end=int(str(item).strip().split('\t')[9])

                    if query in contigs_pool.keys() and subject in processed_contigs.keys():
                        # if '--' in str(subject):
                        #     subject_name=str(subject).split('--')[0]+'--'+str(subject).split('--')[-1]
                        # else:
                        #     subject_name=str(subject)

                        if '--' in str(query):
                            query_name=str(query).split('--')[0]+'--'+str(query).split('--')[-1]
                        else:
                            query_name=str(query)

                        f=open(query_name+'_seq.txt','w')
                        f.write('>'+str(query)+'\n'+str(total_seq[query])+'\n')
                        f.close()
                        positive_blast=elongation_sub_contig(merge_seq_name, query_name+'_seq.txt', n, aligned_length, similarity, num_threads)
                        if positive_blast != 0:
                            query_1=str(positive_blast).strip().split('\t')[0]
                            subject_1=str(positive_blast).strip().split('\t')[1]
                            query_start_1=int(str(positive_blast).strip().split('\t')[6])
                            query_end_1=int(str(positive_blast).strip().split('\t')[7])
                            subject_start_1=int(str(positive_blast).strip().split('\t')[8])
                            subject_end_1=int(str(positive_blast).strip().split('\t')[9])
                            num_seq+=1
                            A=seq_merge(total_seq, query_1, subject_1, query_start_1, query_end_1, subject_start_1, subject_end_1, num_seq, target_bin)
                            total_seq.update(A[0])
                            processed_contigs.update(A[1])
                            merge_seq=A[2]
                            merge_seq_name=A[3]
                            del contigs_pool[query]
                            blast_filtration_list.remove(item)
                    elif subject in contigs_pool.keys() and query in processed_contigs.keys():
                        f=open(subject+'_seq.txt','w')
                        f.write('>'+str(subject)+'\n'+str(total_seq[subject])+'\n')
                        f.close()
                        positive_blast=elongation_sub_contig(merge_seq_name, subject+'_seq.txt', n, aligned_len_cutoff, similarity_cutoff, num_threads)
                        if positive_blast != 0:
                            query_1=str(positive_blast).strip().split('\t')[0]
                            subject_1=str(positive_blast).strip().split('\t')[1]
                            query_start_1=int(str(positive_blast).strip().split('\t')[6])
                            query_end_1=int(str(positive_blast).strip().split('\t')[7])
                            subject_start_1=int(str(positive_blast).strip().split('\t')[8])
                            subject_end_1=int(str(positive_blast).strip().split('\t')[9])
                            num_seq+=1
                            A=seq_merge(total_seq, query_1, subject_1, query_start_1, query_end_1, subject_start_1, subject_end_1, num_seq, target_bin)
                            total_seq.update(A[0])
                            processed_contigs.update(A[1])
                            merge_seq=A[2]
                            merge_seq_name=A[3]
                            del contigs_pool[subject]
                            blast_filtration_list.remove(item)
                num2=len(contigs_pool)

def parse_checkm_1(test_bin_folder_checkm_containning_folder, pwd):
    # pwd=os.getcwd()
    # os.system('checkm lineage_wf -t 42 -x fa '+str(test_bin_folder)+' '+str(test_bin_folder)+'_checkm')
    os.chdir(test_bin_folder_checkm_containning_folder)
    # os.chdir(test_bin_folder+'_checkm/storage/')
    bin_checkm={}
    for root, dirs, files in os.walk(pwd+'/'+test_bin_folder_checkm_containning_folder):
        for file in files:
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:                   
                        binID=str(line).strip().split('\t')[0].strip()
                        bin_checkm[binID]={}

                        try:
                            genome_size=str(line).strip().split('\t')[8].strip()
                            completeness=str(line).strip().split('\t')[1].strip()
                            contamination=str(line).strip().split('\t')[2].strip()
                            N50=str(line).strip().split('\t')[6].strip()
                        except:
                            genome_size=str(line).strip().split('\t')[1].strip()
                            completeness=str(line).strip().split('\t')[2].strip()
                            contamination=str(line).strip().split('\t')[3].strip()
                            N50=str(line).strip().split('\t')[4].strip()

                        bin_checkm[binID]['N50']=int(eval(N50))
                        bin_checkm[binID]['Completeness']=float(eval(completeness))
                        bin_checkm[binID]['Genome size']=int(eval(genome_size))
                        bin_checkm[binID]['Contamination']=float(eval(contamination))
    os.chdir(pwd)
    return bin_checkm

def parse_checkm_2(test_bin_folder_checkm_containning_folder):
    pwd=os.getcwd()
    os.system('cp '+pwd+'/'+test_bin_folder_checkm_containning_folder+'/quality_report.tsv '+pwd+'/t_quality_report.tsv')
    bin_checkm, n = {}, 0
    for line in open('t_quality_report.tsv','r'): 
        n+=1
        if n >= 2:
            binID=str(line).strip().split('\t')[0].strip()
            bin_checkm[binID]={}

            try:
                genome_size=str(line).strip().split('\t')[8].strip()
                completeness=str(line).strip().split('\t')[1].strip()
                contamination=str(line).strip().split('\t')[2].strip()
                N50=str(line).strip().split('\t')[6].strip()
            except:
                genome_size=str(line).strip().split('\t')[1].strip()
                completeness=str(line).strip().split('\t')[2].strip()
                contamination=str(line).strip().split('\t')[3].strip()
                N50=str(line).strip().split('\t')[4].strip()

            bin_checkm[binID]['N50']=int(eval(N50))
            bin_checkm[binID]['Completeness']=float(eval(completeness))
            bin_checkm[binID]['Genome size']=int(eval(genome_size))
            bin_checkm[binID]['Contamination']=float(eval(contamination))
    os.system('rm t_quality_report.tsv')
    return bin_checkm

def bin_comparison(bin_checkm):
    pwd=os.getcwd()
    best_bin, best_bin_checkm, n={}, {}, 0
    for item in bin_checkm.keys():
        # print(item)
        n+=1
        if n == 1:
            selected_bin=item
        else:
            delta_cpn_ctn_query=1000000000*(bin_checkm[item]['Completeness']-5*bin_checkm[item]['Contamination'])
            delta_cpn_ctn_subject=1000000000*(bin_checkm[selected_bin]['Completeness']-5*bin_checkm[selected_bin]['Contamination'])
            if delta_cpn_ctn_query > delta_cpn_ctn_subject:
                selected_bin=item
            elif delta_cpn_ctn_query == delta_cpn_ctn_subject:
                selected_bin_num=selected_bin.count('.')
                item_num=item.count('.')
                if int(item_num) > int(selected_bin_num):
                    #    if bin_checkm[item]['Mean scaffold length'] > bin_checkm[selected_bin]['Mean scaffold length']:
                    selected_bin=item
                elif int(item_num) == int(selected_bin_num):
                    try:
                        threholds_num=float(item.split('_')[-1].split('.')[0])
                        if int(threholds_num) == 1:
                            selected_bin=item
                    # elif int(item_num) == int(selected_bin_num):
                    #     item_threhold=float(item.split('_genomes.')[1].split('.fa')[-2].split('_')[-1])
                    #     selected_bin_threhold=float(selected_bin.split('_genomes.')[1].split('.fa')[-2].split('_')[-1])
                    #     if item_threhold < selected_bin_threhold:
                    #         selected_bin=item
                    except:
                        xyzzz=0
    best_bin_cpn=bin_checkm[selected_bin]['Completeness']
    best_bin_ctn=bin_checkm[selected_bin]['Contamination']
    # best_bin_ml=bin_checkm[selected_bin]['Mean scaffold length']
    best_bin_checkm[selected_bin]=bin_checkm[selected_bin]

    #print(selected_bin, 'completeness', best_bin_cpn, 'contamination', best_bin_ctn, 'mean length', best_bin_ml
    print(selected_bin, 'completeness', best_bin_cpn, 'contamination', best_bin_ctn)
    return selected_bin, best_bin_checkm[selected_bin]
                
def OLC_elongation_main(target_bin, eliminated_bin, target_bin_checkm, iteration_num, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder):
    """
    Perform OLC-based elongation of a single target bin using an eliminated bin.

    Parameters
    ----------
    target_bin : str
        FASTA file of the target bin to be elongated.
    eliminated_bin : str
        FASTA file of the eliminated bin providing candidate contigs.
    target_bin_checkm : dict
        CheckM2 quality metrics for the target bin (completeness, contamination,
        etc.), propagated to newly produced bins.
    iteration_num : int
        Current OLC elongation iteration index.
    num : int
        Mode flag controlling whether contig-bin selection is performed
        (0) or an already selected bin is reused (1).
    aligned_len_cutoff : int
        Minimum BLAST alignment length to consider for elongation.
    similarity_cutoff : float
        Minimum percentage identity for BLAST alignments.
    num_threads : int
        Number of threads for BLAST and other CPU-intensive routines.
    pwd : str
        Working directory path.
    eliminated_bin_containing_folder : str
        Folder that contains the eliminated bin FASTA and depth files.

    Returns
    -------
    dict
        Mapping ``bin_name -> dict(checkm_metrics)`` for elongated bins
        retained after this iteration.
    dict
        Mapping ``bin_name -> dict(contig_id -> sequence)`` representing
        sequences of elongated bins.
    dict
        Mapping ``bin_name -> dict(contig_id -> int)`` indicating which
        contigs come from original bins and which were recruited.
    """
    # def OLC_elongation_main(target_bin, eliminated_bin, target_bin_checkm, iteration_num, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder, all_kmer):
    pwd=os.getcwd()
    clean_stale_olc_intermediates(pwd, target_bin)
    # try:
    merged_bin_recorded={}
    if num == 0:
        threshold=['1', '1.5', '3']
        # contig_bin_pool=elongate_contig_selector(eliminated_bin, threshold, pwd, eliminated_bin_containing_folder)
        # contig_bin_pool=elongate_contig_selector(eliminated_bin, threshold, pwd, eliminated_bin_containing_folder, all_kmer)
        contig_bin_pool=[eliminated_bin]
    elif num == 1:
        contig_bin_pool=[eliminated_bin]

    try:
        os.mkdir(target_bin+'_merged')
    except:
        print(target_bin+'_merged folder exist')
        # print(target_bin+'_merged folder exist. Re-create the folder')
        # os.system('rm -rf '+target_bin+'_merged')
        # os.mkdir(target_bin+'_merged')

    xyz, bin_checkm, bin_contigs, bin_contigs2, mn, recorded_bins =0, {}, {}, {}, 0, {}
    name_lis=target_bin.split('.')
    name_lis.remove(name_lis[-1])
    split_name='.'.join(name_lis)
    bin_checkm[split_name]=target_bin_checkm.copy()
    
    for bin_item in contig_bin_pool:
        merged_files = {}
        try:
            xyz+=1
            bin_item_head=bin_item.split('.')[0]
            target_bin_head=target_bin.split('.')[0]
            blast_name=str(target_bin_head)+'_'+str(bin_item_head)+'_'+str(xyz)+'.txt'

            if num == 0:
                threshold_item=str(bin_item).split('_')[0]
            elif num == 1:
                threshold_item='re'
            
            folder_name='Merged_seqs_'+target_bin+'_'+bin_item
            try:
                os.mkdir(folder_name)
            except:
                print(folder_name+' exists')

            print('Processing '+target_bin+' with '+bin_item)
            A=record_seq(target_bin, bin_item)
            target_contig_seq, target_contig_len, vs_contig_seq, vs_contig_len, Merged_seq, total_seq = A[0], A[1], A[2], A[3], A[4], A[5]
            merged_bin_recorded.update(A[6])
            blast_group=blast_1(target_bin, bin_item, target_contig_seq, target_contig_len, vs_contig_seq, vs_contig_len, aligned_len_cutoff, similarity_cutoff, num_threads, blast_name, folder_name)
            # print(len(blast_group))
            n=0
            for item in blast_group.keys():
                n+=1
                print('--------------------------')
                print('Processing contig group', item, 'in', str(target_bin), 'with', str(bin_item))
                # if n == 78:
                elongation_main(blast_group[item], total_seq, target_contig_seq, aligned_len_cutoff, similarity_cutoff, num_threads, target_bin)

            name, name_select = {}, {}
            for root, dirs, files in os.walk(pwd):
                for file in files:
                    if str(target_bin) in file and '_merged_seq.txt' in file and 'blast_' not in file:
                        merged_files[file]=0
                        orig=str(file).split('--')[1]
                        for record in SeqIO.parse(file, 'fasta'):
                            length=len(record.id)
                        if orig not in name.keys():
                            name[orig]=length
                            name_select[orig]=file
                        elif length > name[orig]:
                            name[orig]=length
                            name_select[orig]=file
        
            selected=[]
            for item in name_select.keys():
                selected.append(name_select[item])

            for root, dirs, files in os.walk(pwd):
                for file in files:
                    if str(target_bin) in file and '_merged_seq.txt' in file and file not in selected:
                        os.system('rm '+file)
                        del merged_files[file]

            if len(merged_files) != 0:
                fxy=open('Merged_'+str(threshold_item)+'_'+target_bin, 'w')
                for item in merged_files.keys():
                    for line in open(item,'r'):
                        fxy.write(line)
                fxy.close()
                # os.system('cat *_merged_seq.txt > Merged_'+str(threshold_item)+'_'+target_bin)

                blast_2(target_bin, target_contig_seq, 'Merged_'+str(threshold_item)+'_'+target_bin, total_seq, threshold_item, iteration_num, num_threads, folder_name)
                os.system('mv Merged_'+str(threshold_item)+'_'+target_bin+' '+folder_name)
                os.system('tar zcvf '+folder_name+'.tar.gz '+folder_name)
                os.system('rm -rf '+folder_name)
                name_lis=target_bin.split('.')
                name_lis.remove(name_lis[-1])
                split_name='.'.join(name_lis)
                new_name='.'.join(name_lis)+'.'+str(iteration_num)+'_'+str(threshold_item)+'.fa'
                bin_checkm[split_name]=target_bin_checkm.copy()
                
                mn+=1
                bin_contigs[mn]={}
                bin_contigs2[mn]=str(new_name)
                for record in SeqIO.parse(new_name, 'fasta'):
                    bin_contigs[mn][record.id]=record.seq
                recorded_bins[new_name]=0
                os.system('mv '+folder_name+'.tar.gz '+target_bin+'_merged')
            os.system('rm -rf '+target_bin+'_split_blast_output')
            os.system('rm -rf '+folder_name)
        except:
            xyzzz=0

    del_bin, record_del_index={}, {}
    for i in range(2, len(bin_contigs)+1):
        if bin_contigs[i-1]==bin_contigs[i]:
            record_del_index[i-1]=''
            # del bin_contigs[i-1]
            del_bin[bin_contigs2[i-1]]=''
            # del bin_contigs2[i-1]
    
    for i in record_del_index.keys():
        del bin_contigs[i]
        del bin_contigs2[i]

    recorded_bin2={}
    for bins in recorded_bins.keys():
        if bins in del_bin.keys(): 
            os.system('rm '+bins)
        else:
            recorded_bin2[bins]=0

    # try:
    #     if len(merged_files) != 0 and len(recorded_bin2) != 0:
    #         for bins in recorded_bin2.keys():
    #             os.system('mv '+str(bins)+' '+target_bin+'_merged')
    #         os.system('cp '+target_bin+' '+target_bin+'_merged')

    #         os.system('checkm2 predict -t '+str(num_threads)+' -i '+target_bin+'_merged -x fa -o '+str(target_bin)+'_checkm')
    #         test_checkm=parse_checkm_2(target_bin+'_checkm')
    #         # print(test_checkm)
    #         bin_checkm.update(test_checkm)
    # except:
    #     xxyyzz=0

    try:
        if len(merged_files) != 0 and len(recorded_bin2) != 0:
            xline2, xt=0, 0
            for line in open(target_bin,'r'):
                xline2+=1
            os.system('cp '+target_bin+' '+target_bin+'_merged')

            for bins in recorded_bin2.keys():
                xline=0
                for line in open(bins,'r'):
                    xline+=1

                if xline < xline2:
                    os.system('mv '+str(bins)+' '+target_bin+'_merged')
                    xt+=1

            if xt > 0:
                test_checkm=evaluate_bins_checkm2(target_bin+'_merged', num_threads)
                if test_checkm is None:
                    os.system('checkm2 predict -t '+str(positive_thread_count(num_threads))+' -i '+target_bin+'_merged -x fa -o '+str(target_bin)+'_checkm --force')
                    test_checkm=parse_checkm_2(target_bin+'_checkm')
                if test_checkm:
                    bin_checkm.update(test_checkm)
    except:
        xxyyzz=0

    best_bin=bin_comparison(bin_checkm)
    select_bin_checkm=best_bin[1]
    if best_bin[0] != split_name:
        best_bin_name=best_bin[0]+'.fa'
        # print(best_bin_name
        os.system('mv '+pwd+'/'+target_bin+'_merged/'+best_bin_name+' '+pwd)
    else:
        best_bin_name=str(target_bin)
        # name_lis=str(target_bin).split('.')
        # name_lis.remove(name_lis[-1])
        # target_bin_rename='.'.join(name_lis)+'.fa'
        # os.system('cp '+str(target_bin)+' '+str(target_bin_rename))
    # except:
    #     best_bin_name=str(target_bin)
    #     name_lis=target_bin.split('.')
    #     name_lis.remove(name_lis[-1])
    #     split_name='.'.join(name_lis)
    #     select_bin_checkm={}
    #     select_bin_checkm[split_name]=target_bin_checkm.copy()
    #     print('Skipping eliminated bin '+str(eliminated_bin)+', error in recording '+str(eliminated_bin)+' information')
    #     fbin_record_error=open('Bin_record_error.txt','a')
    #     fbin_record_error.write('Recorded '+str(eliminated_bin)+' error'+'\n')
    #     fbin_record_error.close()
    return best_bin_name, select_bin_checkm, merged_bin_recorded

def merge(target_bin_folder, target_bin, eliminated_bin_list, target_bin_checkm, num, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder):
    # def merge(target_bin_folder, target_bin, eliminated_bin_list, target_bin_checkm, num, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder, all_kmer):
    #pwd=os.getcwd()
    n=0
    for item in eliminated_bin_list:
        merged_bin_recorded={}
        print('Utilization of', item, ' to elongate', target_bin)
        if num == 0:
            qz_list=str(item).split('_genomes.')[0].split('_')
            qz_list.remove(qz_list[-1])
            qz_list.remove(qz_list[-1])
            eliminated_bin_containing_folder='_'.join(qz_list)+'_BestBinsSet'
            os.system('cp '+pwd+'/'+eliminated_bin_containing_folder+'/'+str(item)+' '+pwd)
            if len(mod_bin_folder) != 0:
                for folder_name in mod_bin_folder:
                    os.system('cp '+pwd+'/'+folder_name+'/'+str(item)+' '+pwd)
        elif num == 1:
            if '_re-assembly_contigs.fa' in item:
                os.system('cp '+pwd+'/'+str(target_bin_folder)+'/'+item+' '+pwd)
            else:
                mod_folder=target_bin_folder.split('_re-assembly')[0]+'_mod'
                os.system('cp '+pwd+'/'+str(mod_folder)+'/'+item+' '+pwd)

        n+=1
        if n == 1:
            best_bin=OLC_elongation_main(target_bin, item, target_bin_checkm, n, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder)
            # best_bin=OLC_elongation_main(target_bin, item, target_bin_checkm, n, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder, all_kmer)
        else:
            best_bin=OLC_elongation_main(best_bin[0], item, best_bin[1], n, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder)
            # best_bin=OLC_elongation_main(best_bin[0], item, best_bin[1], n, num, aligned_len_cutoff, similarity_cutoff, num_threads, pwd, eliminated_bin_containing_folder, all_kmer)
        os.system('rm '+str(item))

        merged_bin_recorded.update(best_bin[2])
        for item2 in merged_bin_recorded.keys():
            os.system('rm '+str(item2))
        print('-------------------')
    return best_bin, item

def finding_similar_bins(target_bin_folder, bin_comparison_folder):
    pwd=os.getcwd()
    print('Forming similar bins file')
    binset, bestbinset_sim_bin, m = {}, {}, 0
    os.chdir(pwd+'/'+target_bin_folder)
    for root, dirs, files in os.walk(pwd+'/'+target_bin_folder):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('.')[-1]
                if 'fa' in hz:
                    name_lis=str(file).split('.')
                    name_lis.remove(name_lis[-1])
                    new_name='.'.join(name_lis)+'.'
                    binset[file]=new_name
            
            if 'Selected_bins_best_binset.txt' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2 and '---' in line:
                        x=0
                        bin1=str(line).strip().split('\t')[1].split('---')[0]
                        bin2=str(line).strip().split('\t')[1].split('---')[1]
                        bin1_ctn=float(eval(str(line).strip().split('Contamination')[1].split(':')[1].split('}')[0].strip()))
                        bin2_ctn=float(eval(str(line).strip().split('Contamination')[2].split(':')[1].split('}')[0].strip()))
                        if len(bestbinset_sim_bin) != 0:
                            for item in bestbinset_sim_bin.keys():
                                if bin1 in bestbinset_sim_bin[item].keys() and bin2 not in bestbinset_sim_bin[item].keys():
                                    bestbinset_sim_bin[item][bin2]=bin2_ctn
                                    x=1
                                elif bin2 in bestbinset_sim_bin[item].keys() and bin1 not in bestbinset_sim_bin[item].keys():
                                    bestbinset_sim_bin[item][bin1]=bin1_ctn
                                    x=1
                                elif bin1 in bestbinset_sim_bin[item].keys() and bin2 in bestbinset_sim_bin[item].keys():
                                    x=1
                                
                            if x == 0:
                                m+=1
                                bestbinset_sim_bin[m]={}
                                bestbinset_sim_bin[m][bin1]=bin1_ctn
                                bestbinset_sim_bin[m][bin2]=bin2_ctn
                        else:
                            m+=1
                            bestbinset_sim_bin[m]={}
                            bestbinset_sim_bin[m][bin1]=bin1_ctn
                            bestbinset_sim_bin[m][bin2]=bin2_ctn

    os.chdir(pwd+'/'+bin_comparison_folder)
    for root, dirs, files in os.walk(pwd+'/'+bin_comparison_folder):
        for file in files:
            if 'Selected_bins_' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2 and '---' in line:
                        x=0
                        bin1=str(line).strip().split('\t')[1].split('---')[0].strip()
                        bin2=str(line).strip().split('\t')[1].split('---')[1].strip()
                        bin1_ctn=float(eval(str(line).strip().split('Contamination')[1].split(':')[1].split('}')[0].strip()))
                        bin2_ctn=float(eval(str(line).strip().split('Contamination')[2].split(':')[1].split('}')[0].strip()))
                        if len(bestbinset_sim_bin) != 0:
                            for item in bestbinset_sim_bin.keys():
                                if bin1 in bestbinset_sim_bin[item].keys() and bin2 not in bestbinset_sim_bin[item].keys():
                                    bestbinset_sim_bin[item][bin2]=bin2_ctn
                                    x=1
                                elif bin2 in bestbinset_sim_bin[item].keys() and bin1 not in bestbinset_sim_bin[item].keys():
                                    bestbinset_sim_bin[item][bin1]=bin1_ctn
                                    x=1
                                elif bin1 in bestbinset_sim_bin[item].keys() and bin2 in bestbinset_sim_bin[item].keys():
                                    x=1
                                
                            if x == 0:
                                m+=1
                                bestbinset_sim_bin[m]={}
                                bestbinset_sim_bin[m][bin1]=bin1_ctn
                                bestbinset_sim_bin[m][bin2]=bin2_ctn
                        else:
                            m+=1
                            bestbinset_sim_bin[m]={}
                            bestbinset_sim_bin[m][bin1]=bin1_ctn
                            bestbinset_sim_bin[m][bin2]=bin2_ctn

    os.chdir(pwd+'/'+target_bin_folder)
    num, bestbinset_sim_bin2=0, {}
    for root, dirs, files in os.walk(pwd+'/'+target_bin_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                num+=1
                if '_genomes.' in file:
                    qz=str(file).split('_genomes.')[0]
                    bin_num=str(file).split('_genomes.')[1].split('.')[0]
                    bin_id=qz+'_genomes.'+bin_num
                for item in bestbinset_sim_bin.keys():
                    for item2 in bestbinset_sim_bin[item].keys():
                        if bin_id +'.' in item2:
                            #del bestbinset_sim_bin[item][item2]
                            bestbinset_sim_bin2[file]=[]
                            for itemx in bestbinset_sim_bin[item].keys():
                                itemx_name_list=itemx.split('.')
                                itemx_name_list.pop()
                                itemx_name='.'.join(itemx_name_list)
                                file_name_list=str(file).split('.')
                                file_name_list.pop()
                                file_name='.'.join(file_name_list)
                                if itemx_name != file_name:
                                    if bestbinset_sim_bin[item][itemx] <= 5:
                                        bestbinset_sim_bin2[file].append(itemx)
                            #del bestbinset_sim_bin[item]

    os.chdir(pwd)
    del_bin={}
    for item in bestbinset_sim_bin2.keys():
        for i in range(0,len(bestbinset_sim_bin2[item])):
            if item == bestbinset_sim_bin2[item][i]:
                del_bin[item]=i

    for item in del_bin.keys():
        del bestbinset_sim_bin2[item][del_bin[item]]

    f=open('Similar_bins.txt','w')
    del_bin={}
    for item in bestbinset_sim_bin2.keys():
        if len(bestbinset_sim_bin2[item]) != 0:
            f.write(str(item)+'\t'+str(bestbinset_sim_bin2[item])+'\n')
        else:
            del_bin[item]=''
    f.close()

    for item in del_bin.keys():
        del bestbinset_sim_bin2[item]

    return bestbinset_sim_bin2

def mapping(total_fa, datasets_list, fq, num_threads):
    """
    Map short reads to a combined bin reference and prepare bin-specific FASTQ.

    Parameters
    ----------
    total_fa : str
        Path to the concatenated FASTA file containing all bin sequences.
    datasets_list : dict
        Mapping ``sample_id -> [r1_fastq, r2_fastq]`` with paired-end reads.
    fq : dict
        Mapping ``bin_id -> dict(read_id -> int)`` used to track mates
        when writing bin-specific FASTQ files.
    num_threads : int
        Number of Bowtie2 threads to use.

    Returns
    -------
    None
        Writes bin-specific FASTQ files to the current working directory.
    """
    os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))
    from S9_Reassembly_10262023 import parse_sam
    pair = {bin_id: {} for bin_id in fq}
    n = 0
    for item in datasets_list.keys():
        n+=1
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets_list[item][0])+' -2 '+str(datasets_list[item][1])+' -S '+str(item)+'.sam -q --no-unal')
        parse_sam(str(item)+'.sam', fq, pair, n)

def reassembly_paired_bins(target_bin_folder):
    """
    Group reassembled bins that originate from the same original bin.

    Parameters
    ----------
    target_bin_folder : str
        Folder containing original and reassembled bin FASTA files.

    Returns
    -------
    dict
        Mapping ``bin_filename -> list(similar_bin_filenames)`` indicating
        which reassembled bins should be compared together in downstream
        OLC steps.
    """
    pwd=os.getcwd()
    print('Forming similar bins file')
    binset, bestbinset_sim_bin, total_seq, m = {}, {}, {}, 0
    os.chdir(pwd+'/'+target_bin_folder)

    reassembly_bins=[]
    for root, dirs, files in os.walk(pwd+'/re-assembly_binset/'):
        for file in files:
            if '.fa' in file:
                reassembly_bins.append(file)
            else:
                continue

    for root, dirs, files in os.walk(pwd+'/'+target_bin_folder):
        for file in files:
            bestbinset_sim_bin[file]=[]
            if '_re-assembly_contigs.fa' in file:
                org_bin=str(file).split('_')[0]+'.fa'
                bestbinset_sim_bin[file].append(org_bin)
                if 'SPAdes' in file:
                    IDBA_bin=org_bin.split('.fa')[0]+'_IDBA_re-assembly_contigs.fa'
                    megahit_bin=org_bin.split('.fa')[0]+'_megahit_re-assembly_contigs.fa'
                    if IDBA_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(IDBA_bin)
                    
                    if megahit_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(megahit_bin)

                elif 'IDBA' in file:
                    SPAdes_bin=org_bin.split('.fa')[0]+'_SPAdes_re-assembly_contigs.fa'
                    megahit_bin=org_bin.split('.fa')[0]+'_megahit_re-assembly_contigs.fa'
                    if SPAdes_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(SPAdes_bin)

                    if megahit_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(megahit_bin)

                elif 'megahit' in file:
                    IDBA_bin=org_bin.split('.fa')[0]+'_IDBA_re-assembly_contigs.fa'
                    SPAdes_bin=org_bin.split('.fa')[0]+'_SPAdes_re-assembly_contigs.fa'
                    if IDBA_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(IDBA_bin)

                    if SPAdes_bin in reassembly_bins:
                        bestbinset_sim_bin[file].append(SPAdes_bin)

                for record in SeqIO.parse(file, 'fasta'):
                    total_seq[record.id]=record.seq
            elif '.fa' in file:
                IDBA_bin=str(file).split('.fa')[0]+'_IDBA_re-assembly_contigs.fa'
                SPAdes_bin=str(file).split('.fa')[0]+'_SPAdes_re-assembly_contigs.fa'
                megahit_bin=str(file).split('.fa')[0]+'_megahit_re-assembly_contigs.fa'
                xxx=0
                if IDBA_bin in reassembly_bins:
                    bestbinset_sim_bin[file].append(IDBA_bin)
                    xxx=1
                
                if SPAdes_bin in reassembly_bins:
                    bestbinset_sim_bin[file].append(SPAdes_bin)
                    xxx=1

                if megahit_bin in reassembly_bins:
                    bestbinset_sim_bin[file].append(megahit_bin)
                    xxx=1

                if xxx == 1:
                    for record in SeqIO.parse(file, 'fasta'):
                        total_seq[record.id]=record.seq
            else:
                continue

    os.chdir(pwd)
    f=open('Similar_bins.txt','w')
    for item in bestbinset_sim_bin.keys():
        f.write(str(item)+'\t'+str(bestbinset_sim_bin[item])+'\n')
    f.close()

    print('Obtaining depth file')
    f=open('Total_contigs_after_OLC_reassembly.fa','w')
    for item in total_seq.keys():
        f.write('>'+str(item)+'\n'+str(total_seq[item])+'\n')
    f.close()

    return bestbinset_sim_bin

def mul_threads(item, bestbinset_sim_bin, bestbinset_checkm, step, pwd, aligned_len_cutoff, similarity_cutoff, mod_bin_folder, target_bin_folder, num_threads):
    """
    Worker function for parallel OLC processing of a single bin.

    Parameters
    ----------
    item : str
        File name of the current bin to be processed.
    bestbinset_sim_bin : dict
        Mapping ``bin_filename -> list(similar_bin_filenames)`` generated
        by :func:`reassembly_paired_bins`.
    bestbinset_checkm : dict
        Mapping ``bin_id -> dict(checkm_metrics)`` with quality statistics
        for each bin.
    step : {'assemblies_OLC', 'OLC_after_reassembly'}
        Pipeline step indicating whether OLC is being run on assemblies
        or on reassembled bins.
    pwd : str
        Working directory path.
    aligned_len_cutoff : int
        Minimum alignment length for BLAST-based comparisons.
    similarity_cutoff : float
        Minimum percentage identity for BLAST-based comparisons.
    mod_bin_folder : str
        Folder containing the modified binset used for OLC.
    target_bin_folder : str
        Folder containing the bins being processed.
    num_threads : int
        Number of worker threads for BLAST and other operations.

    Returns
    -------
    dict
        Mapping ``selected_bin_filename -> selected_bin_id`` for the best
        bin chosen for this item.
    dict
        Mapping ``selected_bin_id -> dict(checkm_metrics)`` with updated
        CheckM statistics.
    dict
        Mapping ``bin_filename -> 0`` for bins that should be removed.
    """
    total_selected_bin, total_selected_bin_checkm, remove_bin={}, {}, {}
    checkm_name_list=item.split('.')
    checkm_name_list.remove(checkm_name_list[-1])
    checkm_name='.'.join(checkm_name_list)
    target_bin_rename=checkm_name+'.fa'
    os.system('cp '+pwd+'/'+target_bin_folder+'/'+str(item)+' '+pwd+'/'+target_bin_rename)
    os.system('mv '+pwd+'/'+target_bin_folder+'/'+str(item)+' '+pwd+'/'+target_bin_folder+'_temp')
    if step == 'assemblies_OLC':
        best_bin, rm_bin=merge(target_bin_folder, target_bin_rename, bestbinset_sim_bin[item], bestbinset_checkm[checkm_name], 0, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder)
        # best_bin, rm_bin=merge(target_bin_folder, target_bin_rename, bestbinset_sim_bin[item], bestbinset_checkm[checkm_name], 0, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder, all_kmer)
    elif step == 'OLC_after_reassembly':
        best_bin, rm_bin=merge(target_bin_folder, target_bin_rename, bestbinset_sim_bin[item], bestbinset_checkm[checkm_name], 1, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder)
        # best_bin, rm_bin=merge(target_bin_folder, target_bin_rename, bestbinset_sim_bin[item], bestbinset_checkm[checkm_name], 1, pwd, aligned_len_cutoff, similarity_cutoff, num_threads, mod_bin_folder, all_kmer)
    ###
    remove_bin[rm_bin]=0

    try:
        selected_bin=best_bin[0]
        selected_bin_checkm=best_bin[1]
        name_lis=str(selected_bin).split('.')
        name_lis.remove(name_lis[-1])
        selected_bin_checkm_name='.'.join(name_lis) ### checkm name is different from actual name
        total_selected_bin_checkm[selected_bin_checkm_name]=selected_bin_checkm
        print('Selected', selected_bin_checkm_name)
        total_selected_bin[selected_bin_checkm_name]=1
    except:
        selected_bin=item
        name_lis=str(item).split('.')
        name_lis.remove(name_lis[-1])
        selected_bin_checkm_name='.'.join(name_lis)
        total_selected_bin_checkm[selected_bin_checkm_name]=bestbinset_checkm[selected_bin_checkm_name]
    os.system('cp '+selected_bin+' '+target_bin_folder+'_OLC')

    return total_selected_bin, total_selected_bin_checkm, remove_bin

def OLC_main(target_bin_folder, step, bin_comparison_folder,
             aligned_len_cutoff, similarity_cutoff, coverage_extension,
             num_threads, mod_bin_folder):
    """
    High-level driver for the S8 OLC reassembly step (CheckM2 branch).

    Parameters
    ----------
    target_bin_folder : str
        Folder containing the bins to be processed by OLC.
    step : {'assemblies_OLC', 'OLC_after_reassembly'}
        Pipeline step indicating whether OLC is applied to assemblies
        or to previously reassembled bins.
    bin_comparison_folder : str
        Folder containing CheckM-derived bin comparison results.
    aligned_len_cutoff : int
        Minimum alignment length for sequence comparisons.
    similarity_cutoff : float
        Minimum percentage identity for sequence comparisons.
    coverage_extension : float
        Minimum coverage extension threshold used when assessing bin
        improvements after OLC.
    num_threads : int
        Number of parallel worker processes to spawn.
    mod_bin_folder : str
        Folder containing the modified binset used for OLC.

    Returns
    -------
    None
        Writes updated bins, CheckM reports, and comparison summaries
        to disk; the function coordinates per-bin workers and does not
        return Python objects.
    """
    pwd=os.getcwd()
    ferror=open('OLC_merged_error_blast_results.txt','w')
    ferror.close()
    fblasterror=open('BLAST_output_error.txt','w')
    fblasterror.close()
    fbin_record_error=open('Bin_record_error.txt','w')
    fbin_record_error.close()
    accomplished_bins={}
    try:
        os.mkdir(target_bin_folder+'_OLC')
    except:
        print(target_bin_folder+'_OLC Exists. Reading bins in the folder.')
        # os.system('rm -rf '+target_bin_folder+'_OLC')
        # os.mkdir(target_bin_folder+'_OLC')
        os.chdir(target_bin_folder+'_OLC')
        x,y=0,0
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_OLC'):
            for file in files:
                try:
                    hz=str(file).split('.')[-1]
                    if 'fa' in hz:
                        if 'genomes' in file:
                            fasta_path=os.path.join(pwd, target_bin_folder+'_OLC', file)
                            if not fasta_has_sequence(fasta_path):
                                print('Ignoring incomplete OLC FASTA '+str(file))
                                continue
                            file_name_qz=str(file).split('_genomes.')[0]
                            file_name_num=str(file).split('_genomes.')[1].split('.')[0]
                            file_name=file_name_qz+'_genomes.'+str(file_name_num)
                            accomplished_bins[file_name+'.fa']=0
                            accomplished_bins[file_name+'.fasta']=0
                        x+=1
                except:
                    y=0
        print(str(x),'bin(s) finished OLC')
        os.chdir(pwd)

    try:
        os.mkdir(target_bin_folder+'_temp')
    except:
        # os.system('rm -rf '+target_bin_folder+'_temp')
        # os.mkdir(target_bin_folder+'_temp')
        os.chdir(target_bin_folder+'_temp')
        x,y=0,0
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_temp'):
            for file in files:
                try:
                    hz=str(file).split('.')[-1]
                    if 'fa' in hz:
                        os.system('mv '+str(file)+' '+pwd+'/'+target_bin_folder)
                        x+=1
                except:
                    y=0
        print(str(x),' were moved back to the target bin folder')
        os.chdir(pwd)

    bestbinset_checkm_org={}
    bestbinset_checkm=parse_checkm_1(target_bin_folder, pwd)
    bestbinset_checkm_org.update(bestbinset_checkm)
    # print(bestbinset_checkm

    if step == 'assemblies_OLC':
        bestbinset_sim_bin=finding_similar_bins(target_bin_folder, bin_comparison_folder)
    elif step == 'OLC_after_reassembly':
        bestbinset_sim_bin=reassembly_paired_bins(target_bin_folder)
    else:
        print('step parameter error!')
    
    bestbinset_sim_bin2={}
    bestbinset_sim_bin2=copy.deepcopy(bestbinset_sim_bin)
    for item in bestbinset_sim_bin2.keys():
        if len(bestbinset_sim_bin2[item]) == 0:
            del bestbinset_sim_bin[item]

    # print('Reading kmer files')
    # kmer_list=glob.glob(r'*.kmer.txt')
    # all_kmer={}
    # for item in kmer_list:
    #     n=0
    #     for line in open(item, 'r'):
    #         n+=1
    #         if n >= 2:
    #             contig=str(line).strip().split('\t')[0]
    #             all_kmer[contig]=str(line).strip()

    pending_count=0
    for item in bestbinset_sim_bin:
        if item not in accomplished_bins.keys():
            pending_count+=1
    xyzzz=0
    if pending_count >= 2:
        pool_size=min(positive_thread_count(num_threads), pending_count)
        threads_each=worker_thread_budget(num_threads, pool_size)
        print('Multiple threads started: '+str(pool_size)+' workers, '+str(threads_each)+' threads each')
        pool=Pool(processes=pool_size)
        result, remove_bin = {}, {}
        for item in bestbinset_sim_bin:
            if item not in accomplished_bins.keys():
                print('Processing '+str(item))
                # result[item]=pool.apply_async(mul_threads, args=(item, bestbinset_sim_bin, bestbinset_checkm, step, pwd, aligned_len_cutoff, similarity_cutoff, mod_bin_folder, target_bin_folder, num_threads, all_kmer))
                result[item]=pool.apply_async(mul_threads, args=(item, bestbinset_sim_bin, bestbinset_checkm, step, pwd, aligned_len_cutoff, similarity_cutoff, mod_bin_folder, target_bin_folder, threads_each,))
                # total_selected_bin, total_selected_bin_checkm, remove_bin=mul_threads(item, bestbinset_sim_bin, bestbinset_checkm, step, pwd, aligned_len_cutoff, similarity_cutoff, mod_bin_folder, target_bin_folder, num_threads, all_kmer)
        pool.close()
        pool.join()

        result2={}
        for item in result:
            result2[item]=result[item].get()

        total_selected_bin, total_selected_bin_checkm, remove_bin={}, {}, {}
        for item in result2.keys():
            total_selected_bin.update(result2[item][0])
            total_selected_bin_checkm.update(result2[item][1])
            remove_bin.update(result2[item][2])
        
        for item in remove_bin.keys():
            os.system('rm '+item)
        print('Multiple threads ended')
        xyzzz=1
    elif len(bestbinset_sim_bin) != 0:
        remove_bin = {}
        for item in bestbinset_sim_bin:
            if item not in accomplished_bins.keys():
                print('Processing '+str(item))
                total_selected_bin, total_selected_bin_checkm, remove_bin=mul_threads(item, bestbinset_sim_bin, bestbinset_checkm, step, pwd, aligned_len_cutoff, similarity_cutoff, mod_bin_folder, target_bin_folder, num_threads)

        for item in remove_bin.keys():
            os.system('rm '+item)
        xyzzz=1
    else:
        print('There is no qualified redundent bin')

    if xyzzz == 1:
        OLC_bins={}
        os.chdir(pwd+'/'+target_bin_folder+'_OLC/')
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_OLC/'):
            for file in files:
                hz=str(file).split('.')[-1]
                if 'fa' in  hz:
                    n, x=0, 0
                    for line in open(file,'r'):
                        n+=1
                        if n == 2:
                            x=1
                            break
                    
                    if x == 1:
                        qz=str(file).split('_genomes.')[0]
                        bin_num=str(file).split('_genomes.')[1].split('.')[0]
                        target_bin_rename=qz+'_genomes.'+str(bin_num)+'.fa'
                        target_bin_rename2=qz+'_genomes.'+str(bin_num)+'.fasta'
                        OLC_bins[target_bin_rename]=0
                        OLC_bins[target_bin_rename2]=0


        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_temp'):
            for file in files:
                if '.fa' in  file:
                    os.system('mv '+pwd+'/'+target_bin_folder+'_temp/'+file+' '+pwd+'/'+target_bin_folder)

        os.chdir(pwd+'/'+target_bin_folder)
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder):
            for file in files:
                hz=str(file).split('.')[-1]
                if 'fa' in  hz:
                    if file not in OLC_bins.keys():
                        name_lis=str(file).split('.')
                        name_lis.remove(name_lis[-1])
                        target_bin_checkm_name='.'.join(name_lis)
                        target_bin_rename=target_bin_checkm_name+'.fa'
                        print('Moved', file, 'to the selected bins pool')
                        os.system('cp '+file+' '+pwd+'/'+target_bin_folder+'_OLC/'+target_bin_rename)
                        total_selected_bin_checkm[target_bin_checkm_name]=bestbinset_checkm[target_bin_checkm_name]

        # print(str(len(total_selected_bin_checkm)), 'recored in quality tsv file')
        os.chdir(pwd)
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+target_bin_folder+'_OLC -x fa -o '+target_bin_folder+'_OLC_checkm')
        os.system('cp '+pwd+'/'+target_bin_folder+'_OLC_checkm/quality_report.tsv '+pwd+'/'+target_bin_folder+'_OLC/OLC_quality_report.tsv')

        os.chdir(pwd+'/'+target_bin_folder+'_OLC')
        max_num=0
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_OLC'):
            for file in files:            
                hz=str(file).split('.')[-1]
                if 'fa' in hz:
                    iteration_num=int(str(file).split('_')[0])
                    if iteration_num > max_num:
                        max_num=iteration_num
        
        os.chdir(pwd+'/'+target_bin_folder)
        assemblies_name={}
        for root, dirs, files in os.walk(pwd+'/'+target_bin_folder+'_OLC'):
            for file in files:            
                hz=str(file).split('.')[-1]
                if 'fa' in hz:
                    name_list=str(file).split('.')
                    name_list.remove(name_list[-1])
                    name_list.remove(name_list[-1])
                    name_join='.'.join(name_list)
                    assemblies_name[name_join]=0
        os.chdir(pwd)
    else:
        os.chdir(pwd)
        os.system('rm -rf '+target_bin_folder+'_OLC')
        os.system('cp -r '+target_bin_folder+'_backup '+target_bin_folder+'_OLC')

    os.system('rm *_db.txt')
    os.system('rm *.nsq')
    os.system('rm *.nin')
    os.system('rm *.nsi')
    os.system('rm *.nsd')
    os.system('rm *.nog')
    os.system('rm *.nhr')
    os.system('rm *.nhi')
    os.system('rm *.nhd')
    try:
        for i in range(1, max_num+1):
            os.system('rm Merged_'+str(i)+'*')
            os.system('rm Merged_seqs_'+str(i)+'*')
            os.system('rm *'+str(i)+'_merged_seq.txt')
            os.system('rm Filtrated_blast_'+str(i)+'_*')
            os.system('rm blast_'+str(i)+'_*')
            thresholds=[1, 1.5, 3]
            for item in thresholds:
                os.system('rm '+str(item)+'_'+str(i)+'_*')
    except:
        xyz=0
    
    try:
        for item in assemblies_name.keys():
            os.system('rm '+str(item)+'.*')
    except:
        xyz=0

    os.system('rm -rf *.fa_checkm')
    os.system('rm -rf *.fa_merged')
    os.system('rm -rf *.kmer.txt')
    os.system('rm -rf *.perf')
    os.system('rm -rf *_split_blast_output')
    os.system('rm *_split_blast_output.txt')
    os.system('rm -rf *_outlier')
    os.system('rm *_merged_seq.txt')
    os.system('rm -rf Merged_seqs_*')
    os.system('rm -rf *.ntf *.not *.nos *.njs *.ndb *.nto')

    print('Done!')

if __name__ == '__main__': 
    ### 1st OLC folder
    # target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'
    target_bin_folder='BestBinset_outlier_refined_filtrated_retrieved'
    bin_comparison_folder='BestBinset_comparison_files'
    step='assemblies_OLC'

    ### re-assembly folder
    #target_bin_folder='BestBinset_outlier_refined_t14_retrieved_total'
    #step='OLC_after_reassembly' ### 'assemblies_OLC' or 'OLC_after_reassembly' ; 'assemblies_OLC' means 1st OLC; 'OLC_after_reassembly' means OLC after reassembly
    #bin_comparison_folder='' ### if 'assemblies_OLC', folder bin_comparison_folder is needed; else if 'assemblies_OLC', the folder is not needed, just let it blank

    ### refinement mode needs to provide this folder list; if basalt start from autobinner, you shall let the following list blank
    # mod_bin_folder=['1_cat_bins','2_S001_bins','3_S002_bins']
    mod_bin_folder=[]
    
    num_threads=30
    aligned_len_cutoff=200
    similarity_cutoff=99
    coverage_extension=95
    OLC_main(target_bin_folder, step, bin_comparison_folder, aligned_len_cutoff, similarity_cutoff, coverage_extension, num_threads, mod_bin_folder)
