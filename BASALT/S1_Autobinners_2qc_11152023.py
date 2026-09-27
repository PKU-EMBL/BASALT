#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S1: Autobinning and initial quality control for BASALT.

This module converts read formats, prepares paired-end tracking
information and runs multiple binners to generate initial binsets
and PE-based contig connection files.
"""

from lib2to3.fixes import fix_buffer
from Bio import SeqIO
import sys, os, shlex, time, gc
from collections import Counter
from multiprocessing import Pool
import shutil
from basalt_runtime import bundled_executable

def fq2fa_conversion(filename):
    """
    Convert a FASTQ file into a FASTA file using Biopython.

    Parameters
    ----------
    filename : str
        Path to the input FASTQ file ('.fq' or '.fastq').

    Returns
    -------
    str
        Basename of the generated FASTA file (without path),
        with extension '.fasta'.
    """
    start=time.time()
    wrerr = sys.stderr.write
    print('---')
    print("Coverting "+str(filename)+" file to FA file")

    if '.fq' in str(filename):
        filename1=filename.replace('.fq','')
    elif '.fastq' in str(filename):
        filename1=filename.replace('.fastq','')

    file1Array = filename1.split('/')
    disName = file1Array[-1]
    pwd=os.getcwd()
    filename1 = os.path.join(pwd, disName)

    count = SeqIO.convert(filename, "fastq", filename1+".fasta", "fasta")
    print("Converted %i records" % count)
    end=time.time()
    wrerr("OK, Converting Finished in %3.2f secs\n" % (end-start))
    return disName+".fasta"

def ModifyEnd_fa(filename, n):
    """
    Add PE tracking suffix to FASTA record identifiers for later mapping.

    Parameters
    ----------
    filename : str
        Path to the FASTA file whose record IDs will be modified.
    n : int
        Index of the read mate (1 or 2) used to encode PE information.

    Returns
    -------
    str
        Name of the newly generated FASTA file with modified identifiers.
    """
    print('---')
    print('Adding end for PE-tracking of '+str(filename))
    print('--- pwd: ', os.getcwd())

    fout=open('PE_r'+str(n)+'_'+filename, 'w')
    m=0
    for record in SeqIO.parse(filename, 'fasta'):
        m+=1
        ids=str(m)+'_'+str(n)+'/'+str(n)
        fout.write('>'+str(ids)+'\n')
        fout.write(str(record.seq)+'\n')
    fout.close()
    print('Accomplished of adding end to '+str(filename))
    print('-----------------------------')
    return 'PE_r'+str(n)+'_'+filename

def ModifyEnd(filename, n):
    """
    Add PE tracking suffix to FASTQ record identifiers for later mapping.

    Parameters
    ----------
    filename : str
        Path to the FASTQ file whose record IDs will be modified.
    n : int
        Index of the read mate (1 or 2) used to encode PE information.

    Returns
    -------
    str
        Name of the newly generated FASTQ file with modified identifiers.
    """
    print('---')
    print('Adding end for PE-tracking of '+str(filename))
    print('--- pwd: ', os.getcwd())

    fout=open('PE_r'+str(n)+'_'+filename, 'w')
    m, m2=0, 0 
    for line in open(filename, 'r'):
        m+=1
        if '@' in line and m%4==1:
            m2+=1
            if ' ' in line:
                ids='@Seq'+str(m2)+'_modifiedID_'+str(n)+'/'+str(n)
            else:
                ids='@Seq'+str(m2)+'_modifiedID_'+str(n)+'/'+str(n)
            fout.write(str(ids)+'\n')
        else:
            fout.write(str(line))
    fout.close()
    print('Accomplished of adding end to '+str(filename))
    print('-----------------------------')
    return 'PE_r'+str(n)+'_'+filename

def PE_tracker(sam_file, output_name):
    """
    Parse a SAM file to infer paired-end contig connections.

    Parameters
    ----------
    sam_file : str
        Path to the SAM file containing PE mappings.
    output_name : str
        Output filename for the aggregated PE connections.

    Returns
    -------
    None
        Results are written to ``output_name``.
    """
    contig_pe, contig_pe_mock, n={}, {}, 0
    for line in open(sam_file,'r'):
        n+=1
        if '\t' in line:
            sam_list=str(line).strip().split('\t')
            if len(sam_file) >= 7:
                contig=str(line).strip().split('\t')[2]
                if '_1/1' in line:
                    reads=str(line).strip().split('\t')[0].split('_1/1')[0]
                elif '_2/2' in line:
                    reads=str(line).strip().split('\t')[0].split('_2/2')[0]
                else:
                    reads=str(line).strip().split('\t')[0]
                num=len(contig_pe_mock)
                contig_pe_mock[reads]={}
                if len(contig_pe_mock) > num:
                    contig_pe[reads]={}
                    contig_pe[reads][contig]=1
                else:
                    contig_pe[reads][contig]=1
        
        if n % 1000000 == 0 :
            print('Parsed '+str(n)+' lines')
    
    print('Parsed contigs connections')
    contig_connection_list, connections, connections_mock=[], {}, 'test'
    for item in contig_pe.keys():
        if len(contig_pe[item]) == 2:
            contig_connection_list.append(str(contig_pe[item]))
    contig_connection_list.sort()

    for item in contig_connection_list:
        if item != connections_mock:
            connections_mock=item
            connecting_contigs=str(item).replace('{','').replace('}','').replace('\'','').replace(': 1','')
            connections[connecting_contigs]=1
        else:
            connections[connecting_contigs]+=1
    
    f=open(output_name,'w')
    f.write('node1'+'\t'+'inter'+'\t'+'node2'+'\t'+'connections'+'\n')
    for item in connections:
        f.write(str(item).split(',')[0]+'\t'+'0'+'\t'+str(item).split(',')[1]+'\t'+str(connections[item])+'\n')
    f.close()

def cal_connections(connections):
    """
    Aggregate PE connection counts from multiple connection files.

    Parameters
    ----------
    connections : list of str
        List of filenames, each containing PE connection records
        produced by ``PE_tracker``.

    Returns
    -------
    dict
        Mapping 'id1\\tinter\\tid2' -> aggregated connection count.
    """
    PEC={}
    for item in connections:
        f=open(item, 'r')
        n=0
        for line in f:
            n+=1
            if n >= 2:
                id1=str(line).strip().split('\t')[0]
                id2=str(line).strip().split('\t')[2]
                inter=str(line).strip().split('\t')[1]
                cnt=int(str(line).strip().split('\t')[3])
                PE_ids=str(id1)+'\t'+str(inter)+'\t'+str(id2)
                if str(PE_ids) not in PEC.keys():
                    PEC[str(PE_ids)]=cnt
                else:
                    PEC[str(PE_ids)]+=cnt
    return PEC

def parse_lr_sam_hifi_connecting_contigs(sam_file):
    """
    Parse HiFi long-read SAM to find contigs jointly covered by the same read.

    Parameters
    ----------
    sam_file : str
        Path to the SAM file generated from HiFi long-read mapping.

    Returns
    -------
    None
        Appends results to ``Long_reads_connecting_contigs.txt``.
    """
    print('Processing '+sam_file+' for finding connecting contigs')
    lr_contig, lr_contig2 = {}, {}
    m, m1 = 0, 0
    for line in open(sam_file,'r'):
        m+=1
        flist=str(line).strip().split('\t')
        if len(flist) >= 12:
            m1 += 1
            read_id=flist[0].split('-')[0]
            contig_id=flist[2]
            # if len(read_id) != 0 and len(contig_id) != 0:
            try:
                lr_contig2[read_id]+=1
                try:
                    lr_contig[read_id][contig_id]+=1
                except:
                    lr_contig[read_id][contig_id]=1
            except:
                lr_contig2[read_id]=1
                lr_contig[read_id]={}
                lr_contig[read_id][contig_id]=1
            
        if m % 1000000 == 0:
            print('Read', m,'lines')

    # print(str(lr_contig))
    try:
        f=open('Long_reads_connecting_contigs.txt','a')
    except:
        f=open('Long_reads_connecting_contigs.txt','w')

    for read_id in lr_contig2.keys():
        if len(lr_contig[read_id]) >= 2:
            l={}
            # l, s={}, []
            for contig_id in lr_contig[read_id].keys():
                if lr_contig[read_id][contig_id] >= 14: ## 14 * 75 = 1050 bp
                    # s.append(lr_contig[read_id][contig_id])
                    l[contig_id]=lr_contig[read_id][contig_id]
            if len(l) >= 2:
                f.write(str(read_id)+'\t'+str(l)+'\n')
    f.close()

def parse_lr_sam_connecting_contigs(sam_file):
    """
    Parse generic long-read SAM to find contigs jointly covered by the same read.

    Parameters
    ----------
    sam_file : str
        Path to the SAM file generated from long-read mapping.

    Returns
    -------
    None
        Appends results to ``Long_reads_connecting_contigs.txt``.
    """
    print('Processing '+sam_file+' for finding connecting contigs')
    lr_contig, lr_contig2 = {}, {}
    m, m1 = 0, 0
    for line in open(sam_file,'r'):
        m+=1
        flist=str(line).strip().split('\t')
        if len(flist) >= 12:
            m1 += 1
            read_id=flist[0]
            contig_id=flist[2]
            try:
                lr_contig2[read_id]+=1
                try:
                    lr_contig[read_id][contig_id]+=1
                except:
                    lr_contig[read_id][contig_id]=1
            except:
                lr_contig2[read_id]=1
                lr_contig[read_id]={}
                lr_contig[read_id][contig_id]=1
            
        if m % 1000000 == 0:
            print('Read', m,'lines')

    # print(str(lr_contig))
    try:
        f=open('Long_reads_connecting_contigs.txt','a')
    except:
        f=open('Long_reads_connecting_contigs.txt','w')

    for read_id in lr_contig2.keys():
        if len(lr_contig[read_id]) >= 2:
            l={}
            # l, s={}, []
            for contig_id in lr_contig[read_id].keys():
                if lr_contig[read_id][contig_id] >= 14: ## 14 * 75 = 1050 bp
                    # s.append(lr_contig[read_id][contig_id])
                    l[contig_id]=lr_contig[read_id][contig_id]
            if len(l) >= 2:
                f.write(str(read_id)+'\t'+str(l)+'\n')
    f.close()

def mapping_lr_o(assembly, group, datasets, num_threads, pwd, data_type):
    """
    Map long-read datasets (ONT/PB) to the assembly and generate BAM files.

    Parameters
    ----------
    assembly : str
        Assembly filename (without group prefix).
    group : str
        Group identifier used to prefix output files.
    datasets : list of str
        List of long-read FASTQ/FASTA files.
    num_threads : int
        Number of threads to use for minimap2 and samtools.
    pwd : str
        Working directory.
    data_type : {'ont', 'pb'}
        Long-read platform, used to select minimap2 preset.

    Returns
    -------
    str
        Space-separated string of sorted BAM file names.
    """
    print('Mapping '+str(datasets)+' to contigs/scaffolds')
    try:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'a')
    except:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'w')

    for i in range(1, len(datasets)+1):
        if data_type == 'ont': 
            os.system('minimap2 -t '+str(num_threads)+' -ax map-ont '+str(group)+'_'+assembly+' '+str(datasets[i-1])+' > '+str(group)+'_lr'+str(i)+'.sam')
        elif data_type == 'pb': 
            os.system('minimap2 -t '+str(num_threads)+' -ax map-pb '+str(group)+'_'+assembly+' '+str(datasets[i-1])+' > '+str(group)+'_lr'+str(i)+'.sam')

        ### Parsing connections
        # parse_lr_sam_connecting_contigs(str(group)+'_lr'+str(i)+'.sam')

        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(group)+'_lr'+str(i)+'.sam -o '+str(group)+'_lr'+str(i)+'.bam')
        os.system('rm '+str(group)+'_lr'+str(i)+'.sam')

        print('Sorting bam file')
        ### py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(group)+'_lr'+str(i)+'.bam '+str(group)+'_lr'+str(i)+'_sorted') 

        try:
            with open(str(group)+'_lr'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(group)+'_lr'+str(i)+'.bam failed. Redoing')
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_lr'+str(i)+'_sorted.bam '+str(group)+'_lr'+str(i)+'.bam' )
    
        f_coverage_matrix.write('Coverage_list_lr'+str(i)+'.txt'+'\n')

        if i == 1:
            bam_sorted=str(group)+'_lr1_sorted.bam'
        else:
            bam_sorted+=' '+str(group)+'_lr'+str(i)+'_sorted.bam'
    f_coverage_matrix.close()
    print('Mapping Done!')
    print('-------------')
    return bam_sorted

def mapping_hifi_split(assembly, group, long_read_split_fa, num_threads, pwd):
    """
    Map split HiFi read pairs to the assembly using bowtie2.

    Parameters
    ----------
    assembly : str
        Assembly filename (without group prefix).
    group : str
        Group identifier used to prefix output files.
    long_read_split_fa : dict
        Mapping index -> [r1_fasta, r2_fasta] representing split HiFi pairs.
    num_threads : int
        Number of threads for bowtie2 and samtools.
    pwd : str
        Working directory.

    Returns
    -------
    str
        Space-separated string of sorted BAM file names.
    """
    print('Mapping datasets to contigs/scaffolds')
    # f_coverage_matrix=open('Coverage_list_'+assembly+'.txt', 'w')
    try:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'a')
    except:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'w')

    for i in range(1, len(long_read_split_fa)+1):
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -1 '+str(long_read_split_fa[i][0])+' -2 '+str(long_read_split_fa[i][1])+' -S '+str(group)+'_LR-'+str(i)+'-bw2.sam -f --no-unal')
        # os.system('bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -U '+str(long_read_split_fa[i-1])+' -S '+str(group)+'_LR-'+str(i)+'-bw2.sam -q --no-unal')
        # os.system('bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -U '+str(long_read_split_fa[i-1])+' -S '+str(group)+'_LR-'+str(i)+'-bw2.sam -f --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(group)+'_LR-'+str(i)+'-bw2.sam -o '+str(group)+'_LR-'+str(i)+'-bw2.bam')

        parse_lr_sam_hifi_connecting_contigs(str(group)+'_LR-'+str(i)+'-bw2.sam')
        os.system('rm '+str(group)+'_LR-'+str(i)+'-bw2.sam')
        ### py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(group)+'_LR-'+str(i)+'-bw2.bam '+str(group)+'_LR-'+str(i)+'-bw2_sorted') 

        try:
            with open(str(group)+'_LR-'+str(i)+'-bw2_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(group)+'_LR-'+str(i)+'-bw2.bam failed. Re doing')
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_LR-'+str(i)+'-bw2_sorted.bam '+str(group)+'_LR-'+str(i)+'-bw2.bam' )
        os.system('rm '+str(group)+'_LR-'+str(i)+'-bw2.bam')
       
        f_coverage_matrix.write('Coverage_list_LR-'+str(i)+'-bw2.txt'+'\n')

        if i == 1:
            bam_sorted=str(group)+'_LR-1-bw2_sorted.bam'
        else:
            bam_sorted+=' '+str(group)+'_LR-'+str(i)+'-bw2_sorted.bam'
    
    f_coverage_matrix.close()

    print('Mapping Done!')
    print('-------------')
    return bam_sorted

def mapping_hifi_minimap(assembly, group, long_read_split_fa, num_threads, pwd):
    """
    Map HiFi reads to the assembly using minimap2 (map-hifi preset).

    Parameters
    ----------
    assembly : str
        Assembly filename (without group prefix).
    group : str
        Group identifier used to prefix output files.
    long_read_split_fa : list of str
        List of HiFi read files to be mapped.
    num_threads : int
        Number of threads for minimap2 and samtools.
    pwd : str
        Working directory.

    Returns
    -------
    str
        Space-separated string of sorted BAM file names.
    """
    print('Mapping datasets to contigs/scaffolds')
    # f_coverage_matrix=open('Coverage_list_'+assembly+'.txt', 'w')
    try:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'a')
    except:
        f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'w')

    os.system('minimap2 -x map-hifi -t '+str(num_threads)+' -d '+str(group)+'_'+assembly+'.mmi '+str(group)+'_'+assembly)
    for i in range(1, len(long_read_split_fa)+1):
        print(str(long_read_split_fa[i-1])+' mapped against '+str(assembly))
        os.system('minimap2 -a --sam-hit-only --secondary=no -t '+str(num_threads)+' --split-prefix temp_PeopleGutHiFi1-5_cat '+str(group)+'_'+assembly+'.mmi '+str(long_read_split_fa[i-1])+' > '+str(group)+'_LR-'+str(i)+'.sam')

        # os.system('bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -U '+str(long_read_split_fa[i-1])+' -S '+str(group)+'_LR-'+str(i)+'.sam -f --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(group)+'_LR-'+str(i)+'.sam -o '+str(group)+'_LR-'+str(i)+'.bam')

        # parse_lr_sam_connecting_contigs(str(group)+'_LR-'+str(i)+'.sam')
        os.system('rm '+str(group)+'_LR-'+str(i)+'.sam')
        ### py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(group)+'_LR-'+str(i)+'.bam '+str(group)+'_LR-'+str(i)+'_sorted') 

        try:
            with open(str(group)+'_LR-'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(group)+'_LR-'+str(i)+'.bam failed. Re doing')
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_LR-'+str(i)+'_sorted.bam '+str(group)+'_LR-'+str(i)+'.bam' )
        
        os.system('rm '+str(group)+'_LR-'+str(i)+'.bam')
       
        f_coverage_matrix.write('Coverage_list_LR-'+str(i)+'.txt'+'\n')

        if i == 1:
            bam_sorted=str(group)+'_LR-1_sorted.bam'
        else:
            bam_sorted+=' '+str(group)+'_LR-'+str(i)+'_sorted.bam'
    
    f_coverage_matrix.close()

    # print('Scorting SAM file(s)')
    # print('CMD: jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+ str(bam_sorted))
    # try:
    #     os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
    #     nxxyy=0
    #     for line in open(str(group)+'_assembly.depth.txt','r'):
    #         nxxyy+=1
    #         if nxxyy == 2:
    #             break
    # except:
    #     os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))

    # n=0
    # for line in open(str(group)+'_assembly.depth.txt', 'r'):
    #     n+=1
    #     if n == 1:
    #         num=len(str(line).strip().split('\t'))-3
    #         num2=int(num/2)

    # for i in range(1, num2+1):
    #     n=0
    #     for line in open(str(group)+'_assembly.depth.txt', 'r'):
    #         n+=1
    #         if n == 1:
    #             nn=str(line).strip().split('\t')[2*int(i)+1].split('_sorted')[0].split('_')[1]
    #             f=open('Coverage_list_'+str(nn)+'.txt', 'w')

    #         if n > 1:
    #             ids=str(line).strip().split('\t')[0]
    #             coverage=str(line).strip().split('\t')[2*int(i)+1]
    #             f.write(str(ids)+'\t'+str(coverage)+'\n')
    #         else:
    #             continue
    #     f.close()
    # print('Done with generation of depth file!')
    # print('-------------')
    
    print('Mapping Done!')
    print('-------------')
    return bam_sorted

def mapping(assembly, group, datasets, num_threads, pwd):
    """
    Map paired-end short reads to the assembly using bowtie2.

    Parameters
    ----------
    assembly : str
        Assembly filename (without group prefix).
    group : str
        Group identifier used to prefix output files.
    datasets : dict
        Mapping dataset_id -> [R1, R2] FASTQ filenames.
    num_threads : int
        Number of threads for bowtie2 and samtools.
    pwd : str
        Working directory.

    Returns
    -------
    tuple
        (coverage_list_file, sorted_bam_string)
        where ``sorted_bam_string`` is a space-separated list of BAM files.
    """
    n=0
    logfile=open('Mapping_log_'+str(group)+'_'+assembly+'.txt', 'w')

    # print('Building Bowtie2 index')
    # os.system('bowtie2-build '+str(group)+'_'+assembly+' '+str(group)+'_'+assembly)
    # print('Done!')
    # print('-------------')

    print('Mapping datasets to contigs/scaffolds')
    f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'w')
    connections=[]
    for i in range(1, len(datasets)+1):
        logfile.write(str('Command: bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -1 '+str(datasets[str(i)][0])+' -2 '+str(datasets[str(i)][1])+' -S '+str(group)+'_DNA-'+str(i)+'.sam -q --no-unal')+'\n')
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(group)+'_'+assembly+' -1 '+str(datasets[str(i)][0])+' -2 '+str(datasets[str(i)][1])+' -S '+str(group)+'_DNA-'+str(i)+'.sam -q --no-unal')
        logfile.write(str('Command: samtools view -b -S '+str(group)+'_DNA-'+str(i)+'.sam -o '+str(group)+'_DNA-'+str(i)+'.bam')+'\n')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(group)+'_DNA-'+str(i)+'.sam -o '+str(group)+'_DNA-'+str(i)+'.bam')
        os.system('perl '+shlex.quote(bundled_executable('Cytoscapeviz.pl', __file__))+' -i '+shlex.quote(str(group)+'_DNA-'+str(i)+'.sam')+' -f 2 -a 150 -e 500 -m 3000 -c')
        os.system('mv condensed.cytoscape.connections.tab condensed.cytoscape.connections_'+str(group)+'_DNA-'+str(i)+'.tab')
        connections.append('condensed.cytoscape.connections_'+str(group)+'_DNA-'+str(i)+'.tab')
        os.system('rm '+str(group)+'_DNA-'+str(i)+'.sam')
        ### py2
        logfile.write(str('Command: samtools sort -@ '+str(num_threads)+' '+str(group)+'_DNA-'+str(i)+'.bam '+str(group)+'_DNA-'+str(i)+'_sorted')+'\n')
        os.system('samtools sort -@ '+str(num_threads)+' '+str(group)+'_DNA-'+str(i)+'.bam '+str(group)+'_DNA-'+str(i)+'_sorted') 

        try:
            with open(str(group)+'_DNA-'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(group)+'_DNA-'+str(i)+'.bam failed. Re doing')
            ### py3
            logfile.write(str('Command: samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_DNA-'+str(i)+'_sorted.bam '+str(group)+'_DNA-'+str(i)+'.bam')+'\n')
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_DNA-'+str(i)+'_sorted.bam '+str(group)+'_DNA-'+str(i)+'.bam' )
       
        f_coverage_matrix.write('Coverage_list_DNA-'+str(i)+'.txt'+'\n')

        if i == 1:
            bam_sorted=str(group)+'_DNA-1_sorted.bam'
        else:
            bam_sorted+=' '+str(group)+'_DNA-'+str(i)+'_sorted.bam'
    
    f_coverage_matrix.close()
    print('Mapping Done!')
    print('-------------')
    return bam_sorted

def bin_filtration(folder_name, pwd):
    """
    Filter bins shorter than specified thresholds from a bin folder.

    Parameters
    ----------
    folder_name : str
        Name of the bin folder containing FASTA bin files.
    pwd : str
        Working directory.

    Returns
    -------
    None
        Operates in-place on files in ``folder_name``.
    """
    os.chdir(pwd+'/'+str(folder_name))
    for root,dirs,files in os.walk(pwd+'/'+str(folder_name)):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                try:
                    seq_len=0
                    for record in SeqIO.parse(file,'fasta'):
                        seq_len+=int(len(record.seq))
                except:
                    seq_len=0
                
                if seq_len >= 30000000:
                    os.system('rm '+str(file))
    os.chdir(pwd)

def metabat(assembly_file, pwd, depth_file, threshold, num_threads):
    """
    Run MetaBAT2 binning on an assembly with a given depth file.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    pwd : str
        Working directory.
    depth_file : str
        Depth file in MetaBAT2-compatible format.
    threshold : list of float
        List of minCV thresholds to use.
    num_threads : int
        Number of threads for MetaBAT2.

    Returns
    -------
    list of str
        List of output bin folder names produced by MetaBAT2.
    """
    metabat_genome=str(assembly_file)+'_'+str(threshold)+'_metabat_genomes'
    # print(assembly_file)
    os.chdir(pwd+'/'+str(metabat_genome))
    print('Starting metabat autobinner in '+str(threshold))
    print('metabat2 -i '+str(assembly_file)+' -a '+str(depth_file)+' -o '+str(metabat_genome)+' --unbinned --maxEdges '+str(threshold))
    print('----------')
    os.system('metabat2 -t '+str(num_threads)+' -i '+str(assembly_file)+' -a '+str(depth_file)+' -o '+str(metabat_genome)+' --unbinned --maxEdges '+str(threshold))
    ###os.system('metabat2 -t '+str(num_threads)+' -i '+str(assembly_file)+' -a '+str(depth_file)+' -o '+str(metabat_genome)+' --unbinned --maxEdges '+str(threshold))
    os.system('mv '+str(metabat_genome)+'.unbinned.fa '+str(metabat_genome)+'.unbinned.txt')
    os.system('rm '+str(assembly_file))
    
    os.chdir(pwd)
    f=open('Autobinner_checkpoint.txt','a')
    f.write('Binning: '+str(metabat_genome)+'\n') #EMA: end modification accomplished
    f.close()
    bin_filtration(metabat_genome, pwd)
    # print('checking metabat bins with checkM')
    # print('----------')
    # metabat_checkm=str(assembly_file)+'_'+str(threshold)+'_metabat_checkm'
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(metabat_genome)+' '+str(metabat_checkm))
    print(assembly_file+' at '+str(threshold)+' Metabat autobinning done!')
    print('----------')

    # os.chdir(pwd+'/'+str(metabat_checkm)+'/storage')
    # p_bin_num=0
    # for line in open('bin_stats_ext.tsv','r'):
    #     p_bin_num+=1

def maxbin2(assembly_file, pwd, depth_file, threshold, Coverage_list_file, num_threads):
    """
    Run MaxBin2 binning on an assembly with a given depth file.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    pwd : str
        Working directory.
    depth_file : str
        Depth file in MaxBin2-compatible format.
    threshold : list of float
        List of min_prob thresholds to use.
    Coverage_list_file : str
        Output file listing coverage matrix files.
    num_threads : int
        Number of threads for MaxBin2.

    Returns
    -------
    list of str
        List of output bin folder names produced by MaxBin2.
    """
    maxbin2_genome=str(assembly_file)+'_'+str(threshold)+'_maxbin2_genomes'
    # maxbin2_checkm=str(assembly_file)+'_'+str(threshold)+'_maxbin2_checkm'
    print('Starting maxbin2 autobinner in'+str(threshold))
    print('run_MaxBin.pl -abund_list '+Coverage_list_file+' -thread '+str(num_threads)+' -contig '+str(assembly_file)+' -out '+str(maxbin2_genome)+' -prob_threshold '+str(threshold))
    fb=open('Basalt_log.txt','a')
    fb.write('Starting maxbin2 autobinner in'+str(threshold)+'\n')
    fb.write('run_MaxBin.pl -abund_list '+Coverage_list_file+' -thread '+str(num_threads)+' -contig '+str(assembly_file)+' -out '+str(maxbin2_genome)+' -prob_threshold '+str(threshold)+'\n')
    fb.close()
    os.system('cp '+str(assembly_file)+' Coverage_list* '+str(pwd)+'/'+str(maxbin2_genome))
    os.chdir(pwd+ '/'+str(maxbin2_genome))
    os.system('run_MaxBin.pl -abund_list '+Coverage_list_file+' -thread '+str(num_threads)+' -contig '+str(assembly_file)+' -out '+str(maxbin2_genome)+' -prob_threshold '+str(threshold))

def concoct_mod_file(assembly_file, depth_file):
    """
    Prepare CONCOCT input files (coverage and assembly) from given depth file.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    depth_file : str
        Depth file used to derive coverage matrix for CONCOCT.

    Returns
    -------
    tuple
        (concoct_assembly_path, concoct_depth_path)
    """
    concoct_assembly=open('Concoct_'+assembly_file,'w')
    record_ids={}
    for record in SeqIO.parse(assembly_file, 'fasta'):
        record_ids['X'+str(record.id)]=str(record.id)
        concoct_assembly.write('>X'+str(record.id)+'\n'+str(record.seq)+'\n')
    concoct_assembly.close()

    concoct_depth_file=open('Concoct_'+depth_file,'w')
    concoct_depth={}
    n=0
    for line in open(depth_file, 'r'):
        n+=1
        if n == 1:
            lis=str(line).strip().split('\t')
            n1=0
            for i in range(0, len(lis)):
                if i == 0:
                    title_line='contig'
                elif i%2==1 and i != 1:
                    n1+=1
                    title_line+='\t'+'DNA'+str(n1)
                else:
                    continue
            concoct_depth_file.write(str(title_line)+'\n')
        else:
            ids=str(line).strip().split('\t')[0]
            lis=str(line).strip().split('\t')
            n1=0
            for i in range(0, len(lis)):
                if i == 0:
                    coverage_line='X'+str(lis[0])
                elif i%2==1 and i != 1:
                    n1=+1
                    coverage_line+='\t'+str(lis[i])
                else:
                    continue
                n+=1
            concoct_depth_file.write(str(coverage_line)+'\n')
    concoct_depth_file.close()
    return str('Concoct_'+assembly_file), str('Concoct_'+depth_file)
    
def concoct(assembly_file, pwd, depth_file, threshold, num_threads):
    """
    Run CONCOCT binning on an assembly.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    pwd : str
        Working directory.
    depth_file : str
        Depth file used to derive coverage matrix.
    threshold : list
        Placeholder list of thresholds (kept for API compatibility).
    num_threads : int
        Number of threads for CONCOCT.

    Returns
    -------
    list of str
        List of output bin folder names produced by CONCOCT.
    """
    org_assembly=str(assembly_file).split('Concoct_')[1]
    concoct_genome=str(org_assembly)+'_'+str(threshold)+'_concoct_genomes'
    concoct_checkm=str(org_assembly)+'_'+str(threshold)+'_concoct_checkm'
    print('Starting concoct autobinner in '+str(threshold))
    print('concoct -c '+str(threshold)+' --coverage_file '+str(depth_file)+' --composition_file '+str(assembly_file)+' -b '+str(concoct_genome)+' --threads '+str(num_threads))
    fb=open('Basalt_log.txt','a')
    fb.write('Starting concoct autobinner in '+str(threshold)+'\n')
    fb.write('concoct -c '+str(threshold)+' --coverage_file '+str(depth_file)+' --composition_file '+str(assembly_file)+' -b '+str(concoct_genome)+' --threads '+str(num_threads)+'\n')
    fb.close()
    os.system('concoct -c '+str(threshold)+' --coverage_file '+str(depth_file)+' --composition_file '+str(assembly_file)+' -b '+str(concoct_genome)+' --threads '+str(num_threads))
    print('----------')
    ids_seq={}
    for record in SeqIO.parse(assembly_file, 'fasta'):
        ids_seq[str(record.id).split('X')[1]]=str(record.seq)

    os.chdir(pwd+'/'+str(concoct_genome))
    ids, n={}, 0
    for line in open('clustering_gt1000.csv', 'r'):
        n+=1
        if n >= 2:
            contigs_ids=str(line).strip().split(',')[0].split('X')[1]
            bin_ids=str(line).strip().split(',')[1]
            if str(bin_ids) not in ids.keys():
                ids[str(bin_ids)]=[str(contigs_ids)]
            else:
                ids[str(bin_ids)].append(str(contigs_ids))
    
    for bins in ids.keys():
        bin_name=str(concoct_genome)+'.'+str(bins)+'.fasta'
        f=open(bin_name, 'w')
        seq_len=0
        for contigs in ids[bins]:
            f.write('>'+str(contigs)+'\n'+str(ids_seq[contigs])+'\n')
            seq_len+=len(ids_seq[contigs])
        f.close()

        if seq_len >= 30000000:
            os.system('rm '+str(bin_name))

    os.chdir(pwd)
    # print('checking concoct bins with checkM')
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fasta '+str(concoct_genome)+' '+str(concoct_checkm))
    f=open('Autobinner_checkpoint.txt','a')
    f.write('Binning: '+str(concoct_genome)+'\n') #EMA: end modification accomplished
    f.close()
    fb=open('Basalt_log.txt','a')
    fb.write('Binning: '+str(concoct_genome)+'\n')
    fb.close()
    print(org_assembly+' in '+str(threshold)+' Concoct autobinning done!')
    print('----------')
    try:
        del ids_seq
        gc.collect()
    except:
        xyo=0


def lorbin(assembly_file, pwd, bam_sorted, num_threads):
    """
    Run LorBin binning given a sorted BAM file list.

    Parameters
    ----------
    assembly_file : str
        Assembly FASTA file.
    pwd : str
        Working directory.
    bam_sorted : str
        Space-separated sorted BAM files used for coverage.
    num_threads : int
        Number of threads for LorBin.

    Returns
    -------
    str
        Name of the LorBin output bin folder.
    """
    """
    LorBin autobinner：
      - assembly_file: 在 autobinners 里传入的组装文件名，比如 '1_assembly.fasta'
      - bam_sorted   : 字符串，多个 BAM 用空格分隔，比如 '1_DNA-1_sorted.bam 1_DNA-2_sorted.bam'
      - pwd          : 工作目录

    这里只负责跑 LorBin + 规范化输出，不做 checkm/checkm2，
    和 metabat/maxbin2/concoct 的风格一致。
    """
    os.chdir(pwd)
    lorbin_genome = str(assembly_file) + '_100_LorBin_genomes'

    # 输出目录
    # os.system('mkdir -p '+str(lorbin_genome))

    # BAM 的绝对路径
    bam_list = str(bam_sorted).split()
    bam_abs_list = [os.path.join(pwd, b) for b in bam_list]
    bam_abs = ' '.join(bam_abs_list)

    # 多样本时加 --multi
    multi_flag = ''
    if len(bam_list) > 1:
        multi_flag = ' --multi'

    cmd = (
        'LorBin bin'
        ' -o ' + lorbin_genome +
        ' -fa ' + str(assembly_file) +
        ' -b ' + bam_abs +
        '--num_process'+ num_threads + 
        multi_flag
    )

    print('Starting LorBin autobinner')
    print(cmd)
    os.system(cmd)
    print('LorBin binning finished for ' + lorbin_genome)

    # 规范化 bin 文件名为： <assembly_file>_100_LorBin_genomes.N.fa
    bin_index = 0
    outdir = os.path.join(pwd, lorbin_genome)
    for root, dirs, files in os.walk(outdir):
        for file in files:
            if file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fasta'):
                bin_index += 1
                old_path = os.path.join(root, file)
                new_name = lorbin_genome + '.' + str(bin_index) + '.fa'
                new_path = os.path.join(outdir, new_name)
                if old_path != new_path:
                    os.system('mv ' + old_path + ' ' + new_path)

    # trash_dir = os.path.join(outdir, 'output_bins')
    # shutil.rmtree(trash_dir)  # 递归删除文件夹及其内容
    # 跟其他 autobinner 一样做长度过滤
    bin_filtration(lorbin_genome, pwd)

    return lorbin_genome


def checkm_mul(num_threads, binset, binset_checkm_folder, checkm_done_f):
    """
    Run CheckM lineage workflow for a given binset if not already done.

    Parameters
    ----------
    num_threads : int
        Number of threads for CheckM.
    binset : str
        Path to the binset folder.
    binset_checkm_folder : str
        Output folder for CheckM results.
    checkm_done_f : dict
        Dictionary recording binsets that have already been processed.

    Returns
    -------
    None
        Writes CheckM outputs to ``binset_checkm_folder`` and updates
        checkpoint/log files.
    """
    if binset_checkm_folder not in checkm_done_f.keys():
        print('CheckM processing folder: '+str(binset))
        if 'concoct' in str(binset) or 'maxbin' in str(binset):
            # os.system('python S1_Checkm.py -i '+str(binset)+' -o '+str(binset_checkm_folder)+' -f fasta -t '+str(num_threads))
            os.system('checkm lineage_wf -t '+str(num_threads)+' -x fasta '+str(binset)+' '+str(binset_checkm_folder))
            os.system('rm -rf '+str(binset_checkm_folder)+'/bins')
        elif 'metabat' in str(binset):
            os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset)+' '+str(binset_checkm_folder))
            # os.system('python S1_Checkm.py -i '+str(binset)+' -o '+str(binset_checkm_folder)+' -f fa -t '+str(num_threads))
            os.system('rm -rf '+str(binset_checkm_folder)+'/bins')
        print(str(binset_checkm_folder)+' done!')
        fb=open('Basalt_log.txt','a')
        fb.write(str(binset_checkm_folder)+' done!'+'\n')
        fb.close()
        f=open('Autobinner_checkpoint.txt','a')
        f.write('Checkm: '+str(binset_checkm_folder)+'\n') #EMA: end modification accomplished
        f.close()
    else:
        print('Checkm of '+str(binset_checkm_folder)+' already been done in last run.')

def autobinners(softwares, assembly_file, depth_file, depth_file_list, Coverage_list_file, sensitive, binning_ds, checkm_done_f, QC, num_threads, ram, bam_sorted):
    """
    Run one or more autobinners (MetaBAT2, MaxBin2, CONCOCT, LorBin).

    Parameters
    ----------
    softwares : str
        Either 'all' or a specific binner selection (currently 'all' used).
    assembly_file : str
        Assembly FASTA file.
    depth_file : str
        Depth file for the assembly.
    depth_file_list : list
        List of depth files (kept for compatibility).
    Coverage_list_file : str
        Path to coverage list file to be updated.
    sensitive : {'quick', 'sensitive', 'more-sensitive'}
        Sensitivity level controlling threshold sets.
    binning_ds : dict
        Dictionary tracking which binsets have already been generated.
    checkm_done_f : dict
        Dictionary tracking which binsets have already run CheckM.
    QC : {'checkm', 'checkm2'}
        Quality control backend selection.
    num_threads : int
        Number of threads to use.
    ram : int
        Available RAM in gigabytes.
    bam_sorted : str
        Space-separated list of sorted BAM files for LorBin (optional).

    Returns
    -------
    dict
        Mapping from binset folder name to corresponding CheckM output folder.
    """
    pwd=os.getcwd()
    os.chdir(pwd)

    genome_folders=[]
    if str(softwares) == 'all':
        if sensitive == 'quick':
            # maxbin2_threshold=[0.3]
            metabat_threshold=[200, 300, 400, 500]
        elif sensitive == 'sensitive':
            # maxbin2_threshold=[0.3, 0.9]
            metabat_threshold=[200, 300, 400, 500]
        elif sensitive == 'more-sensitive':      
            maxbin2_threshold=[0.3, 0.5, 0.7, 0.9]
            # maxbin2_threshold=[0.3, 0.9]
            metabat_threshold=[200, 300, 400, 500]
            
        concoct_threshold=[]
        concoct_threshold_dict={}
        # concoct_threshold=[400, 200]
        c_m=concoct_mod_file(assembly_file, depth_file)
        concoct_assembly=c_m[0]
        concoct_depth=c_m[1]
        # maxbin2_threshold=[0.9]
        # metabat_threshold=[200]
        binset_checkm={}

        # lorbin_genome = str(assembly_file) + '_100_LorBin_genomes'
        # lorbin_checkm = str(assembly_file) + '_100_LorBin_checkm'
        # genome_folders.append(lorbin_genome)

        # 只在没跑过且有 BAM 的情况下跑
        # if lorbin_genome not in binning_ds.keys():
        #     if bam_sorted is None or str(bam_sorted).strip() == '':
        #         print('Warning: LorBin is enabled but bam_sorted is empty. Skip LorBin.')
        #     else:
        #         print('Running LorBin as the first autobinner')
        #         lorbin(assembly_file, pwd, bam_sorted,num_threads)

        #         # 记录到 checkpoint / log，和其他 binner 一致
        #         f = open('Autobinner_checkpoint.txt','a')
        #         f.write('Binning: ' + str(lorbin_genome) + '\n')
        #         f.close()
        #         fb = open('Basalt_log.txt','a')
        #         fb.write('Binning: ' + str(lorbin_genome) + '\n')
        #         fb.close()

        #         binning_ds[lorbin_genome] = ''
        # else:
        #     print('LorBin binning ' + str(lorbin_genome) + ' already done in last run.')

        # # 放进 binset_checkm，等后面统一跑 checkm/checkm2
        # binset_checkm[lorbin_genome] = lorbin_checkm

        ### maxbin2
        if sensitive == 'more-sensitive':
            mbn_num=len(maxbin2_threshold)
            xyz=0
            for line in open(depth_file, 'r'):
                xyz+=1

            max_num_project=1
            depth_num=[4000000, 3000000, 2000000, 1500000, 1000000, 500000, 250000]
            ram_each_list=[100, 75, 50, 35, 25, 10, 5]
            if xyz < 250000:
                max_num_project=mbn_num
            elif xyz >= 4000000:
                max_num_project=int(ram/100)
            else:
                for i in range(1,len(depth_num)):
                    if xyz >= depth_num[i] and xyz < depth_num[i-1]:
                        ram_each=ram_each_list[i]
                        max_num_project=int(ram/ram_each)

            if max_num_project == 0:
                max_num_project = 1
            
            if mbn_num <= max_num_project:
                num_threads_per_project=int(num_threads/mbn_num)
                project_num=mbn_num
            else:
                num_threads_per_project=int(num_threads/max_num_project)
                project_num=max_num_project

            pool=Pool(processes=project_num)
            for item in maxbin2_threshold:
                maxbin2_folder_name=str(assembly_file)+'_'+str(item)+'_maxbin2_genomes'
                genome_folders.append(maxbin2_folder_name)
                if maxbin2_folder_name not in binning_ds.keys():
                    try:
                        os.system('mkdir '+str(maxbin2_folder_name))
                        os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(maxbin2_folder_name))
                    except:
                        print(str(maxbin2_folder_name)+' already exist')
                        fb=open('Basalt_log.txt','a')
                        fb.write(str(maxbin2_folder_name)+' already exist'+'\n') #EMA: end modification accomplished
                        fb.close()
                    pool.apply_async(maxbin2, args=(assembly_file, pwd, depth_file, item, Coverage_list_file, num_threads_per_project))
            pool.close()
            pool.join()

            for item in maxbin2_threshold:
                maxbin2_genome=str(assembly_file)+'_'+str(item)+'_maxbin2_genomes'
                maxbin2_checkm=str(assembly_file)+'_'+str(item)+'_maxbin2_checkm'
                binset_checkm[str(assembly_file)+'_'+str(item)+'_maxbin2_genomes']=maxbin2_checkm
                os.chdir(pwd+'/'+str(maxbin2_genome))
                os.system('rm '+str(assembly_file))
                os.chdir(pwd)
                f=open('Autobinner_checkpoint.txt','a')
                f.write('Binning: '+str(maxbin2_genome)+'\n') #EMA: end modification accomplished
                f.close()
                fb=open('Basalt_log.txt','a')
                fb.write('Binning: '+str(maxbin2_genome)+'\n') #EMA: end modification accomplished
                fb.close()
                bin_filtration(maxbin2_genome, pwd)
            
            # print('checking maxbin2 bins with checkM')
            # print('----------')
            # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fasta '+str(maxbin2_genome)+' '+str(maxbin2_checkm))
            # print(assembly_file+' in '+str(item)+' Maxbin2 checkM done!')

        # for item in maxbin2_threshold:
        #     maxbin2_folder_name=str(assembly_file)+'_'+str(item)+'_maxbin2_genomes'
        #     genome_folders.append(maxbin2_folder_name)
        #     os.system('mkdir '+str(maxbin2_folder_name))
        #     os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(maxbin2_folder_name))
        #     maxbin2(assembly_file, pwd, depth_file, item, Coverage_list_file, num_threads)

        for item in metabat_threshold:
            metabat_folder_name=str(assembly_file)+'_'+str(item)+'_metabat_genomes'
            genome_folders.append(metabat_folder_name)
            if metabat_folder_name not in binning_ds.keys():
                # try:
                os.system('mkdir '+str(metabat_folder_name))
                # os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(metabat_folder_name))
                os.system('cp '+str(assembly_file)+' '+str(depth_file)+' '+pwd+'/'+str(metabat_folder_name))
                metabat(assembly_file, pwd, depth_file, item, num_threads)
                binset_checkm[metabat_folder_name]=str(assembly_file)+'_'+str(item)+'_metabat_checkm'
                # except:
                #     print(str(metabat_folder_name)+' already exist')
                #     fb=open('Basalt_log.txt','a')
                #     fb.write(str(metabat_folder_name)+' already exist'+'\n') #EMA: end modification accomplished
                #     fb.close()
                if len(depth_file_list) == 3:
                    for i in range(1,3):
                        threshold=str(i)+str(item)
                        metabat_folder_name2=str(assembly_file)+'_'+str(threshold)+'_metabat_genomes'
                        os.system('mkdir '+str(metabat_folder_name2))
                        genome_folders.append(metabat_folder_name2)
                        os.system('cp '+str(assembly_file)+' '+str(depth_file_list[i])+' '+pwd+'/'+str(metabat_folder_name2))
                        metabat(assembly_file, pwd, str(depth_file_list[i]), threshold, num_threads)
                        binset_checkm[metabat_folder_name2]=str(assembly_file)+'_'+str(threshold)+'_metabat_checkm'
                    
            p_bin_num=0
            os.chdir(pwd+'/'+str(metabat_folder_name))
            for root, dirs, files in os.walk(pwd+'/'+str(metabat_folder_name)):
                for file in files:
                    hz=file.split('.')[-1]
                    if 'fa' in hz:
                        p_bin_num+=1
            os.chdir(pwd)

            bin_num=100*(int(p_bin_num/100)+1)
            concoct_threshold_dict[bin_num]=''

        for bin_num in concoct_threshold_dict.keys():
            if bin_num not in concoct_threshold:
                concoct_threshold.append(bin_num)

            bin_num2=bin_num+100
            if bin_num2 not in concoct_threshold:
                concoct_threshold.append(bin_num2)
            
            # bin_num2=bin_num-100
            # if bin_num2 > 0 and bin_num2 not in concoct_threshold:
            #     concoct_threshold.append(bin_num2)

        # dl=[100,200,300,400,500]
        # for item in dl:
        #     if item not in concoct_threshold:
        #         concoct_threshold.append(item)


        if sensitive != 'quick':
            cct_num=len(concoct_threshold)
            num_threads_per_project=int(num_threads/cct_num)
            pool=Pool(processes=len(concoct_threshold))        
            for item in concoct_threshold:
                concoct_folder_name=str(assembly_file)+'_'+str(item)+'_concoct_genomes'
                genome_folders.append(concoct_folder_name)
                if concoct_folder_name not in binning_ds.keys():
                    os.system('mkdir '+str(concoct_folder_name))
                    os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(concoct_folder_name))
                    # concoct(concoct_assembly, pwd, concoct_depth, item, num_threads)
                    pool.apply_async(concoct, args=(concoct_assembly, pwd, concoct_depth, item, num_threads_per_project))
            pool.close()
            pool.join()

            ids_seq={}
            for record in SeqIO.parse(concoct_assembly, 'fasta'):
                ids_seq[str(record.id).split('X')[1]]=str(record.seq)

            for threshold in concoct_threshold:
                org_assembly=str(assembly_file)
                concoct_genome=str(org_assembly)+'_'+str(threshold)+'_concoct_genomes'
                concoct_checkm=str(org_assembly)+'_'+str(threshold)+'_concoct_checkm'
                binset_checkm[concoct_genome]=concoct_checkm

                os.chdir(pwd+'/'+str(concoct_genome))
                ids, n={}, 0
                for line in open('clustering_gt1000.csv', 'r'):
                    n+=1
                    if n >= 2:
                        contigs_ids=str(line).strip().split(',')[0].split('X')[1]
                        bin_ids=str(line).strip().split(',')[1]
                        if str(bin_ids) not in ids.keys():
                            ids[str(bin_ids)]=[str(contigs_ids)]
                        else:
                            ids[str(bin_ids)].append(str(contigs_ids))
                
                for bins in ids.keys():
                    bin_name=str(concoct_genome)+'.'+str(bins)+'.fasta'
                    f=open(bin_name, 'w')
                    seq_len=0
                    for contigs in ids[bins]:
                        f.write('>'+str(contigs)+'\n'+str(ids_seq[contigs])+'\n')
                        seq_len+=len(ids_seq[contigs])
                    f.close()

                    if seq_len >= 30000000:
                        os.system('rm '+str(bin_name))

                os.chdir(pwd)
                # print('checking concoct bins with checkM')
                # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fasta '+str(concoct_genome)+' '+str(concoct_checkm))
                print(org_assembly+' in '+str(threshold)+' Concoct autobinning done!')
                print('----------')
                fb=open('Basalt_log.txt','a')
                fb.write(org_assembly+' in '+str(threshold)+' Concoct autobinning done!'+'\n')
                fb.close()
            try:
                del ids_seq
                gc.collect()
            except:
                xyo=0


        print('Running checkm with binsets gernerated from autobinner')
        fb=open('Basalt_log.txt','a')
        fb.write('Running checkm with binsets gernerated from autobinner'+'\n')
        fb.close()

        ### checkm2
        for binset in binset_checkm.keys():
            # database_path='--database_path ~/databases/CheckM2_database/uniref100.KO.1.dmnd'
            if 'concoct' in str(binset) or 'maxbin' in str(binset):
                if QC == 'checkm2':
                    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(binset)+' -x fasta -o '+str(binset_checkm[binset]))
                elif QC == 'checkm':
                    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fasta '+str(binset)+' '+str(binset_checkm[binset]))

                # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(metabat_genome)+' '+str(metabat_checkm))
            elif 'metabat' in str(binset):                
                if QC == 'checkm2':
                    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(binset)+' -x fa -o '+str(binset_checkm[binset]))
                elif QC == 'checkm':
                    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset)+' '+str(binset_checkm[binset]))

        ### checkm
        # checkm_proj_num=int(ram/40)
        # checkm_folder_num=len(binset_checkm)
        # if checkm_proj_num > checkm_folder_num:
        #     checkm_proj_num=checkm_folder_num
            
        # num_threads_per_project=int(num_threads/checkm_proj_num)+1

        # pool=Pool(processes=checkm_proj_num)
        # for binset in binset_checkm.keys():
        #     pool.apply_async(checkm_mul, args=(num_threads_per_project, binset, binset_checkm[binset], checkm_done_f))
        # pool.close()
        # pool.join()

    elif str(softwares) == 'metabat':
        metabat_threshold=[200, 300, 400, 500]
        for item in metabat_threshold:
            metabat_folder_name=str(assembly_file)+'_'+str(item)+'_metabat_genomes'
            genome_folders.append(metabat_folder_name)
            os.system('mkdir '+str(metabat_folder_name))
            # os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(metabat_folder_name))
            os.system('cp '+str(assembly_file)+' '+str(depth_file)+' '+pwd+'/'+str(metabat_folder_name))
            metabat(assembly_file, pwd, depth_file, item, num_threads)

    elif str(softwares) == 'maxbin2':
        maxbin2_threshold=[0.3, 0.5, 0.7, 0.9]
        for item in maxbin2_threshold:
            maxbin2_folder_name=str(assembly_file)+'_'+str(item)+'_maxbin2_genomes'
            genome_folders.append(maxbin2_folder_name)
            os.system('mkdir '+str(maxbin2_folder_name))
            os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(maxbin2_folder_name))
            maxbin2(assembly_file, pwd, depth_file, item, Coverage_list_file, num_threads)

    elif str(softwares) == 'concoct':
        concoct_threshold=[400, 200]
        c_m=concoct_mod_file(assembly_file, depth_file)
        concoct_assembly=c_m[0]
        concoct_depth=c_m[1]
        for item in concoct_threshold:
            concoct_folder_name=str(assembly_file)+'_'+str(item)+'_concoct_genomes'
            genome_folders.append(concoct_folder_name)
            os.system('mkdir '+str(concoct_folder_name))
            os.system('cp '+str(depth_file)+' '+str(pwd)+'/'+str(concoct_folder_name))
            concoct(concoct_assembly, pwd, concoct_depth, item, num_threads)

    else:
        print('Error! Make sure you wrote the right name of binning software')
        fb=open('Basalt_log.txt','a')
        fb.write('Error! Make sure you wrote the right name of binning software'+'\n')
        fb.close()
    return genome_folders

def split_reads(reads, read_rename1, read_rename2, insert_size):
    """
    Split interleaved paired-end reads into two separate files.

    Parameters
    ----------
    reads : str
        Input FASTQ file with interleaved pairs.
    read_rename1 : str
        Output filename for mate 1 reads.
    read_rename2 : str
        Output filename for mate 2 reads.
    insert_size : int
        Estimated insert size (kept for compatibility, not used here).

    Returns
    -------
    None
        Writes two FASTQ files with split read pairs.
    """
    print('Splitting '+str(reads))
    f=open(read_rename1,'w')
    f2=open(read_rename2,'w')
    # even={0:0, 2:2, 4:4, 6:6, 8:8}
    n, m = 0, 0
    for line in open(reads,'r'):
        n+=1
        if n%4 == 1:
            m+=1
            # ids=str(line).strip().replace(' ','_')
            # ids=str(line).strip().replace(' ','_').replace('@','>')
            # ids='@Seq'+str(m)
            ids='>Seq'+str(m)
        elif n%4 == 2:
            lines=line.strip()
            ll=len(lines)
            if ll%insert_size != 0:
                sn=int(ll/insert_size)+1
            else:
                sn=int(ll/insert_size)

            if sn % 2 == 0:
                sn=sn
            else:
                sn=sn-1
            for i in range(1,sn+1):
                s, e=insert_size*(i-1), insert_size*i
                lines=line.strip()
                paired_num=int((i+1)/2)
                if i%2 == 0:
                    f2.write(str(ids)+'-pair'+str(paired_num)+'_modifiedID 2_frag'+str(i)+'\n'+str(str(lines)[s:e])+'\n')
                else:
                    f.write(str(ids)+'-pair'+str(paired_num)+'_modifiedID 1_frag'+str(i)+'\n'+str(str(lines)[s:e])+'\n')
        #     seq_l=[]
        #     if sn % 2 == 0:
        #         sn=sn
        #     else:
        #         sn=sn-1
        #     for i in range(1,sn+1):
        #         s, e=insert_size*(i-1), insert_size*i
        #         seq_l.append(str(lines)[s:e])
        #         # f.write(str(ids)+'_'+str(i)+'\n'+str(str(lines)[s:e])+'\n')
        # elif n%4 == 0:
        #     for i in range(1,sn+1):
        #         s, e=insert_size*(i-1), insert_size*i
        #         lines=line.strip()
        #         paired_num=int((i+1)/2)
        #         if i%2 == 0:
        #             f2.write(str(ids)+'-pair'+str(paired_num)+'_modifiedID 2_frag'+str(i)+'\n'+str(seq_l[i-1])+'\n'+'+'+'\n'+str(str(lines)[s:e])+'\n')
        #         else:
        #             f.write(str(ids)+'-pair'+str(paired_num)+'_modifiedID 1_frag'+str(i)+'\n'+str(seq_l[i-1])+'\n'+'+'+'\n'+str(str(lines)[s:e])+'\n')

    f.close()
    f2.close()
    print('Split '+str(reads)+' done!')

def autobinner_main(assembly_list, datasets, lr, hifi_list, insert_size, num_threads, ram, sensitive, QC, pwd):
    """
    Top-level entry point for S1 autobinning.

    Parameters
    ----------
    assembly_list : list of str
        List of assembly FASTA filenames.
    datasets : dict
        Mapping dataset_id -> [R1, R2] short-read FASTQ pairs.
    lr : list of str
        Long-read datasets (ONT/PB).
    hifi_list : list of str
        HiFi read datasets.
    insert_size : int
        Estimated insert size for PE reads.
    num_threads : int
        Number of CPU threads to use.
    ram : int
        Available RAM in gigabytes.
    sensitive : {'quick', 'sensitive', 'more-sensitive'}
        Autobinning sensitivity level.
    QC : {'checkm', 'checkm2'}
        Quality control backend to use.
    pwd : str
        Working directory.

    Returns
    -------
    tuple
        (bins_folders_dic, connections_total_dict, depth_total, assembly_MoDict)
        where each component summarises bin folders, PE connections,
        depth files and modified assembly names.
    """
    bins_folders, datasets_fq, connections_total_dict, depth_total, assembly_MoDict={}, {}, {}, {}, {}

    try:
        f=open('Autobinner_checkpoint.txt','a')
    except:
        f=open('Autobinner_checkpoint.txt','w')
    f.close()

    try:
        fb=open('Basalt_log.txt','a')
    except:
        fb=open('Basalt_log.txt','w')
    fb.close()
    
    end_mo_acc, ab_acc_asb, mapping_ds, binning_ds, checkm_done_f, binned_proj = {}, {}, {}, {}, {}, {}
    for line in open('Autobinner_checkpoint.txt', 'r'):
        if 'EMA: ' in line: #EMA: end modification accomplished
            processed_dataset=line.strip().split('EMA: ')[1]
            end_mo_acc[processed_dataset]=''
        elif 'ABA: ' in line: #ABA: auto-binning accomplished
            processed_dataset=line.strip().split('ABA: ')[1]
            ab_acc_asb[processed_dataset]=''
        elif 'Mapping: ' in line: #ABA: auto-binning accomplished
            processed_dataset=line.strip().split('Mapping: ')[1]
            mapping_ds[processed_dataset]=''
        elif 'Binning: ' in line: #ABA: auto-binning accomplished
            processed_dataset=line.strip().split('Binning: ')[1].strip()
            binning_ds[processed_dataset]=''
            binned_proj[processed_dataset]=''
            assembly_name_list=str(processed_dataset).split('_')
            assembly_name='_'.join(assembly_name_list[1:-3])
            try:
                bins_folders[str(assembly_name)].append(processed_dataset)
            except:
                bins_folders[str(assembly_name)]=[processed_dataset]
        elif 'Checkm: ' in line: #ABA: auto-binning accomplished
            processed_dataset=line.strip().split('Checkm: ')[1]
            checkm_done_f[processed_dataset]=''

    if len(datasets) != 0:
        for item in datasets.keys():
            if str(item) not in end_mo_acc.keys():
                datasets_fq[item]=[]
                datasets_fq[item].append(ModifyEnd(datasets[item][0], 1))
                # os.system('rm '+str(datasets[item][0]))
                datasets_fq[item].append(ModifyEnd(datasets[item][1], 2))
                # os.system('rm '+str(datasets[item][1]))
                print('End modification accomplished: '+str(item))
                f=open('Autobinner_checkpoint.txt','a')
                f.write('EMA: '+str(item)+'\n') #EMA: end modification accomplished
                f.close()
                fb=open('Basalt_log.txt','a')
                fb.write('EMA: '+str(item)+'\n') #EMA: end modification accomplished
                fb.close()
            else:
                datasets_fq[item]=[]
                datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
                datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))
            # datasets_fq[item].append(str(datasets[item][0]))
            # datasets_fq[item].append(str(datasets[item][1]))

    for item in range(0, len(assembly_list)):
        mo_assembly=str(int(item)+1)+'_'+str(assembly_list[item])
        mo_assembly_depth=str(int(item)+1)+'_assembly.depth.txt'
        if str(item) not in ab_acc_asb.keys():
            PEC={}
            group=int(item)+1
            assembly=str(assembly_list[item])
            f_filtrated=open(str(group)+'_'+assembly, 'w')
            print('Filtration of the contigs/scaffolds. The process keeps Contig/scaffold with larger than 1000 bp.')
            n=0
            for record in SeqIO.parse(assembly,'fasta'):
                if len(record.seq) >= 1000:
                    n+=1
                    f_filtrated.write('>'+str(group)+'-'+str(n)+'\n'+str(record.seq)+'\n')
            f_filtrated.close()
            print('Done!')
            print('-------------')

            depth_file_list=[]
            if str(group) not in mapping_ds.keys():
                tok=0
                if len(datasets_fq) != 0:
                    print('Building Bowtie2 index')
                    os.system('bowtie2-build '+str(group)+'_'+assembly+' '+str(group)+'_'+assembly)
                    print('Done!')
                    tok=1

                    print('-------------')
                    bam_sorted1=mapping(str(assembly_list[item]), int(item)+1, datasets_fq, num_threads, pwd)
                    # bam_sorted1='1_DNA-1_sorted.bam'

                connections=[]
                for i in range(1, len(datasets_fq)+1):
                    connections.append('condensed.cytoscape.connections_'+str(int(item)+1)+'_DNA-'+str(i)+'.tab')

                PEC=cal_connections(connections)

                connections_total=open('condense_connections_'+str(assembly_list[item])+'.txt', 'w')
                connections_total.write('node1'+'\t'+'interaction'+'\t'+'node2'+'\t'+'connections'+'\n')
                for item2 in PEC.keys():
                    connections_total.write(str(item2)+'\t'+str(PEC[item2])+'\n')
                connections_total.close()

                sort_bam={}
                if len(lr) != 0 or len(hifi_list) != 0:
                    print('Start mapping long-reads')
                    if len(hifi_list) != 0:
                        if tok == 0:
                            print('Building Bowtie2 index')
                            os.system('bowtie2-build '+str(group)+'_'+assembly+' '+str(group)+'_'+assembly)
                            print('Done!')
                            tok=1

                        long_read_split_fa={}
                        pool=Pool(processes=num_threads)
                        n=0
                        for reads in hifi_list:
                            n+=1
                            if '.fq' in reads:
                                read_rename1=reads.replace('.fq','_split_R1.fa')
                                read_rename2=reads.replace('.fq','_split_R2.fa')
                                # read_rename1=reads.replace('.fq','_split_R1.fq')
                                # read_rename2=reads.replace('.fq','_split_R2.fq')
                                # read_rename=reads.replace('.fq','_split.fa')
                            elif '.fastq' in reads:
                                read_rename1=reads.replace('.fastq','_split_R1.fa')
                                read_rename2=reads.replace('.fastq','_split_R2.fa')
                                # read_rename1=reads.replace('.fastq','_split_R1.fq')
                                # read_rename2=reads.replace('.fastq','_split_R2.fq')
                                # read_rename=reads.replace('.fastq','_split.fa')

                            long_read_split_fa[n]=[]
                            long_read_split_fa[n].append(read_rename1)
                            long_read_split_fa[n].append(read_rename2)
                            # read_rename=split_reads(reads)
                            pool.apply_async(split_reads,args=(reads,read_rename1,read_rename2, insert_size))
                        pool.close()
                        pool.join()
                        
                        if len(hifi_list) == 1:
                            bam_hifi_bw2_sorted=mapping_hifi_split(assembly, group, long_read_split_fa, num_threads, pwd)
                            bam_hifi_minimap2_sorted=mapping_hifi_minimap(assembly, group, hifi_list, num_threads, pwd)
                            bam_sorted2=bam_hifi_bw2_sorted+' '+bam_hifi_minimap2_sorted
                            sort_bam[bam_hifi_bw2_sorted]=0
                            sort_bam[bam_hifi_minimap2_sorted]=0
                            
                        else:
                            bam_hifi_bw2_sorted=mapping_hifi_split(assembly, group, long_read_split_fa, num_threads, pwd)
                            bam_sorted2=bam_hifi_bw2_sorted
                            sort_bam[bam_hifi_bw2_sorted]=0

                    if len(lr) != 0:
                        data_type='ont'
                        bam_sorted2=mapping_lr_o(assembly, group, lr, num_threads, pwd, data_type)

                if len(datasets_fq) == 0:
                    bam_sorted=bam_sorted2
                elif len(datasets_fq) != 0:
                    if len(lr) == 0 and len(hifi_list) == 0:
                        bam_sorted=bam_sorted1
                    else:
                        bam_sorted=bam_sorted1+' '+bam_sorted2
                
                # print(str(bam_sorted))

                depth_file_list=[]
                if len(sort_bam) != 2: ### Generate total depth; len(sort_bam) == 0: there is only HTS data; len(sort_bam) == 1: there is only bw2 mapping result; 
                    # logfile=open('Mapping_log_'+str(group)+'_'+assembly+'.txt', 'a')
                    print('Scorting SAM file(s)')
                    print('CMD: jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+ str(bam_sorted))
                    # logfile.write(str('Command: jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))+'\n')
                    try:
                        os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
                        nxxyy=0
                        for line in open(str(group)+'_assembly.depth.txt','r'):
                            nxxyy+=1
                            if nxxyy == 2:
                                break
                    except:
                        os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
                    # logfile.close()
                    depth_file_list.append(str(group)+'_assembly.depth.txt')
                else:
                    print('Scorting SAM file(s)')
                    print('CMD: jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+ str(bam_sorted))
                    try:
                        os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
                        nxxyy=0
                        for line in open(str(group)+'_assembly.depth.txt','r'):
                            nxxyy+=1
                            if nxxyy == 2:
                                break
                    except:
                        os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
                    
                    depth_file_list.append(str(group)+'_assembly.depth.txt')
                    depth_file_list.append(str(group)+'_assembly.depth_1.txt')
                    depth_file_list.append(str(group)+'_assembly.depth_2.txt')

                    fd1=open(str(group)+'_assembly.depth_1.txt','w')
                    fd2=open(str(group)+'_assembly.depth_2.txt','w')
                    for line in open(str(group)+'_assembly.depth.txt', 'r'):
                        content_l=str(line).strip().split('\t')
                        fd1.write(str(content_l[0])+'\t'+str(content_l[1])+'\t'+str(content_l[2])+'\t'+str(content_l[3])+'\t'+str(content_l[4])+'\n')
                        fd2.write(str(content_l[0])+'\t'+str(content_l[1])+'\t'+str(content_l[2])+'\t'+str(content_l[5])+'\t'+str(content_l[6])+'\n')
                    fd1.close()
                    fd2.close()
                        
                n=0
                for line in open(str(group)+'_assembly.depth.txt', 'r'):
                    n+=1
                    if n == 1:
                        num=len(str(line).strip().split('\t'))-3
                        num2=int(num/2)

                for i in range(1, num2+1):
                    n=0
                    for line in open(str(group)+'_assembly.depth.txt', 'r'):
                        n+=1
                        if n == 1:
                            nn=str(line).strip().split('\t')[2*int(i)+1].split('_sorted')[0].split('_')[1]
                            f=open('Coverage_list_'+str(nn)+'.txt', 'w')

                        if n > 1:
                            ids=str(line).strip().split('\t')[0]
                            coverage=str(line).strip().split('\t')[2*int(i)+1]
                            f.write(str(ids)+'\t'+str(coverage)+'\n')
                        else:
                            continue
                    f.close()
                print('Done with generation of depth file!')
                print('-------------')
                f=open('Autobinner_checkpoint.txt','a')
                f.write('Mapping: '+str(group)+'\n') #EMA: end modification accomplished
                f.close()
            else:
                if len(datasets_fq) != 0:
                    for i in range(1, len(datasets_fq)+1):
                        if i == 1:
                            bam_sorted1='1_DNA-1_sorted.bam'
                        else:
                            bam_sorted1+=' 1_DNA-'+str(i)+'_sorted.bam'
                    
                if len(lr) != 0:
                    for i in range(1, len(lr)+1):
                        if i == 1:
                            bam_sorted2='1_LR-1_sorted.bam 1_LR-1-bw2_sorted.bam'
                        else:
                            bam_sorted2+=' 1_LR-'+str(i)+'_sorted.bam 1_LR-'+str(i)+'-bw2_sorted.bam'

                if len(datasets_fq) != 0 and len(lr) == 0:
                    bam_sorted=bam_sorted1
                elif len(datasets_fq) == 0 and len(lr) != 0:
                    bam_sorted=bam_sorted2
                elif len(datasets_fq) != 0 and len(lr) != 0:
                    bam_sorted=bam_sorted1+' '+bam_sorted2
                    
            ### metabat etc.
            # Coverage_list_file='Coverage_list_'+str(int(item)+1)+'_'+str(assembly_list[item])+'.txt'
            # bins_folders[str(assembly_list[item])]=autobinners('all', mo_assembly, mo_assembly_depth, depth_file_list, Coverage_list_file, sensitive, binning_ds, checkm_done_f, QC, num_threads, ram)
            # os.system('rm Coverage_list_*')
            # bam_sorted = 
            ### metabat / maxbin2 / concoct 
            Coverage_list_file='Coverage_list_'+str(int(item)+1)+'_'+str(assembly_list[item])+'.txt'
            bins_folders[str(assembly_list[item])] = autobinners(
                'all',
                mo_assembly,
                mo_assembly_depth,
                depth_file_list,
                Coverage_list_file,
                sensitive,
                binning_ds,
                checkm_done_f,
                QC,
                num_threads,
                ram,
                bam_sorted1  
            )
            os.system('rm Coverage_list_*')

            # ### SemiBin2
            print('Performing Semibin2')
            # if len(depth_file_list) == 1:
            semibin2_list=[]
            # bam_sorted='x y'
            if len(lr) != 0 or len(hifi_list) != 0:
                semibin_folder_name=str(group)+'_'+assembly+'_100_semibin_genomes'
                if semibin_folder_name not in binning_ds.keys():
                    os.system('SemiBin2 single_easy_bin -i '+str(group)+'_'+assembly+' -b '+str(bam_sorted)+' -o '+str(group)+'_'+assembly+'_100_semibin_genomes --sequencing-type=long_read --processes '+str(num_threads))
                    semibin2_list.append(str(group)+'_'+assembly+'_100_semibin_genomes')

                    f=open('Autobinner_checkpoint.txt','a')
                    f.write('Binning: '+str(group)+'_'+assembly+'_100_semibin_genomes'+'\n')
                    f.close()

                # if len(hifi_list) == 1:
                #     bam_list=str(bam_sorted).split(' ')
                #     xyz=100
                #     semibin_folder_name=str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes'
                #     for sorted_bam in bam_list:
                #         xyz+=100
                #         if semibin_folder_name not in binning_ds.keys():
                #             os.system('SemiBin2 single_easy_bin -i '+str(group)+'_'+assembly+' -b '+str(sorted_bam)+' -o '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes --sequencing-type=long_read')
                #             semibin2_list.append(str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes')

                #             f=open('Autobinner_checkpoint.txt','a')
                #             f.write('Binning: '+str(semibin_folder_name)+'\n')
                #             f.close()
            else:
                semibin_folder_name=str(group)+'_'+assembly+'_100_semibin_genomes'
                if semibin_folder_name not in binning_ds.keys():
                    os.system('SemiBin2 single_easy_bin -i '+str(group)+'_'+assembly+' -b '+str(bam_sorted)+' -o '+str(group)+'_'+assembly+'_100_semibin_genomes --processes '+str(num_threads))
                    semibin2_list.append(str(group)+'_'+assembly+'_100_semibin_genomes')

                    f=open('Autobinner_checkpoint.txt','a')
                    f.write('Binning: '+str(group)+'_'+assembly+'_100_semibin_genomes'+'\n')
                    f.close()

            for item_xyz in semibin2_list:
                ### item = str(group)+'_'+assembly+'_100_semibin_genomes'
                xyz=item_xyz.split('_semibin')[0].split('_')[-1]

                os.system('mv '+pwd+'/'+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes/output_bins/* '+pwd+'/'+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes')
                os.chdir(pwd+'/'+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes')
                os.system('rm -rf output_bins')
                os.system('gzip -d *.gz')

                os.chdir(pwd+'/'+str(item_xyz))
                for root, dirs, files in os.walk(pwd+'/'+str(item_xyz)):
                    for file in files:
                        if '.fa' in file:
                            bin_num=str(file).split('SemiBin_')[1].split('.fa')[0]
                            os.system('mv '+str(file)+' '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes.'+str(bin_num)+'.fa')

                os.chdir(pwd)
                if QC == 'checkm2':
                    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes -x fa -o '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_checkm')
                elif QC == 'checkm':
                    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes '+str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_checkm')

                try:
                    bins_folders[str(assembly_list[item])].append(str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes')
                except:
                    bins_folders[str(assembly_list[item])]=[]
                    bins_folders[str(assembly_list[item])].append(str(group)+'_'+assembly+'_'+str(xyz)+'_semibin_genomes')

            # Collecting single contig bin
            
            os.system('mkdir '+str(group)+'_'+assembly+'_1_SingleContig_genomes')
            f2=open('Potential_bins.txt','w')
            p_b, n, y={}, 0, 0
            for record in SeqIO.parse(str(group)+'_'+assembly, 'fasta'):
                n+=1
                if len(record.seq) >= 500000:
                    y+=1
                    f2.write(str(record.id)+'\t'+str(len(record.seq))+'\n')
                    f=open(str(group)+'_'+assembly+'_1_SingleContig_genomes.'+str(n)+'.fa','w')
                    f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    f.close()
                    p_b[str(group)+'_'+assembly+'_1_SingleContig_genomes.'+str(n)+'.fa']=0
            f2.close()

            if y >= 1:
                for itemx in p_b.keys():
                    os.system('mv '+str(itemx)+' '+pwd+'/'+str(group)+'_'+assembly+'_1_SingleContig_genomes')
                
                if QC == 'checkm2':
                    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(group)+'_'+assembly+'_1_SingleContig_genomes -x fa -o '+str(group)+'_'+assembly+'_1_SingleContig_checkm')
                elif QC == 'checkm':
                    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(group)+'_'+assembly+'_1_SingleContig_genomes '+str(group)+'_'+assembly+'_1_SingleContig_checkm')

                bins_folders[str(assembly)].append(str(group)+'_'+assembly+'_1_SingleContig_genomes')
            else:
                os.system('rm -rf '+str(group)+'_'+assembly+'_1_SingleContig_genomes')

            f=open('Autobinner_checkpoint.txt','a')
            f.write('ABA: '+str(item)+'\n') #EMA: end modification accomplished
            f.close()
            fb=open('Basalt_log.txt','a')
            fb.write('ABA: '+str(item)+'\n') #EMA: end modification accomplished
            fb.close()
        connections_total_dict[str(assembly_list[item])]='condense_connections_'+str(assembly_list[item])+'.txt'
        depth_total[str(assembly_list[item])]=mo_assembly_depth
        assembly_MoDict[str(assembly_list[item])]=mo_assembly

    # os.system('rm *.bam')
    print('Auto-binning done!')
    return bins_folders, connections_total_dict, depth_total, assembly_MoDict, datasets_fq

if __name__ == '__main__': 
    # assembly_list=['HumanGut_HiFi1-5_cat_assembly.fasta']
    # assembly_list=['HumanGut_cat.fasta']
    assembly_list=['assembly.fasta']
    num_threads=60
    ram=250
    pwd=os.getcwd()
    datasets={}
    # datasets={'1':['SRR10988543_1.fastq','SRR10988543_2.fastq']}
    lr_list=[]
    hifi_list=['sample1.R1.fq', 'sample1.R2.fq'] ### put all the long reads into this list, including ont, pb dataset, hifi dataset etc.
    insert_size=100
    QC='checkm2' ### checkm2 or checkm
    # lr=['SRR15275210.fastq', 'SRR15275211.fastq', 'SRR15275212.fastq', 'SRR15275213.fastq', 'SRR17687125.fastq']
    # binning_mode='more-sensitive' ### 'quick', 'sensitive', 'more-sensitive'
    sensitive='sensitive' ### default: 'sensitive'; option: 'quick', 'more-sensitive'
    autobinner_main(assembly_list, datasets, lr_list, hifi_list, insert_size, num_threads, ram, sensitive, QC, pwd)
