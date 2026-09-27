#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Core routines for BASALT data feeding.

These functions merge external binsets with existing BASALT results,
recalculate coverage matrices and generate inputs for downstream
dereplication and refinement modules.
"""

from Bio import SeqIO
import sys
import os
import shlex
import time
from collections import Counter
from basalt_runtime import bundled_executable, coerce_count, depth_sample_count

def merge_bin(binset_folder, bs_id, pwd):
    bins_checkm={}
    try:
        os.mkdir(str(bs_id)+'_'+str(binset_folder))
    except:
        print(str(bs_id)+'_'+str(binset_folder)+' is exist. Re-create the folder')
        os.system('rm -rf '+str(bs_id)+'_'+str(binset_folder))
        os.mkdir(str(bs_id)+'_'+str(binset_folder))
    
    os.chdir(pwd+'/'+binset_folder)
    f=open('Mod_contig_id.txt','w')
    f1=open('Bin_name_mod.txt','w')
    # f2=open('Bin_name_mod_bin_stats_ext.tsv','w')
    n, record_seq, mod_bin_list = 0, {}, []
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            # print(hz
            if 'fa' in hz:
                n+=1
                m=0
                record_seq[str(bs_id)+'_bin'+str(n)]={}
                mod_bin_list.append(str(bs_id)+'_bin'+str(n))
                f1.write(str(file)+'\t'+str(bs_id)+'_bin'+str(n)+'.fa'+'\n')
                for record in SeqIO.parse(file, 'fasta'):
                    m+=1
                    f.write(str(bs_id)+'_bin'+str(n)+'_'+str(m)+'\t'+str(record.id)+'\n')
                    record_seq[str(bs_id)+'_bin'+str(n)][str(bs_id)+'_bin'+str(n)+'_'+str(m)]=str(record.seq)
    f.close()
    f1.close()

    os.system('mv Mod_contig_id.txt Bin_name_mod.txt '+pwd+'/'+str(binset_folder))
    
    os.chdir(pwd+'/'+str(bs_id)+'_'+str(binset_folder))
    # f1=open(str(bs_id)+'_'+binset_folder+'_seq.fa','w')
    for bin_id in record_seq.keys():
        f=open(str(bin_id)+'.fa','w')
        for contigs in record_seq[bin_id].keys():
            f.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
            # f1.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
        f.close()
    # f1.close()

    os.system('cat *.fa > '+str(pwd)+'/'+str(bs_id)+'_'+binset_folder+'_seq.fa')
    os.system('cat *.fa > '+str(pwd)+'/'+str(bs_id)+'_'+binset_folder+'_seq_test2.fa')
    os.chdir(pwd)
    return str(bs_id)+'_'+str(binset_folder), str(bs_id)+'_'+binset_folder+'_seq.fa', mod_bin_list, record_seq

def intervalue(Xmin, Xmax, Y, Z):
    delta=Xmax-Xmin
    tim=int(delta/Y)
    for i in range(1, tim+2):
        if Z > Xmin:
            Xmin += Y
        else:
            return Xmin

def covrange(X):
    if float(X)==0:
        return '0000'
    elif float(X) > 0 and  float(X) < 9:
        X+=1
        return '000'+str(int(X))
    elif float(X) == 9:
        return '0009'
    elif float(X) > 9 and  float(X) < 10:
        return '0010'
    elif float(X) >= 10 and float(X) < 20:
        return '00'+str(intervalue(10, 20, 2, X))
    elif float(X) >= 20 and float(X) < 50:
        return '00'+str(intervalue(20, 50, 3, X))
    elif float(X) >= 50 and float(X) < 100:
        return '00'+str(intervalue(50, 100, 5, X))
    elif float(X) >= 100 and  float(X) < 300:
        return '0'+str(intervalue(100, 300, 10, X))
    elif float(X) >= 300 and  float(X) < 700:
        return '0'+str(intervalue(300, 700, 20, X))
    elif float(X) >= 700 and  float(X) < 1000:
        return '0'+str(intervalue(700, 1000, 30, X))
    elif float(X) >= 1000 and  float(X) < 2000:
        return str(intervalue(1000, 2000, 50, X))
    elif float(X) >= 2000 and  float(X) < 5000:
        return str(intervalue(2000, 5000, 100, X))
    else:
        return '10000'

def dcovrange(X):
    if float(X)==0:
        return '0000'
    elif float(X) > 0 and float(X) < 9:
        X+=1
        return '000'+str(int(X))
    elif float(X) == 9:
        return '0009'
    elif float(X) > 9 and float(X) < 20:
        X+=1
        return '00'+str(int(X))
    elif float(X) >= 20 and float(X) < 100:
        return '00'+str(intervalue(20, 100, 2, X))
    elif float(X) >= 100 and  float(X) < 200:
        return '0'+str(intervalue(100, 200, 2, X))
    elif float(X) >= 200 and  float(X) < 1000:
        return '0'+str(intervalue(200, 1000, 5, X))
    elif float(X) >= 1000 and  float(X) < 2000:
        return str(intervalue(1000, 2000, 10, X))
    elif float(X) >= 2000 and  float(X) < 5000:
        return str(intervalue(2000, 5000, 50, X))
    else:
        return '10000'

def PE_tracker(sam_file, output_name):
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

def parse_checkm(checkm_containing_folder, pwd):
    bins_checkm={}
    os.chdir(pwd+'/'+checkm_containing_folder)
    for root, dirs, files in os.walk(pwd+'/'+checkm_containing_folder):
        for file in files:        
            if 'bin_stats_ext.tsv' in file: 
                for line in open(file,'r'):
                    genome_ids=str(line).strip().split('\t')[0]
                    bins_checkm[genome_ids]={}
                    try:
                        marker_lineage=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                        bins_checkm[genome_ids]['marker lineage']=marker_lineage
                    except:
                        print('marker lineage error')
                        bins_checkm[genome_ids]['marker lineage']=root

                    try:    
                        completeness=str(line).strip().split('\'Completeness\': ')[1].split(', ')[0]
                        bins_checkm[genome_ids]['Completeness']=float(completeness)
                    except:
                        bins_checkm[genome_ids]['Completeness']=0

                    try:
                        genome_size=str(line).strip().split('\'Genome size\':')[1].strip().split(', ')[0]
                        bins_checkm[genome_ids]['Genome size']=coerce_count(genome_size)
                    except:
                        bins_checkm[genome_ids]['Genome size']=0

                    bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0])
                    bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\': ')[1].split(',')[0].split('}')[0].strip())
                    bins_checkm[genome_ids]['Connections']=0 ### Be careful
    os.chdir(pwd)
    return bins_checkm

def mapping(assembly, group, datasets, num_threads, pwd):
    print('Building Bowtie2 index')
    os.system('bowtie2-build '+assembly+' '+assembly)
    print('Done!')
    print('-------------')

    print('Mapping datasets to contigs/scaffolds')
    f_coverage_matrix=open('Coverage_list_'+str(group)+'_'+assembly+'.txt', 'w')
    connections=[]
    for i in range(1, len(datasets)+1):
        print('Mapping '+str(assembly)+' with dataset '+str(i))
        os.system('bowtie2 -p '+str(num_threads)+' -x '+assembly+' -1 '+str(datasets[str(i)][0])+' -2 '+str(datasets[str(i)][1])+' -S '+str(group)+'_DNA-'+str(i)+'.sam -q --no-unal')
        os.system('samtools view -@ '+str(num_threads)+' -b -S '+str(group)+'_DNA-'+str(i)+'.sam -o '+str(group)+'_DNA-'+str(i)+'.bam')
        cytoscape=bundled_executable('Cytoscapeviz.pl', __file__)
        os.system('perl '+shlex.quote(cytoscape)+' -i '+shlex.quote(str(group)+'_DNA-'+str(i)+'.sam')+' -f 2 -a 150 -e 500 -m 3000 -c')
        # PE_tracker(str(group)+'_DNA-'+str(i)+'.sam', 'condensed.cytoscape.connections_'+str(group)+'_DNA-'+str(i)+'.tab')
        os.system('mv condensed.cytoscape.connections.tab condensed.cytoscape.connections_'+str(group)+'_DNA-'+str(i)+'.tab')
        connections.append('condensed.cytoscape.connections_'+str(group)+'_DNA-'+str(i)+'.tab')
        os.system('rm '+str(group)+'_DNA-'+str(i)+'.sam')
        ### py2
        os.system('samtools sort -@ '+str(num_threads)+' '+str(group)+'_DNA-'+str(i)+'.bam '+str(group)+'_DNA-'+str(i)+'_sorted') 
        try:
            with open(str(group)+'_DNA-'+str(i)+'_sorted.bam', 'r') as fh:
                pass
        except FileNotFoundError:
            print('samtools sorting '+str(group)+'_DNA-'+str(i)+'.bam failed. Re-do')
            ### py3
            os.system('samtools sort -@ '+str(num_threads)+' -o '+str(group)+'_DNA-'+str(i)+'_sorted.bam '+str(group)+'_DNA-'+str(i)+'.bam' )        
        f_coverage_matrix.write(str(group)+'_'+assembly+'_coverage_list_DNA-'+str(i)+'.txt'+'\n')

        if i == 1:
            bam_sorted=str(group)+'_DNA-1_sorted.bam'
        else:
            bam_sorted+=' '+str(group)+'_DNA-'+str(i)+'_sorted.bam'
    
    f_coverage_matrix.close()
    print('Mapping Done!')
    print('-------------')

    print('Scorting SAM file(s)')
    print('CMD: jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+ str(bam_sorted))
    os.system('jgi_summarize_bam_contig_depths --outputDepth '+str(group)+'_assembly.depth.txt '+str(bam_sorted))
    # os.system('rm '+str(bam_sorted))
    # os.system('cp '+str(group)+'_assembly.depth.txt '+str(pwd)+'/'+str(assembly))

    fout=open('Coverage_matrix_for_binning_'+str(assembly)+'.txt', 'w')
    n, title, cov=0, {}, {}
    for line in open(str(group)+'_assembly.depth.txt', 'r'):
        n+=1
        if n == 1:
            num_cov_groups=depth_sample_count(line)+2
            title['Name']='Length'+'\t'+'totalCoverage'+'\t'+'avgCoverage'
            for i in range(2, num_cov_groups):
                m=title['Name']
                title['Name']=m+'\t'+'Coverage'+str(i-1)+'\t'+'Cov'+str(i-1)+'range'+'\t'+'Cov'+str(i-1)+'drange'
                if i == 2:
                    covgs='Length'+'\t'+'totalCoverage'+'\t'+'avgCoverage'+'\t'+'Coverage1'+'\t'+'Cov1range'+'\t'+'Cov1drange'
                else:
                    covgs+='\t'+'Coverage'+str(i-1)+'\t'+'Cov'+str(i-1)+'range'+'\t'+'Cov'+str(i-1)+'drange'
            fout.write('Name'+'\t'+str(title['Name'])+'\n')
        else:
            totalcov=str(line).strip().split('\t')[2]
            num=float(num_cov_groups)-2
            avgcov=round(float(totalcov)/float(num), 3)
            cov[str(line).strip().split('\t')[0]]=str(line).strip().split('\t')[1]+'\t'+str(totalcov)+'\t'+str(avgcov)
            for i in range(2, num_cov_groups):
                m=cov[str(line).strip().split('\t')[0]]
                covi=str(line).strip().split('\t')[i*2-1]
                covr=covrange(float(covi))
                covdr=dcovrange(float(covi))
                cov[str(line).strip().split('\t')[0]]=m+'\t'+str(covi)+'\t'+str(covr)+'\t'+str(covdr)

    for item in cov.keys():
        fout.write(str(item)+'\t'+str(cov[str(item)])+'\n')
    fout.close()

    for i in range(1, len(datasets)+1):
        f=open(assembly+'_coverage_list_DNA-'+str(i)+'.txt', 'w')
        n=0
        for line in open(str(group)+'_assembly.depth.txt', 'r'):
            n+=1
            if n > 1:
                ids=str(line).strip().split('\t')[0]
                coverage=str(line).strip().split('\t')[2*int(i)+1]
                f.write(str(ids)+'\t'+str(coverage)+'\n')
            else:
                continue
        f.close()
    print('Done with generation of depth file!')
    print('-------------')
    # return 'Coverage_matrix_for_binning_'+str(assembly)+'.txt'
    return 'Coverage_matrix_for_binning_'+str(assembly)+'.txt', connections, str(group)+'_assembly.depth.txt'

def ModifyEnd(filename, n): ### Adding end to fasta file    
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

def modification_fa(extra_binset, start_index, pwd):
    print('Start feeding')

    fmod=open(str(start_index)+'_mod_bin.txt','w')
    for i in range(len(extra_binset)):
        binset=extra_binset[i]
        assembly_num=i+int(start_index)
        try:
            os.system('mkdir '+str(assembly_num)+'_'+binset)
        except:
            print(str(assembly_num)+'_'+binset, 'existed. Re-created the folder')
            os.system('rm -rf '+str(assembly_num)+'_'+binset)
            os.system('mkdir '+str(assembly_num)+'_'+binset)

    mod_assembly_list, mod_extra_binset = [], []
    frc=open('Contig_rename.txt','w')
    for i in range(len(extra_binset)):
        binset=extra_binset[i]
        mod_seq_id, n = {}, 0
        assembly_num=i+int(start_index)
        mod_assembly_list.append(str(assembly_num)+'_assembly.fa')
        mod_extra_binset.append(str(assembly_num)+'_'+binset)

        os.chdir(pwd+'/'+binset)
        x, assembly_seq, bin_id_rec = 0, {}, {}
        for root, dirs, files in os.walk(pwd+'/'+binset):
            for file in files:
                file_name_list=str(file).split('.')
                hz=file_name_list[-1]
                file_name_list.pop()
                file_name='.'.join(file_name_list)
                if 'fa' in hz or 'fasta' in hz or 'fna' in hz:
                    x+=1
                    bin_id=str(assembly_num)+'_genomes.'+str(x)
                    bin_id_rec[bin_id+'.fa']=''
                    fmod.write(str(file)+'\t'+str(bin_id)+'.fa'+'\n')
                    xxx=0
                    for record in SeqIO.parse(file,'fasta'):
                        xxx+=1
                        try:
                            if len(record.seq) != 0:
                                assembly_seq[str(record.seq)]+='||'+str(bin_id)+'_'+str(xxx)
                        except:
                            if len(record.seq) != 0:
                                assembly_seq[str(record.seq)]=str(bin_id)+'_'+str(xxx)

        os.chdir(pwd+'/'+str(assembly_num)+'_'+binset)
        for item in bin_id_rec.keys():
            f=open(item,'w')
            f.close()

        for seq in assembly_seq.keys():
            seq_id=assembly_seq[seq]
            if '||' in seq_id:
                bin_id_list=str(seq_id).split('||')
                for bin_id_seq in bin_id_list:
                    bin_id_list2=bin_id_seq.split('_')
                    bin_id_list2.pop()
                    bin_id_name='_'.join(bin_id_list2)
                    bin_id=bin_id_name+'.fa'
                    f=open(bin_id,'a')
                    f.write('>'+str(seq_id)+'\n'+str(seq)+'\n')
                    f.close()
            else:
                bin_id_list=str(seq_id).split('_')
                bin_id_list.pop()
                bin_id_name='_'.join(bin_id_list)
                bin_id=bin_id_name+'.fa'
                f=open(bin_id,'a')
                f.write('>'+str(seq_id)+'\n'+str(seq)+'\n')
                f.close()
        os.chdir(pwd)

        f_ass=open(str(assembly_num)+'_assembly.fa','w')
        for item in assembly_seq.keys():
            f_ass.write('>'+str(assembly_seq[item])+'\n'+str(item)+'\n')
        f_ass.close()
            
    fmod.close()
    print(mod_extra_binset)
    print(mod_assembly_list)
    return mod_extra_binset, mod_assembly_list

def data_feeding(extra_binset, datasets, start_index, num_threads,
                 output_folder_name, qc, pe):
    """
    Integrate external binsets into BASALT and recompute coverage/quality.

    Parameters
    ----------
    extra_binset : list of str
        List of paths to external binset folders.
    datasets : dict
        Mapping ``sample_id -> [r1_fastq, r2_fastq]`` with paired-end reads.
    start_index : int
        Starting index used to enumerate new assemblies and binsets.
    num_threads : int
        Number of threads for mapping and depth computation.
    output_folder_name : str
        Name of the target folder where fed data and metrics are stored.
    qc : {'checkm', 'checkm2'}
        Quality checking tool to use for the new binsets.
    pe : {'y', 'n'}
        Flag indicating whether PE-tracking modifications should be applied.

    Returns
    -------
    None
        Writes coverage matrices, connection files, CheckM outputs and
        renamed binsets into ``output_folder_name``.
    """
    pwd=os.getcwd()
    os.system('mkdir '+output_folder_name)

    datasets_fq={}
    for item in datasets.keys():
        datasets_fq[item]=[]
        if pe == 'y':
            datasets_fq[item].append(ModifyEnd(datasets[item][0], 1))
            datasets_fq[item].append(ModifyEnd(datasets[item][1], 2))
        else:
            datasets_fq[item].append('PE_r1_'+str(datasets[item][0]))
            datasets_fq[item].append('PE_r2_'+str(datasets[item][1]))

    mod_extra_binset, mod_assembly_list=modification_fa(extra_binset, start_index, pwd)

    # mod_binset_folder_list, mod_bin_seq, mod_bin_list, i=[], [], [], 0
    # for binset in extra_binset:
    #     i+=1
    #     bs_id='e'+str(i)
    #     A=merge_bin(binset, bs_id, pwd)
    #     mod_binset_folder_list.append(A[0])
    #     mod_bin_seq.append(A[1])
    #     mod_bin_list.append(A[2])
    # print('----------')

    # coverage_list, bins_checkm_total = [], {}
    for i in range(0, len(mod_assembly_list)):
        print('Calculating depth of bins from '+str(mod_extra_binset[i]))
        x=i+int(start_index)
        bs_id=str(x)
        mapping_output=mapping(mod_assembly_list[i], bs_id, datasets_fq, num_threads, pwd)
        # coverage_list.append(mapping_output[0])
        connections=mapping_output[1]
        os.system('rm Coverage_list_*')
        PEC=cal_connections(connections)
        connections_total=open('condense_connections_'+str(mod_extra_binset[i])+'.txt', 'w')
        connections_total.write('node1'+'\t'+'interaction'+'\t'+'node2'+'\t'+'connections'+'\n')
        for item2 in PEC.keys():
            connections_total.write(str(item2)+'\t'+str(PEC[item2])+'\n')
        connections_total.close()

        print('----------')
        if qc == 'checkm':
            print('checking '+str(mod_extra_binset[i])+' bins with checkM')
            os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(mod_extra_binset[i])+' '+str(mod_extra_binset[i])+'_checkm')
            print('----------')
            bins_checkm=parse_checkm(str(mod_extra_binset[i])+'_checkm/storage', pwd)
            os.chdir(pwd+'/'+str(mod_extra_binset[i]))
            f=open(str(bs_id)+'_bin_stats_ext.tsv','w')
            for item in bins_checkm.keys():
                f.write(str(item)+'\t'+str(bins_checkm[item])+'\n')
            f.close()
            os.chdir(pwd)
        elif qc == 'checkm2':
            print('checking '+str(mod_extra_binset[i])+' bins with checkM2')
            os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(mod_extra_binset[i])+' -x fa -o '+str(mod_extra_binset[i])+'_checkm')
            print('----------')

            os.chdir(str(mod_extra_binset[i])+'_checkm')
            print('Parsing '+str(mod_extra_binset[i])+' checkm output')
            
            checkm, n = {}, 0
            for line in open('quality_report.tsv','r'):
                n+=1
                if n >= 2:
                    binID=str(line).strip().split('\t')[0].strip()
                    genome_size=str(line).strip().split('\t')[8].strip()
                    completeness=str(line).strip().split('\t')[1].strip()
                    contamination=str(line).strip().split('\t')[2].strip()
                    N50=str(line).strip().split('\t')[6].strip()

                    checkm[str(binID)]={}
                    checkm[str(binID)]['N50']=int(N50)
                    checkm[str(binID)]['Completeness']=float(completeness)
                    checkm[str(binID)]['Genome size']=coerce_count(genome_size)
                    checkm[str(binID)]['Contamination']=float(contamination)
            
            os.chdir(pwd+'/'+str(mod_extra_binset[i]))
            f=open(str(bs_id)+'_quality_report.tsv','w')
            f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
            for binID in checkm.keys():
                f.write(str(binID)+'\t'+str(checkm[str(binID)]['Genome size'])+'\t'+str(checkm[str(binID)]['Completeness'])+'\t'+str(checkm[str(binID)]['Contamination'])+'\t'+str(checkm[str(binID)]['N50'])+'\n')
            f.close()
            os.chdir(pwd)

        os.system('mv '+str(mapping_output[0])+' '+output_folder_name)
        os.system('mv '+str(connections)+' '+output_folder_name)
        os.system('mv condense_connections_'+str(mod_extra_binset[i])+'.txt '+output_folder_name)
        os.system('mv '+str(mapping_output[2])+' '+output_folder_name)
        os.system('mv '+str(mod_extra_binset[i])+' '+output_folder_name)
        os.system('mv '+str(mod_extra_binset[i])+'_checkm '+output_folder_name)
        os.system('mv '+str(x)+'_mod_bin.txt '+output_folder_name)
        os.system('mv '+str(x)+'_assembly.fa '+output_folder_name)

    for item in datasets.keys():
        os.system('mv '+str('PE_r1_'+str(datasets[item][0]))+' '+output_folder_name)
        os.system('mv '+str('PE_r2_'+str(datasets[item][1]))+' '+output_folder_name)
    
    os.system('rm *.bam *.bt2')
    print('Done!')
    # return bestbinset

if __name__ == '__main__': 
    extra_binset=['BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_12092023']
    datasets={'1':['SRR12358675_sub_1.fastq','SRR12358675_sub_2.fastq'], '2':['SRR12358676_sub_1.fastq','SRR12358676_sub_2.fastq']}
    num_threads=60
    start_index=500 ### default 1000
    output_folder_name='Data_feeded'
    qc='checkm' ### checkm or checkm2
    pe='y' ### y or n
    data_feeding(extra_binset, datasets, start_index, num_threads, output_folder_name, qc, pe)
