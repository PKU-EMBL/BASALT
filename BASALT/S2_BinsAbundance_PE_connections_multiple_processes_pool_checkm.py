#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Step S2 (CheckM branch): Bin abundance estimation and PE connections.

This module mirrors S2 for the CheckM-based pipeline, generating
coverage matrices and PE connection information compatible with
CheckM quality reports.
"""

from lib2to3.fixes import fix_buffer
from Bio import SeqIO
import sys, os, threading, glob
from multiprocessing import Pool
from collections import Counter
from time import ctime,sleep
from basalt_runtime import depth_sample_count


def intervalue(Xmin, Xmax, Y, Z):
    """
    Helper function to compute the upper bound of a coverage bin.

    Parameters
    ----------
    Xmin : float
        Lower bound of the range.
    Xmax : float
        Upper bound of the range.
    Y : float
        Step size within the range.
    Z : float
        Observed coverage value.

    Returns
    -------
    float
        Upper bound of the bin in which Z falls.
    """
    delta=Xmax-Xmin
    tim=int(delta/Y)
    for i in range(1, tim+2):
        if Z > Xmin:
            Xmin += Y
        else:
            return Xmin

def covrange(X):
    """
    Map a coverage value to a coarse-grained coverage range code.

    Parameters
    ----------
    X : float
        Coverage value.

    Returns
    -------
    str
        Four-character string representing the coverage range.
    """
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
    """
    Map a coverage value to a denser coverage range code for derivatives.

    Parameters
    ----------
    X : float
        Coverage value.

    Returns
    -------
    str
        Four-character string representing the coverage range for delta coverage.
    """
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

    
def CoverageMatrix(depth_file, assembly_name):
    """
    Build a coverage matrix for binning from a depth file (CheckM branch).

    Parameters
    ----------
    depth_file : str
        Path to the depth file produced from read mapping.
    assembly_name : str
        Name of the assembly (used in output filenames).

    Returns
    -------
    tuple
        (cov_dict, covgs_header, matrix_filename)
        where cov_dict maps contig_id -> coverage features,
        covgs_header is the header string, and matrix_filename
        is the name of the coverage matrix file.
    """
    path=os.getcwd()

    fout=open('Coverage_matrix_for_binning_'+str(assembly_name)+'.txt', 'w')

    n, title, cov=0, {}, {}
    for line in open(str(depth_file), 'r'):
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
            # if float(num) == 0:
            #     avgcov = 0.0
            #     print(f"Warning: zero length/count detected for a contig in {depth_file}")
            # else:
            #     avgcov = round(float(totalcov)/float(num), 3)            
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
    return cov, covgs, 'Coverage_matrix_for_binning_'+str(assembly_name)+'.txt'

def BinAbundance(depth, cov, covgs, output_format, name_of_the_binning_project,
                 path, genome_summary_dict, genome_summary_dict2, gen_sum):
    """
    Compute bin-level abundance statistics and write summary files (CheckM branch).

    See the non-CheckM S2 implementation for parameter and return value
    semantics; this variant writes CheckM-compatible summaries.
    """
    # os.chdir(path+'_genomes')
    genome_contig, genome_contig2={}, {}
    # try:
    #     os.mkdir(path+'_genomes/Original_bins')
    # except:
    #     print('Folder Original_bins already exist')

    for root,dirs,files in os.walk(path+'_genomes'):
        for file in files:
            # print('Reading', file
            hz=str(file).split('.')[-1]
            if str(hz) == 'fa' or str(hz) == 'fasta' or str(hz) == 'fna':
                if str(depth) not in str(file):
                    file_name_list=str(file).split('.')
                    file_name_list.remove(file_name_list[-1])
                    file_name='.'.join(file_name_list)
                    # file_name_qz=file_name.split('_genomes.')[0]
                    # file_name_hz=int(file_name.split('_genomes.')[1])
                    # file_name=file_name_qz+'_genomes.'+str(file_name_hz)
                    fout=open(str(file_name)+'_contigs_summary.txt', 'w')
                    genome_summary_dict[str(file_name)+'_contigs_summary.txt']={}
                    fout.write('Name'+'\t'+str(covgs)+'\t'+'GC%'+'\n')
                    for record in SeqIO.parse(path+'_genomes/'+str(file),'fasta'):
                        gc=int(str(record.seq).count("G"))+int(str(record.seq).count("C"))
                        gc_ratio=round(float(100*gc/float(len(record.seq))), 1)
                        fout.write(str(record.id)+'\t'+str(cov[str(record.id)])+'\t'+str(gc_ratio)+'%'+'\n')
                        genome_contig[str(record.id)]=str(file_name)
                        genome_contig2[str(record.id)]=str(file_name)
                    fout.close()
            else:
                continue
    
    print('----------------------------------------------------')
    print('Checking '+name_of_the_binning_project+' summary file')
    try:
        fx=open('Basalt_log.txt','a')
        fx.write('Checking '+name_of_the_binning_project+'summary file'+'\n')
        # fy=open('S2_checkpoint.txt','a')
    except:
        x=0

    # f=open(name_of_the_binning_project+'_summary.txt','w')
    for item in genome_summary_dict.keys():
        # f.write(item+'\n')
        total_cov_bin=0
        genomeID=str(item).split('_contigs_summary.txt')[0]
        n=0
        for line in open(str(item), 'r'):
            n+=1
            if n>=2:
                total_cov_bin+=float(str(line).strip().split('\t')[3])
            else:
                continue
            if n>= 2:
                genome_summary_dict2[str(genomeID)]=float(total_cov_bin)/float(n-1)
        os.system('mv '+str(item)+' '+path+'_genomes/')
        # os.system('mv '+str(item)+' '+path+'_genomes/Original_bins/')
    # f.close()

    gen_sum=sorted(genome_summary_dict2.items(), key=lambda genome_summary_dict2:genome_summary_dict2[1])
    # print(str(gen_sum))

    fout=open('Genome_summary_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('Prebin'+'\t'+'PreviousID'+'\t'+'avgCov'+'\n')
    genome_id, genome_avgcov, n={}, {}, 0
    for item in gen_sum:
        n+=1
        genome_id[str(item).split('\'')[1]]=n
        genome_avgcov[str(item).split('\'')[1]]=str(item).split(',')[1].split(')')[0].strip()
        fout.write(str(n)+'\t'+str(item).split('\'')[1].strip()+'\t'+str(item).split(',')[1].split(')')[0].strip()+'\n')
    fout.close()
    os.system('mv Genome_summary_'+str(name_of_the_binning_project)+'.txt '+path+'_genomes')

    for item in genome_contig.keys():
        if  genome_contig[str(item)] in genome_id.keys():
            m=str(genome_contig[str(item)])
            genome_contig[str(item)]=str(genome_id[str(genome_contig[str(item)])])+'\t'+str(m)+'\t'+str(genome_avgcov[str(genome_contig[str(item)])])+'\t'+'---'

    fout=open('Genome_contig_summary_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('ID'+'\t'+'Prebinid'+'\t'+'PreviousID'+'\t'+'avgCov'+'\t'+'EssCompleteness'+'\n')
    for item in genome_contig.keys():
        fout.write(str(item)+'\t'+str(genome_contig[str(item)])+'\n')
    fout.close()
    os.system('mv Genome_contig_summary_'+str(name_of_the_binning_project)+'.txt '+path+'_genomes')

    # os.chdir(path+'_genomes')
    # os.system('pwd')
    new_name_dict={}
    f=open('Bins_change_ID_'+name_of_the_binning_project+'.txt', 'w')
    for item in genome_id.keys():
        if 'metabat' in item:
            if 'tooShort' not in item:
                if 'lowDepth' not in item:
                    genome_name_qz=str(item).split('_genomes.')[0]
                    num=str(item).split('_genomes.')[1]
                    new_name=genome_name_qz+'_genomes.'+str(num)
                    f.write(item+'.fa changed to '+str(new_name)+'.fa'+'\n')
                    new_name_dict[item]=new_name
        elif 'maxbin2' in item or 'concoct' in item:
            # os.system('cp '+path+'_genomes/'+item+'.fasta '+path+'_genomes/Original_bins/')
            genome_name_qz=str(item).split('_genomes.')[0]
            num=str(item).split('_genomes.')[1]
            if num != '0':
                genome_name_hz=int(num)
                new_name=genome_name_qz+'_genomes.'+str(genome_name_hz)
                f.write(item+'.fasta changed to '+str(new_name)+'.fasta'+'\n')
                # f.write(item+'.fa changed to '+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fasta'+'\n')
                new_name_dict[item]=new_name
                os.system('mv '+path+'_genomes/'+item+'.fasta '+path+'_genomes/'+str(new_name)+'.fasta')
            else:
                os.system('mv '+path+'_genomes/'+str(item)+'.fasta '+path+'_genomes/'+str(item)+'_noclass.txt')
        elif 'metabinner' in item or 'vamb' in item  or 'semibin' in item or 'SingleContig' in item:
            genome_name_qz=str(item).split('_genomes.')[0]
            num=str(item).split('_genomes.')[1]
            new_name=genome_name_qz+'_genomes.'+str(num)
            f.write(item+'.fa changed to '+str(new_name)+'.fa'+'\n')
            new_name_dict[item]=new_name
    f.close()
    os.system('mv Bins_change_ID_'+name_of_the_binning_project+'.txt '+path+'_genomes')
    
    # for item in genome_id.keys():
    #     if 'metabat' in item:
    #         os.system('cp '+path+'_genomes/Original_bins/'+item+'.fa '+path+'_genomes/'+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fa')
    #     elif 'maxbin2' in item or 'concoct' in item:
    #         os.system('cp '+path+'_genomes/Original_bins/'+item+'.fasta '+path+'_genomes/'+name_of_the_binning_project+'_genomes.'+str(genome_id[item])+'.fasta')
    #     else:
    #         continue

    # os.system('tar zcvf Original_bins.tar.gz Original_bins')
    # os.system('rm -rf Original_bins')

    # try:
    #     os.mkdir(path+'_genomes/Summary')
    # except:
    #     print('Folder Summary already exist')

    ##### ?
    # print(str(genome_summary_dict))
    # for item in genome_summary_dict.keys():
    #     file_name_qz=item.split('_genomes.')[0]
    #     file_name_hz=int(item.split('_genomes.')[1].split('_contigs_summary.txt')[0])
    #     file_name=file_name_qz+'_genomes.'+str(file_name_hz)+'_contigs_summary.txt'
    #     os.system('mv '+str(item)+' '+path+'_genomes/'+str(file_name))
    ####

    checkm={}
    print('Reading '+name_of_the_binning_project+' checkm output')
    f_bin_checkm=open(name_of_the_binning_project+'_bin_stats_ext.tsv', 'w')
    for line in open(str(name_of_the_binning_project)+'_checkm/storage/bin_stats_ext.tsv','r'):
        binID=str(line).strip().split('{\'')[0].strip()
        genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
        taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
        completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
        contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].strip()
        GC=round(float(str(line).strip().split('\'GC\':')[1].split(',')[0].strip())*100, 1)
        Mean_scaffold_length=str(line).strip().split('Mean scaffold length\': ')[1].split(',')[0].split('}')[0].strip()
        checkm[str(binID)]=str(GC)+'\t'+str(genome_size)+'\t'+str(taxon)+'\t'+str(Mean_scaffold_length)+'\t'+str(completeness)+'\t'+str(contamination)
        # if str(binID) in genome_id.keys():
        #     line=line.replace(str(binID), name_of_the_binning_project+'_genomes.'+str(genome_id[binID]))
        if str(binID) in new_name_dict.keys():
            line=line.replace(str(binID), str(new_name_dict[str(binID)]))
            f_bin_checkm.write(line)

    for item in gen_sum:
        if str(item).split('\'')[1] not in checkm.keys():
            checkm[str(item).split('\'')[1]]='0'+'\t'+'0'+'\t'+'root'+'\t'+'0'+'\t'+'0'+'\t'+'0'
            f_bin_checkm.write(name_of_the_binning_project+'_genomes.'+str(genome_id[str(item).split('\'')[1]])+'\t'+'\'GC\': 0, \'GCN4\''+'\t'+'Genome size\': 0, \'Longest'+'\t'+'\'marker lineage\': root'+'\t'+'\'Completeness\': 0'+'\t'+'\'Contamination\': 0'+'\n')
    f_bin_checkm.close()
    os.system('cp '+name_of_the_binning_project+'_bin_stats_ext.tsv '+str(path)+'_checkm/storage/')
    os.system('mv '+name_of_the_binning_project+'_bin_stats_ext.tsv '+path+'_genomes')

    fout=open('prebinned_genomes_output_for_dataframe_'+str(name_of_the_binning_project)+'.txt', 'w')
    fout.write('ID'+'\t'+'Prebinid'+'\t'+'PreviousID'+'\t'+'avgCov'+'\t'+'EssCompleteness'+'\t'+'AvgGC'+'\t'+'GenomeSize'+'\t'+'Taxon'+'\t'+'MeanScaffoldLength'+'\t'+'Completeness'+'\t'+'Contamination'+'\n')
    for item in genome_contig.keys():
        fout.write(str(item)+'\t'+name_of_the_binning_project+'_genomes.'+str(genome_contig[str(item)])+'\t'+str(checkm[str(genome_contig2[str(item)])])+'\n')

    for item in cov.keys():
        if item not in genome_contig.keys():
            fout.write(str(item)+'\t'+name_of_the_binning_project+'_genomes.0'+'\t'+'unclustered'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'---'+'\t'+'uc'+'\t'+'uc'+'\n')
    fout.close()
    os.system('mv prebinned_genomes_output_for_dataframe_'+str(name_of_the_binning_project)+'.txt '+path+'_genomes')
    # os.chdir(path)
    return 'prebinned_genomes_output_for_dataframe_'+str(name_of_the_binning_project)+'.txt'

def GenerationOfGenomeGroupList(prebin_dataframe, PE_connection_file,
                                name_of_the_binning_project, pwd, path):
    """
    Partition bins into genome groups based on PE connections (CheckM branch).
    """
    print('---------------------------')
    print('Reading PE-connections file')

    contig_genome, m={}, 0
    for line in open(path+'_genomes/'+str(prebin_dataframe), 'r'):
        m+=1
        if m >= 2:
            contig_genome[str(line).strip().split('\t')[0]]=str(line).strip().split('\t')[1]

    genome_connection, m={}, 0
    for line in open(pwd+'/'+str(PE_connection_file), 'r'):
        m+=1
        if m >= 2:
            node1=str(line).strip().split('\t')[0]
            node2=str(line).strip().split('\t')[2]
            num_connections=str(line).strip().split('\t')[3]
            # print(str(node1), str(node2)
            if str(node1) in contig_genome.keys() and str(node2) in contig_genome.keys() and contig_genome[str(node1)] != contig_genome[str(node2)]:
                if contig_genome[str(node1)] not in genome_connection.keys():
                    genome_connection[contig_genome[str(node1)]]={}
                    genome_connection[contig_genome[str(node1)]][contig_genome[str(node2)]]=str(num_connections)
                else:
                    if contig_genome[str(node2)] not in genome_connection[contig_genome[str(node1)]].keys():
                        genome_connection[contig_genome[str(node1)]][str(contig_genome[str(node2)])]=str(num_connections)
                    else:
                        m1=genome_connection[contig_genome[str(node1)]][str(contig_genome[str(node2)])]
                        genome_connection[contig_genome[str(node1)]][contig_genome[str(node2)]]=str(int(m1)+int(num_connections))
        
                if contig_genome[str(node2)] not in genome_connection.keys():
                    genome_connection[contig_genome[str(node2)]]={}
                    genome_connection[contig_genome[str(node2)]][contig_genome[str(node1)]]=str(num_connections)
                else:
                    if contig_genome[str(node1)] not in genome_connection[contig_genome[str(node2)]].keys():
                        genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]=str(num_connections)
                    else:
                        m1=genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]
                        genome_connection[contig_genome[str(node2)]][str(contig_genome[str(node1)])]=str(int(m1)+int(num_connections))

    # print(str(genome_connection)

    genome_group='Genome_group_for_cytoscape_'+str(name_of_the_binning_project)+'.txt'
    genome_group_list='Genome_group_all_list_'+str(name_of_the_binning_project)+'.txt'
    fout=open(str(genome_group), 'w')
    fout2=open(str(genome_group_list), 'w')
    # fout3=open('Genome_group_all_top'+str(topGenome)+'.txt', 'w')
    fout.write('Genome1'+'\t'+'Connections'+'\t'+'Genome2'+'\n')
    fout2.write('Genome'+'\t'+'Connecting genomes'+'\n')
    for item in genome_connection.keys():
        fout2.write(str(item)+'\t'+str(genome_connection[item]).strip()+'\n')
        num=len(genome_connection[item])
        if num == 1:
            fout.write(str(item)+'\t'+str(genome_connection[item]).split(':')[1].split('\'')[1].strip()+'\t'+str(genome_connection[item]).split(':')[0].split('\'')[1].strip()+'\n')        
        else:
            lis=str(genome_connection[item]).split(',')
            for i in range(0, num):
                fout.write(str(item)+'\t'+str(lis[i]).split(':')[1].split('\'')[1].strip()+'\t'+str(lis[i]).split(':')[0].split('\'')[1].strip()+'\n')

    fout.close()
    fout2.close()

    genome_total_connection={}
    genome_total_connection_file='Bins_total_connections_'+str(name_of_the_binning_project)+'.txt'
    f=open(str(genome_total_connection_file), 'w')
    f.write('Bin'+'\t'+'Total_connections'+'\n')
    for item in genome_connection.keys():
        genome_total_connection[item]=0
        if len(genome_connection[item]) != 0:
            for i in genome_connection[item].keys():
                genome_total_connection[item]+=int(genome_connection[item][i])
        f.write(str(item)+'\t'+str(genome_total_connection[item])+'\n')
    f.close()
    os.system('mv '+genome_group+' '+genome_group_list+' '+genome_total_connection_file+' '+pwd+'/'+name_of_the_binning_project+'_genomes')

    print('Generation of Genome Group of '+str(name_of_the_binning_project)+' List Done!')
    try:
        fy=open('S2_checkpoint.txt','a')
    except:
        print('S2_checkpoint.txt did not found')
    
    try:
        fy.write(str(name_of_the_binning_project)+'\t'+'done!'+'\n')
    except:
        xyzt=0

    try:
        fx=open('Basalt_log.txt','a')
    except:
        print('Basalt_log.txt did not found')
    
    try:
        fx.write('Parsed '+str(name_of_the_binning_project)+'\t'+'done!'+'\n')
    except:
        xyzt=0
    return genome_total_connection_file

def multi_threads(pwd, depth_file, cov, covs, bin_format, bin_folder,
                  genome_summary_dict, genome_summary_dict2, gen_sum,
                  PE_connections_file):
    """
    Worker function to run BinAbundance and genome grouping for one bin folder.
    """
    path=pwd+'/'+bin_folder
    a=BinAbundance(depth_file, cov, covs, bin_format, bin_folder, path, genome_summary_dict, genome_summary_dict2, gen_sum)
    GenerationOfGenomeGroupList(a, PE_connections_file, bin_folder, pwd, path)
    # GenerationOfGenomeGroupList(BinAbundance(depth_file, cov, covs, bin_format, bin_folder, path, genome_summary_dict, genome_summary_dict2, gen_sum), PE_connections_file, bin_folder, pwd, path)

def binsabundance_pe_connections(assembly_binning_group, depth_files,
                                 PE_connections_files, assembly_names,
                                 num_threads):
    """
    Entry point for S2 (CheckM branch): compute bin abundance and PE connections.

    Parameters mirror the non-CheckM implementation but refer to
    CheckM-specific outputs and file naming.
    """
    print('-------------------------------')
    print('Processing Step2')
    print(str(assembly_binning_group))
    print(str(depth_files))
    print(str(PE_connections_files))
    print(str(assembly_names))
    print('-------------------------------')
    
    try:
        fx=open('Basalt_log.txt','a')
    except:
        fx=open('Basalt_log.txt','w')
    fx.write('-------------------------------'+'\n')
    fx.write('Processing Step2'+'\n'+str(assembly_binning_group)+'\n'+str(depth_files)+'\n'+str(PE_connections_files)+'\n'+str(assembly_names)+'\n')
    fx.write('-------------------------------'+'\n')
    # fx.close()

    try:
        fy=open('S2_checkpoint.txt','a')
    except:
        fy=open('S2_checkpoint.txt','w')
    # fy.close()
    
    parsed_folder=[]
    for line in open('S2_checkpoint.txt','r'):
        parsed_folder.append(line.strip().split('\t')[0].strip())

    bin_folder, coverage_matrix_list = {}, {}
    pwd=os.getcwd()
    pool=Pool(processes=num_threads)

    for item in assembly_binning_group.keys():
        coverage_matrix_list[item]=[]
        depth_file=depth_files[item]
        PE_connections_file=PE_connections_files[item]
        assembly_name=assembly_names[item]
        bins_folders_name_list=assembly_binning_group[item]

        coverage_matrix=CoverageMatrix(depth_file, assembly_name)
        coverage_matrix_list[item]=coverage_matrix[2]
        
        genome_summary, genome_summary2, gen_sum = {}, {}, {}
        for item2 in bins_folders_name_list:
            if item2 not in parsed_folder:
                fx.write('Parsing '+str(item2)+'\n')
                if 'maxbin2' in item2 or 'concoct' in item2:
                    genome_summary[item2], genome_summary2[item2], gen_sum[item2] ={}, {}, []
                    pool.apply_async(multi_threads,args=(pwd, depth_file, coverage_matrix[0], coverage_matrix[1], 'fasta', item2, genome_summary[item2], genome_summary2[item2], gen_sum[item2], PE_connections_file))
                elif 'metabat' in item2 or 'metabinner' in item2 or 'vamb' in item2 or 'semibin' in item2 or 'SingleContig' in item2:
                    genome_summary[item2], genome_summary2[item2], gen_sum[item2] ={}, {}, []
                    pool.apply_async(multi_threads,args=(pwd, depth_file, coverage_matrix[0], coverage_matrix[1], 'fa', item2, genome_summary[item2], genome_summary2[item2], gen_sum[item2], PE_connections_file))
                else:
                    genome_summary[item2], genome_summary2[item2], gen_sum[item2] ={}, {}, []
                    pool.apply_async(multi_threads,args=(pwd, depth_file, coverage_matrix[0], coverage_matrix[1], 'fa', item2, genome_summary[item2], genome_summary2[item2], gen_sum[item2], PE_connections_file))

            bin_folder[pwd+'/'+item2+'_genomes']=''
    pool.close()
    pool.join()

    fx.close()
    fy.close()

    # for item in bin_folder.keys():
    #     os.chdir(item)
    #     output=os.path.isfile('Original_bins.tar.gz')
    #     if output != 'True':
    #         os.system('tar zcvf Original_bins.tar.gz Original_bins')
    # #    except:
    # #        xxxxxx=0
    #     os.system('rm -rf Original_bins')

    #     os.chdir(pwd)
    genome_summary_list=glob.glob(r'*_contigs_summary.txt')
    for item in genome_summary_list:
        file_name_qz=item.split('_genomes.')[0]
        try:
            file_name_hz=int(item.split('_genomes.')[1].split('_contigs_summary.txt')[0])
            file_name=file_name_qz+'_genomes.'+str(file_name_hz)+'_contigs_summary.txt'
            folder_name=file_name_qz+'_genomes'
            os.system('mv '+str(item)+' '+pwd+'/'+str(folder_name)+'/'+str(file_name))
        except:
            os.system('rm '+str(item))
    # os.system('rm Bins_change_ID_*')
    return coverage_matrix_list
            
if __name__ == '__main__': 
    # assembly_binning_group={'1':['1_adh_dn1_contigs.fasta_0.3_maxbin2',  '1_adh_dn1_contigs.fasta_0.5_maxbin2', '1_adh_dn1_contigs.fasta_0.7_maxbin2', '1_adh_dn1_contigs.fasta_0.9_maxbin2',
    # '1_adh_dn1_contigs.fasta_200_metabat', '1_adh_dn1_contigs.fasta_300_metabat', '1_adh_dn1_contigs.fasta_400_metabat' , '1_adh_dn1_contigs.fasta_500_metabat',
    # '1_adh_dn1_contigs.fasta_300_concoct', '1_adh_dn1_contigs.fasta_400_concoct']
    # }
    assembly_binning_group={'1':['1_assembly_sample1.fa_metabinner']}
    # assembly_binning_group={'1':['1_assembly_sample1.fa_0.9_maxbin2', '1_assembly_sample1.fa_200_concoct', '1_assembly_sample1.fa_200_metabat']}

    # assembly_binning_group={'1':['1_adh_dn1_contigs.fasta_0.3_maxbin2',  '1_adh_dn1_contigs.fasta_0.5_maxbin2', '1_adh_dn1_contigs.fasta_0.7_maxbin2', '1_adh_dn1_contigs.fasta_0.9_maxbin2',
    # '1_adh_dn1_contigs.fasta_200_metabat', '1_adh_dn1_contigs.fasta_300_metabat', '1_adh_dn1_contigs.fasta_400_metabat' , '1_adh_dn1_contigs.fasta_500_metabat',
    # '1_adh_dn1_contigs.fasta_300_concoct', '1_adh_dn1_contigs.fasta_400_concoct'],
    # '2':['2_adh_dn2_contigs.fasta_0.3_maxbin2', '2_adh_dn2_contigs.fasta_0.5_maxbin2', '2_adh_dn2_contigs.fasta_0.7_maxbin2', '2_adh_dn2_contigs.fasta_0.9_maxbin2',
    # '2_adh_dn2_contigs.fasta_200_metabat', '2_adh_dn2_contigs.fasta_300_metabat', '2_adh_dn2_contigs.fasta_400_metabat', '2_adh_dn2_contigs.fasta_500_metabat', '2_adh_dn2_contigs.fasta_500_concoct']
    # }

    depth_files={'1':'1_assembly.depth.txt'}
    PE_connections_files={'1':'condense_connections_assembly_sample1.fa.txt'}
    assembly_names={'1':'1_assembly_sample1.fa'}
    num_threads=2

    # binsabundance_pe_connections(bins_folders_name_list, depth_file, PE_connections_file, assembly_name, num_threads)
    coverage_matrix_list=binsabundance_pe_connections(assembly_binning_group, depth_files, PE_connections_files, assembly_names, num_threads)
