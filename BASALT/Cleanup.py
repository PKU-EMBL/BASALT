#!/usr/bin/env python

"""
Helper utilities for cleaning up intermediate BASALT files.
"""

import os


def cleanup(assembly_list):
    """
    Remove intermediate index/coverage files and archive key matrices.

    Parameters
    ----------
    assembly_list : list
        List of assemblies. Currently unused but kept for API compatibility.
    """
    os.system('rm *.njs *.ndb *.nto *.ntf *.not *.nos')
    os.system('rm -f Checkm2_metrics_cache.tsv')
    os.mkdir('Coverage_depth_connection_SimilarBin_files_backup')
    os.system(
        'mv *.depth.txt Coverage_matrix_* Combat_* condense_connections_* '
        'Connections_* Similar_bins.txt Coverage_depth_connection_SimilarBin_files_backup'
    )
    os.system(
        'tar -zcvf Coverage_depth_connection_SimilarBin_files_backup.tar.gz '
        'Coverage_depth_connection_SimilarBin_files_backup'
    )
    os.system('rm -rf Coverage_depth_connection_SimilarBin_files_backup')
    os.system('rm -rf *_kmer bin_coverage Bin_coverage_after_contamination_removal bin_comparison_folder bin_extract-eleminated-selected_contig Bins_blast_output')
    os.system('tar -zcvf Group_comparison_files.tar.gz *_comparison_files')
    os.system('tar -zcvf Group_Bestbinset.tar.gz *_BestBinsSet')
    # os.system('tar -zcvf Group_checkm.tar.gz *_checkm')
    os.system('tar -zcvf Group_genomes.tar.gz *_genomes')
    os.system('rm -rf *_sr_bins_seq')
    os.system('tar -zcvf Binsets_backup.tar.gz BestBinse*') ###
    os.system('rm -rf *_comparison_files *_checkm *_genomes *_BestBinset *_BestBinsSet BestBinse* Deep_retrieved_bins coverage_deep_refined_bins S6_coverage_filtration_matrix S6_TNF_filtration_matrix split_blast_output TNFs_deep_refined_bins')
    os.system('rm *_checkpoint.txt')
    os.system('rm -rf Merged_seqs_*')
    os.system('rm -rf *.bt2 Outlier_in_threshold* Summary_threshold* Refined_total_bins_contigs.fa Total_bins.fa') 
    for i in range(1,20):
        os.system('rm -rf *_deep_retrieval_'+str(i))
    os.system('rm -rf *_MP_1 *_MP_2 *_gf_lr_polished *_gf_lr *_gf_lr_mod *_gf_lr_checkm *_long_read')
    os.system('rm Bin_reads_summary.txt Depth_total.txt Basalt_log.txt Assembly_mo_list.txt Assembly_MoDict.txt *_gf_lr_blasted.txt Bestbinset_list.txt Bin_extract_contigs_after_coverage_filtration.txt Bin_lw.txt Bin_record_error.txt')
    os.system('rm Bins_folder.txt BLAST_output_error.txt Concoct_* condensed.cytoscape* cytoscape.*')
    os.system('rm Hybrid_re-assembly_status.txt Mapping_log_* OLC_merged_error_blast_results.txt Potential_contaminted_seq_vari.txt Reassembled_bins_comparison.txt Rejudge_clean.txt')
    os.system('rm Remained_seq* Remapped_depth_test.txt Re-mapped_depth.txt Remapping.fasta TNFs_exceptional_contigs.txt Total_contigs_after_OLC_reassembly.fa')
    os.system('rm PE_r1_* PE_r2_*')
    for i in range(1, len(assembly_list)+1):
        os.system('rm '+str(i)+'_'+str(assembly_list[i-1]))

if __name__ == '__main__':
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    cleanup(assembly_list)
