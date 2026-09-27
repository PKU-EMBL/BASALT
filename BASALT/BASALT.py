#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Command-line interface for the BASALT metagenomic binning pipeline.

This script parses user arguments and dispatches to the corresponding
pipeline modules (autobinning, refinement, reassembly, data feeding,
and dereplication) depending on the selected mode and options.
"""

import time
import sys
import os
import argparse


def positive_integer(value):
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def percentage(value):
    parsed = int(value)
    if not 0 <= parsed <= 100:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return parsed


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='BASALT')
parser.add_argument('-a', '--assemblies', type=str, dest='assemblies',
                    help='List of assemblies, e.g.: as1.fa,as2.fa')
parser.add_argument('-s', '--shortreads', type=str, dest='sr_datasets',
                    help='List of paired-end reads, e.g.: r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq (paried_ends reads need \'/\' to seperate)')
parser.add_argument('-l', '--longreads', type=str, dest='long_datasets',
                    help='Including ont and pb dataset, excluding hifi dataset. List of long reads, e.g.: lr1.fq,lr2.fq')
parser.add_argument('-hf', '--hifi', type=str, dest='hifi_datasets',
                    help='Hifi dataset. List of hifi, e.g.: hf1.fq,hf2.fq')
# parser.add_argument('-c','--HIC', type=str, dest='Hi_C_dataset',
                    # help='List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq')
parser.add_argument('-t','--threads', type=positive_integer, dest='threads', default=4,
                    help='Number of threads, e.g.: 64')
parser.add_argument('-m','--ram', type=positive_integer, dest='ram', default=32,
                    help='Number of ram, minimum ram suggested: 32G')
parser.add_argument('-e','--extra_binner', type=str, dest='extra_binner',
                    help='Extra binner for binning: m: metabinner, v: vamb, l: lorbin; for instance: -e m, means BASALT will use metabinner for binning besides metabat2, maxbin2, and concoct')
parser.add_argument('-o','--out', type=str, dest='output_folder_name', default='Final_binset',
                    help='Name of the output folder. For binning, E.g. -o Anammox. BASALT would put those bins into folder Anammox_final_binset; for data feeding, e.g. -o Anammox; output files will under the folder of Anammox_data_feeded')
parser.add_argument('-q','--quality-check', choices=('checkm2', 'checkm'), dest='quality_check', default='checkm2',
                    help='Chance checkm version, default: checkm2; you may use: \'-q checkm\' to specify checkm for quality check when running BASALT')
parser.add_argument('--min-cpn', type=percentage, dest='Min_completeness', default=35,
                    help='Min completeness of kept bins (default: 35)')
parser.add_argument('--max-ctn', type=percentage, dest='Max_contamination', default=20,
                    help='Max contamination of kept bins (default: 20)')
parser.add_argument('--mode', choices=('new', 'continue'), dest='running_mode', default='continue',
                    help='Start a new project (new) or continue to run (continue). e.g. --mode continue / --mode new')
parser.add_argument('--module', choices=('autobinning', 'refinement', 'reassembly', 'all'), dest='functional_module', default=None,
                    help='Pipeline section: autobinning, refinement, reassembly, or all. A normal run defaults to all. With -r, omitting --module keeps standalone outlier screening. Pass --module refinement, reassembly, or all to continue an external binset through contig retrieval, secondary dereplication, rOLC, and reassembly. --module reassembly finishes any missing refinement stages first.')
# parser.add_argument('--autopara', type=str, dest='autobining_parameters', default='more-sensitive',
#                     help='Three parameters to chose: 1. more-sensitive; 2. sensitive; 3. quick. Default: more-sensitive. e.g. --autopara sensitive')
parser.add_argument('--refinepara', choices=('quick', 'deep'), dest='refinement_paramter', default='quick',
                    help='Two refinement parameters to chose: 1. deep; 2. quick. Default: quick. e.g. --refpara deep')
parser.add_argument('--sensitive', choices=('quick', 'sensitive', 'more-sensitive'), dest='binning_sensitive', default='sensitive',
                    help='Three parameters to chose: 1. quick; 2. sensitive; 3. more-sensitive. Default: sensitive. If you want to change the sensitive, use: e.g. --sensitive sensitive')
# parser.add_argument('--hybrid', type=str, dest='hybrid_reassembly', default='no',
#                     help='Use hybrid re-assembly. e.g. --hybrid y / --hybrid n; defalt no')
# parser.add_argument('--compression', type=str, dest='compression', default=0,
#                     help='Two refinement parameters to chose: 1. deep; 2. quick. Default: deep. e.g. --refpara quick')
# parser.add_argument('--hifi-only', action='store_true', default=False, dest='hifi',
#                     help='Only Hifi data')
# parser.add_argument('--ont-only', action='store_true', default=False, dest='ont',
#                     help='Only ont data')
# parser.add_argument('--hybrid', action='store_true', default=False, dest='ont',
#                     help='Hybrid data, including both HTS and ont/pb datasets')
parser.add_argument('-d','--data-feeding-folder', type=str, dest='data_feeding_folder',
                    help='List of folder name of extra binset(s) for data feeding, e.g.: -d binset1_folder_name,binset2_folder_name')
parser.add_argument('--binset-index', type=positive_integer, dest='extra_binset_start_index', default=500,
                    help='Optional parameter for data feeding. The start index of the extra binset, e.g.: -bi 5. BASALT already set a default index, but if you already had 4 assemblies, the binset start index could be 5')
# parser.add_argument('--only-refinement', action='store_true', dest='only_refinement',
#                     help='Only carry out refinement, e.g.: --only-refinement')
parser.add_argument('-r', '--refinement-binset', type=str, dest='refinement_binset', default='',
                    help='Existing binset folder. Without --module, runs outlier screening only and writes a resumable checkpoint. With --module refinement, reassembly, or all, continues through gap filling. Example: -r BestBinset --module all')
parser.add_argument('-c', '--coverage-list', type=str, dest='coverage_list',
                    help='List of depth file for refinement. Coverage file(s) could be generated from data feeding modole. e.g.: -c Coverage_matrix_for_binning_1_assembly.fa.txt,Coverage_matrix_for_binning_2_assembly.fa.txt')
parser.add_argument('-b', '--binsets-list', type=str, dest='binsets_list',
                    help='List of binsets for de-replication. Binset depth file(s) could be generated from data feeding modole. e.g.: -b 1_assembly_BestBinsSet,2_assembly_BestBinsSet')

args = parser.parse_args()
assemblies=args.assemblies
sr_datasets=args.sr_datasets
long_reads_list=args.long_datasets
hifi_list=args.hifi_datasets
# Hi_C_reads_list=args.Hi_C_dataset
QC_software=args.quality_check
num_threads=args.threads
ram=args.ram
extrabinner=args.extra_binner
min_cpn=args.Min_completeness
max_ctn=args.Max_contamination
continue_mode=args.running_mode
module_explicit=args.functional_module is not None
functional_module=args.functional_module or 'all'
# autobining_parameters=args.autobining_parameters
refinement_paramter=args.refinement_paramter
output_folder=args.output_folder_name
sensitivity=args.binning_sensitive
data_feeding_folder=args.data_feeding_folder
binsetindex=args.extra_binset_start_index
# only_refinement=args.only_refinement
refinement_binset=args.refinement_binset
coverage_list=args.coverage_list
binsets_list=args.binsets_list

# ---------------------------------------------------------------------------
# Parse and normalize input lists
# ---------------------------------------------------------------------------
def comma_list(value):
    return [item.strip() for item in value.split(',') if item.strip()] if value else []


assembly_list = comma_list(assemblies)
lr_list = comma_list(long_reads_list)
hifi_list = comma_list(hifi_list)
hic_list = []  # Hi-C is not exposed by the current CLI.
eb_list = comma_list(extrabinner)
data_feeding_folder = comma_list(data_feeding_folder)
coverage_list = comma_list(coverage_list)
binsets_list = comma_list(binsets_list)

invalid_binners = sorted(set(eb_list) - {'m', 'v', 'l'})
if invalid_binners:
    parser.error(
        "--extra_binner accepts only comma-separated m, v, or l values; got {}"
        .format(','.join(invalid_binners))
    )

datasets = {}
if sr_datasets:
    for index, sample in enumerate(sr_datasets.split('/'), start=1):
        mates = [mate.strip() for mate in sample.split(',')]
        if len(mates) != 2 or not all(mates):
            parser.error(
                "--shortreads requires R1,R2 for each sample and '/' between samples"
            )
        datasets[str(index)] = mates

if any(binner in {'v', 'l'} for binner in eb_list) and not datasets:
    parser.error("VAMB and LorBin adapters require paired short reads for BAM coverage")


# ---------------------------------------------------------------------------
# High-level configuration summary
# ---------------------------------------------------------------------------
pwd=os.getcwd()
print('Processing assemblies:', str(assembly_list))
print('Processing short-reads:', str(datasets))
print('Processing long-reads:', str(lr_list))
print('Processing hifi-reads:', str(hifi_list))
print('Processing Hi-C reads:', str(hic_list))
print('Output folder name will be:', str(output_folder))
print('Process with extra binner:', str(eb_list))
print('Quality check software:', str(QC_software))
print('Binning sensitivity:', str(sensitivity))
print('Processing with:', str(num_threads), 'threads')
print('Processing with:', str(ram), 'G')
print('Running status:', str(continue_mode))
if refinement_binset and not module_explicit:
    print('Binning module: outlier-only (pass --module all for gap filling)')
else:
    print('Binning module:', str(functional_module))
print('Min completeness:', str(min_cpn))
print('Max contamination:', str(max_ctn))
# print('Autobinning parameter:', str(autobining_parameters))
print('Refinement parameter:', str(refinement_paramter))
# print('Hybrid reassembly:', str(hybrid_reassembly))
print('Extra binset(s) for data feeding:', str(data_feeding_folder))
# print('Only refinement:', str(only_refinement))
print('Refinement binset:', str(refinement_binset))
print('List of coverage file(s):', str(coverage_list))
print('Binset(s) list:', str(binsets_list))

# Normalized continue mode used across downstream modules
if continue_mode == 'continue':
    continue_mode = 'last'

def dispatch_external_refinement():
    """Run outlier screening, or the checkpointed refinement/gap-filling route.

    Omitting ``--module`` preserves the historical outlier-only ``-r`` command.
    An explicit refinement, reassembly, or all module prepares the state files
    expected after autobinning and enters the same stages as a full run,
    without regenerating candidates. CheckM2 and legacy CheckM both use this
    dispatch. BASALT-Air is a separate codebase and is not invoked here.
    """
    from external_workflow import (
        OUTLIER_ONLY_NOTICE,
        ExternalBinsetError,
        external_module_plan,
        prepare_external_binset_route,
    )
    from basalt_runtime import append_checkpoint

    try:
        plan = external_module_plan(module_explicit, functional_module)
    except ValueError as exc:
        parser.error(str(exc))

    if plan == 'outlier':
        if QC_software == 'checkm2':
            from S5_Outlier_remover_DL_11012023 import outlier_remover_main
            outlier_remover_main(
                refinement_binset, coverage_list, datasets, lr_list, hifi_list,
                assembly_list, pwd, num_threads
            )
        else:
            from S5_Outlier_remover_DL_checkm import outlier_remover_main
            outlier_remover_main(
                refinement_binset, coverage_list, datasets, assembly_list,
                pwd, num_threads
            )
        append_checkpoint('Basalt_checkpoint.txt', '4th outlier removal done!')
        print(OUTLIER_ONLY_NOTICE)
        return

    try:
        state = prepare_external_binset_route(
            refinement_binset, assembly_list, coverage_list, functional_module,
            'new' if args.running_mode == 'new' else 'continue', pwd
        )
    except ExternalBinsetError as exc:
        parser.error(str(exc))

    print(
        'External binset route: '
        + state['plan']
        + ', resuming after checkpoint step '
        + str(state['checkpoint_step'])
    )
    orchestrator_module = state['orchestrator_module']
    if QC_software == 'checkm2':
        from BASALT_main_d import BASALT_main_d
        BASALT_main_d(
            assembly_list, datasets, num_threads, lr_list, hifi_list,
            hic_list, eb_list, ram, 'last', orchestrator_module,
            sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
            QC_software, output_folder, external_route=True
        )
        return

    if state['checkpoint_step'] < 7:
        from BASALT_main_c_refinement import BASALT_main_c_refinement
        BASALT_main_c_refinement(
            assembly_list, datasets, num_threads, lr_list, hifi_list,
            hic_list, eb_list, ram, 'last', orchestrator_module,
            sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
            QC_software, output_folder
        )
    if orchestrator_module == 'all':
        from BASALT_main_c_re_assembly import BASALT_main_c_re_assembly
        BASALT_main_c_re_assembly(
            assembly_list, datasets, num_threads, lr_list, hifi_list,
            hic_list, eb_list, ram, 'last', orchestrator_module,
            sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
            QC_software, output_folder
        )


def main():
    """
    Entry point for the BASALT command-line interface.

    This function parses command-line arguments (already configured in the
    module-level ``parser``), normalises them into Python data structures,
    and dispatches to the appropriate BASALT workflow depending on:

    - quality control backend (CheckM2 vs. CheckM),
    - selected functional module (autobinning / refinement / reassembly / all),
    - whether data feeding or dereplication on existing binsets is requested.

    Returns
    -------
    None
        Side effects include running the end-to-end BASALT pipeline and
        writing results to the specified output folder.
    """
    global output_folder
    if data_feeding_folder:
        if not datasets:
            parser.error("data feeding requires paired short reads via --shortreads")
    elif binsets_list:
        if not assembly_list or not coverage_list or not datasets:
            parser.error(
                "--binsets-list requires matching --assemblies, --coverage-list, "
                "and --shortreads inputs"
            )
        if not (
            len(binsets_list) == len(assembly_list) == len(coverage_list)
        ):
            parser.error(
                "--binsets-list, --assemblies, and --coverage-list must contain "
                "the same number of entries"
            )
    elif refinement_binset:
        if not assembly_list or not coverage_list or not datasets:
            parser.error(
                "--refinement-binset requires --assemblies, --coverage-list, "
                "and paired --shortreads"
            )
    else:
        if not assembly_list:
            parser.error("a normal BASALT run requires at least one assembly")
        if not datasets and not lr_list and not hifi_list:
            parser.error(
                "a normal BASALT run requires short reads, long reads, or HiFi reads"
            )

    # Main execution: select quality check implementation and functional module
    if QC_software == 'checkm2':
        if len(data_feeding_folder) != 0:
            # Data feeding for CheckM2 branch
            from Data_feeding import data_feeding
            if output_folder != 'Final_binset':
                output_folder = output_folder + '_data_feeded'
            else:
                output_folder = 'Data_feeded'
            data_feeding(
                data_feeding_folder, datasets, binsetindex,
                num_threads, output_folder, QC_software, pe='y'
            )

        elif len(binsets_list) != 0:
            from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023 import multiple_assembly_comparitor_main
            from external_workflow import DEREPLICATION_NOTICE
            step='initial_drep'
            # coverage_list: Coverage matrix
            # dataset: OK with both PE dataset or original datasets
            print(DEREPLICATION_NOTICE)
            multiple_assembly_comparitor_main(assembly_list, binsets_list, coverage_list, datasets, step, num_threads)

        elif refinement_binset != '':
            dispatch_external_refinement()

        else:
            pwd = os.getcwd()
            from BASALT_main_d import BASALT_main_d
            BASALT_main_d(
                assembly_list, datasets, num_threads, lr_list, hifi_list,
                hic_list, eb_list, ram, continue_mode, functional_module,
                sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                QC_software, output_folder
            )

    else:
        if len(binsets_list) != 0:
            from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import multiple_assembly_comparitor_main
            from external_workflow import DEREPLICATION_NOTICE
            step='initial_drep'
            # coverage_list: Coverage matrix
            # dataset: OK with both PE dataset or original datasets
            print(DEREPLICATION_NOTICE)
            multiple_assembly_comparitor_main(assembly_list, binsets_list, coverage_list, datasets, step, num_threads)

        elif refinement_binset != '':
            dispatch_external_refinement()

        else:
            if len(assembly_list) != 0:
                if functional_module == 'autobinning' or functional_module == 'all':
                    from BASALT_main_c_autobinning import BASALT_main_c_autobinning
                    BASALT_main_c_autobinning(
                        assembly_list, datasets, num_threads, lr_list, hifi_list,
                        hic_list, eb_list, ram, continue_mode, functional_module,
                        sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                        QC_software, output_folder
                    )

                if functional_module == 'refinement' or functional_module == 'all':
                    from BASALT_main_c_refinement import BASALT_main_c_refinement
                    BASALT_main_c_refinement(
                        assembly_list, datasets, num_threads, lr_list, hifi_list,
                        hic_list, eb_list, ram, continue_mode, functional_module,
                        sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                        QC_software, output_folder
                    )

                if functional_module == 'reassembly' or functional_module == 'all':
                    from BASALT_main_c_re_assembly import BASALT_main_c_re_assembly
                    BASALT_main_c_re_assembly(
                        assembly_list, datasets, num_threads, lr_list, hifi_list,
                        hic_list, eb_list, ram, continue_mode, functional_module,
                        sensitivity, refinement_paramter, max_ctn, min_cpn, pwd,
                        QC_software, output_folder
                    )

                if len(data_feeding_folder) != 0:
                    pwd = os.getcwd()
                    from BASALT_main_c_datafeeding import data_feeding_main
                    data_feeding_main(
                        assembly_list, datasets, num_threads, data_feeding_folder,
                        pwd, QC_software, output_folder, binsetindex, continue_mode
                    )

                from Cleanup import cleanup
                cleanup(assembly_list)
                print('All accomplish!')

            else:
                if len(data_feeding_folder) != 0:
                    from Data_feeding import data_feeding
                    if output_folder != 'Final_binset':
                        output_folder = output_folder + '_data_feeded'
                    else:
                        output_folder = 'Data_feeded'
                    pe = 'y'

                    data_feeding(
                        data_feeding_folder, datasets, binsetindex,
                        num_threads, output_folder, QC_software, pe
                    )


if __name__ == '__main__':
    main()

