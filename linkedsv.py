#!/usr/bin/env python

import subprocess
from my_utils import *
from extract_weird_reads import *
from cluster_reads import *
from get_high_coverage_regions import *
from global_distribution import *
from find_paired_bk import *
from quantify2bkcand import *
from merge_quantified_calls import *
from filter_calls import *


def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    gc.enable()

    args.global_distribution_calculated = False

    extract_barcode_from_bam(args, endpoint_args)

    extract_weird_reads(args, dbo_args, endpoint_args)

    detect_increased_fragment_ends(args, dbo_args, endpoint_args)

    ## find paired breakpoints ##
    if args.global_distribution_calculated == False: 
        estimate_global_distribution(args, dbo_args, endpoint_args, endpoint_args.bcd22_file, is_fast_mode = False)

    task = 'searching for paired breakpoints'
    if args.run_from_begining == False and check_file_exists(args.bk_cand_pair_file) == True:
        myprint('paired breakpoint file existed, skipped %s' % (task))
    else:
        myprint(task)
        find_paired_bk(args, dbo_args, endpoint_args)
        gc.collect()

    ## quantification ##
    task = 'quantifying SV candidates'
    if args.run_from_begining == False and check_file_exists(args.quantified_bk_pair_file) == True:
        myprint('quantified SV file existed, skipped %s' % (task))
    else:
        myprint(task)
        quantify2bkcand(args, dbo_args, endpoint_args)
        gc.collect()

    ## merge calls ##
    task = 'merging SV candidates'
    if args.run_from_begining == False and check_file_exists(args.merged_bedpe_file) == True:
        myprint('merged bedpe file existed, skipped %s' % (task))
    else:
        myprint(task)
        merge_quantified_calls(args, dbo_args, endpoint_args)
        gc.collect()

    ## filter calls ##

    filter_calls(args, dbo_args, endpoint_args)

    return

def detect_decreased_barcode_overlap (args, dbo_args, endpoint_args):

    gc.enable()
    bam = args.bam
    out_prefix = args.out_prefix
    min_mapq = args.min_mapq
    n_thread = args.n_thread
    bin_size = dbo_args.bin_size
    win_size = dbo_args.win_size
    output_bk_candidate_file = dbo_args.bk_file 
    min_distance = dbo_args.min_distance_for_peak_calling 
    
    bcd12_existed = True 
    for bcd12_file in dbo_args.bcd12_file_list:
        if check_file_exists(bcd12_file) == False:
            bcd12_existed = False 
            break

	###  count barcode overlap | output files: bcd11 and bcd12 ###
    if args.run_from_begining == False and bcd12_existed:
        myprint ('bcd12 files existed, skipped counting barcode overlap...')
    else:
        myprint ('start counting barcode overlap...')
        count_barcode_overlap(args, dbo_args, endpoint_args)
        myprint ('finished counting barcode overlap...')

    gc.collect()

	### calculate expected overlap | output files: bcd13 ###
    bcd13_existed = True
    for bcd13_file in dbo_args.bcd13_file_list:
        if check_file_exists(bcd13_file) == False:
            bcd13_existed = False
            break

    if args.run_from_begining == False and bcd13_existed:
        myprint ('bcd13 files existed, skipped calculating expected overlapped barcodes...')
    else:
        myprint ('start calculating expected overlapped barcodes...')
        cal_expected_overlap_bcd_cnt (args, dbo_args, endpoint_args) 
        myprint ('finished calculating expected overlapped barcodes...')

    gc.collect()

	### peak calling for overlap statistic | output files: output_bk_candidate_file ##
    if args.run_from_begining == False and check_file_exists(output_bk_candidate_file):
        myprint ('dbo peak file existed, skipped peak calling')
    else:
        myprint ('start peak calling...')
        peak_calling1 (args, dbo_args, endpoint_args)
        myprint ('finished peak calling...')

    gc.collect()
    return 
    

def detect_increased_fragment_ends(args, dbo_args, endpoint_args):

    gc.enable()

    ### 1 clustering reads | output file: bcd22 file
    task = 'clustering reads'

    if args.run_from_begining == False and check_file_exists (endpoint_args.bcd22_file):
        myprint('bcd22 file existed, skipped %s' % (task))
    else:
        myprint(task)
        cluster_reads(args, dbo_args, endpoint_args) 
    
    gc.collect()

    ### 2 searching for extremely high coverage region  

    task = 'searching for extremely high coverage region'
    if args.run_from_begining == False and check_file_exists (endpoint_args.barcode_cov_file):
        myprint ('high coverage region file existed, skipped %s' % (task))
    else:
        myprint(task)
        get_high_coverage_regions(args, dbo_args, endpoint_args) 

    ### 3 estimating distribution parameters

    if args.global_distribution_calculated == False: 
        estimate_global_distribution(args, dbo_args, endpoint_args, endpoint_args.bcd22_file)

    output_arguments2file(args, dbo_args, endpoint_args)

    gc.collect()

    return

def extract_barcode_from_bam (args, endpoint_args):
    
    ## sort bam by barcode ##

    cmd = '%s %s | %s sort -m 3G -@ %d -t BX -o %s -' % (args.output_bam_coreinfo, args.bam, args.samtools, args.n_thread, args.sortbx_bam)

    if check_file_exists(args.sortbx_bam):
        myprint('File: %s existed, skipped sorting bam by barcode' % args.sortbx_bam)
    else:
        myprint('sorting bam file by barcode')
        myprint(cmd)
        os.system(cmd)


    ## extract barcode info ##
    cmd = args.extract_barcode + ' ' +  args.sortbx_bam + ' ' + endpoint_args.bcd21_file + ' ' + args.stat_file + ' ' + str(args.min_mapq) + endl

    if args.run_from_begining == False and check_file_exists(args.stat_file) and check_file_exists(endpoint_args.bcd21_file):
        myprint('File: %s existed, skipped extracting barcode from bam' % endpoint_args.bcd21_file)
    else:
        myprint('extracting barcode info from bam file')
        myprint(cmd)
        os.system(cmd)

    return

if __name__ == '__main__':
    main()
