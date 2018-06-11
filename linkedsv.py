#!/usr/bin/env python

import subprocess
from my_utils import *
from count_barcode_overlap import *
from peak_calling1 import *
from peak_calling2 import *
from cluster_read import *
from cal_endpoint_density import * 
from cal_expected_overlap_value import *
from global_distribution import *
from find_paired_bk import *
from refine_breakpoints import *
from quantify2bkcand3 import *
from merge_quantified_calls import *


def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    gc.enable()

    ## extract barcode information from bam file ##
    if args.run_from_begining == False and check_file_exists(args.bcd_file) == True:
        myprint('bcd file existed, skipped extracting barcode information from bam')
    else:
        extract_barcode_from_bam(args)

    ## scan for breakpoint candidates ##
    method1(args, dbo_args, endpoint_args)
    gc.collect()

    method2(args, dbo_args, endpoint_args)
    gc.collect()

    ## find paired breakpoints ##
    estimate_global_distribution(args, dbo_args, endpoint_args)
    if args.run_from_begining == False and check_file_exists(args.bk_cand_pair_file) == True:
        myprint('bk_cand_pair file existed, skipped finding paired breakpoints')
    else:
        find_paired_bk(args, dbo_args, endpoint_args)
        gc.collect()
    
    ## quantification ##
    if args.run_from_begining == False and check_file_exists(args.quantified_bk_pair_file) == True:
        myprint('quantified_bk_pair file existed, skipped quantification of sv candidates')
    else:
        quantify2bkcand(args, dbo_args, endpoint_args)
        gc.collect()

    ## refine breakpoints ##
    #if args.run_from_begining == False and check_file_exists(args.refinedbedpe_file) == True:
    #    myprint('refined bedpe file existed, skipped refinement of sv candidates')
    #else:
    #    refine_breakpoints(args, dbo_args, endpoint_args)
    #    gc.collect()

    ## merge calls ##
    if args.run_from_begining == False and check_file_exists(args.merged_bedpe_file) == True:
        myprint('merged bedpe file existed, skipped merging sv candidates')
    else:
        merge_quantified_calls(args, dbo_args, endpoint_args)
        gc.collect()

    if check_file_exists(args.merged_bedpe_file):
        remove_file(args.bcd_file)
        remove_file(endpoint_args.bcd21_file)
        remove_file(endpoint_args.tmpbcd22_file)
        remove_file(endpoint_args.bcd22_file)
        remove_file(args.bcd_file_of_target_region)
        remove_file(args.args_file)
        remove_file(args.global_distribution_file)

        remove_file(args.node33_file) 
        remove_file(args.node35_file) 
        remove_file(args.node53_file) 
        remove_file(args.node55_file) 
        remove_file(args.node_cluster33_file)
        remove_file(args.node_cluster53_file)
        remove_file(args.node_cluster35_file)
        remove_file(args.node_cluster55_file)

        remove_file(args.bk_cand_pair_file)
        remove_file(args.quantified_bk_pair_file)
    return

def method1(args, dbo_args, endpoint_args):

    if args.all_to_all == True: return
    if args.only_method2 == True:
        myprint ('skipped method1')
        return

    myprint ('start running breakpoint detection method 1')
    detect_decreased_barcode_overlap(args, dbo_args, endpoint_args)
    gc.collect()
    myprint ('finished breakpoint detection method 1')
    return

def method2(args, dbo_args, endpoint_args): 

    if args.only_method1 == True:
        myprint ('skipped method2')
        return 

    myprint ('start breakpoint detection method 2')
    detect_increased_fragment_ends (args, dbo_args, endpoint_args)
    myprint ('finished breakpoint detection method 2')
        
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
    ### 1 sorting bcd file | output file: bcd21
    if args.run_from_begining == False and check_file_exists(endpoint_args.bcd21_file):
        myprint ('bcd21 file existed, skipped sorting')
    else:
        myprint ('sorting barcode...')
        cmd = args.sort_barcode + ' ' + args.bcd_file + ' ' + args.out_prefix + ' ' + args.faidx_file
        os.system (cmd)
        myprint ('finished sorting barcode...')

    gc.collect()

    ### 2 clustering fragments | output file: bcd22 file
    if args.run_from_begining == False and check_file_exists(endpoint_args.bcd22_file):
        myprint ('bcd22 file existed, skipped clustering reads')
    else:
        myprint ('clustering reads...')
        cluster_reads(args, dbo_args, endpoint_args) 
        myprint ('finished clustering reads...')
    
    gc.collect()

    ### 3 estimate_global_distribution

    estimate_global_distribution(args, dbo_args, endpoint_args)
    output_arguments2file(args, dbo_args, endpoint_args)
    gc.collect()

    if args.all_to_all == True: return

    ### 4 calculate endpoint density  | output file: endpoints_density file
    if args.run_from_begining == False and check_file_exists (endpoint_args.endpoints_density_file):
        myprint ('endpoints_density file existed, skipped calculating endpoint density')
    else:
        myprint ('calculating endpoint density...')
        cal_endpoint_density(args, dbo_args, endpoint_args)
    
    gc.collect()

    ### 4 peak calling for endpoints | output file: bk.bed file
    if args.run_from_begining == False and check_file_exists(endpoint_args.bk_file):
        myprint ('endpoint bk file existed, skipped peak calling')
    else:
        myprint('start peak calling for endpoints')
        peak_calling2(args, dbo_args, endpoint_args)
        myprint('finished peak calling for endpoints')

    gc.collect()
    return 

def extract_barcode_from_bam (args):
    
    bcd_file = args.bcd_file
    stat_file = args.stat_file

    cmd = args.extract_barcode + ' ' +  args.bam + ' ' + bcd_file + ' ' + stat_file + ' ' + str(args.min_mapq) + endl

    if args.run_from_begining == False and check_file_exists(bcd_file) and check_file_exists(stat_file):
        myprint('bcd file existed: %s, skipped extracting barcode from bam' % bcd_file)
    else:
        myprint('extracting barcode from bam file...')
        myprint(cmd)
        os.system(cmd)
        myprint('finished extracting barcode from bam file...')

    return

if __name__ == '__main__':
    main()
