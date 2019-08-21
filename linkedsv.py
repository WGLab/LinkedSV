#!/usr/bin/env python

import os
import sys
import subprocess
import bisect
import psutil
import gc
import math
import gzip

from scripts import my_utils 
from scripts import extract_weird_reads
from scripts import get_high_coverage_regions 
from scripts import global_distribution 
from scripts import find_paired_bk 
from scripts import quantify2bkcand 
from scripts import merge_quantified_calls 
from scripts import filter_calls
from scripts import arguments
#from scripts import visualize_sv_calls


tab = '\t'
endl = '\n'

FIX_LENGTH = int(1e10)
TimeFormat = '%m/%d/%Y %H:%M:%S'

def main():

    args, dbo_args, endpoint_args = arguments.parse_user_arguments()

    gc.enable()

    args.global_distribution_calculated = False

    extract_barcode_from_bam(args, endpoint_args)

    detect_increased_fragment_ends(args, dbo_args, endpoint_args)

    detect_decreased_barcode_overlap(args, dbo_args, endpoint_args)


    ## quantification ##
    task = 'quantifying SV candidates'
    if args.run_from_begining == False and my_utils.check_file_exists(args.quantified_bk_pair_file) == True:
        my_utils.myprint('quantified SV file existed, skipped %s' % (task))
    else:
        my_utils.myprint(task)
        quantify2bkcand.quantify2bkcand(args, dbo_args, endpoint_args)
    
    gc.collect()

    ## merge calls ##
    task = 'merging SV candidates'
    if args.run_from_begining == False and my_utils.check_file_exists(args.merged_bedpe_file) == True:
        my_utils.myprint('merged bedpe file existed, skipped %s' % (task))
    else:
        my_utils.myprint(task)
        merge_quantified_calls.merge_quantified_calls(args, dbo_args, endpoint_args)
    
    gc.collect()

    filter_calls.filter_calls(args, dbo_args, endpoint_args)

    gc.collect()

    image_out_dir = os.path.join(args.out_dir, 'images')
    my_utils.make_dir(image_out_dir)

    #visualize_sv_calls.visualize_sv_calls (args.filter_bedpe_file, dbo_args.bcd13_file, endpoint_args.bcd21_file, args.faidx_file, args.cal_2d_overlapping_barcodes, args.cal_read_depth_from_bcd21, image_out_dir, args.bam_name)

    ## remove temp files ##
    if args.rm_temp_files: 
        args.temp_file_list.append(endpoint_args.barcode_cov_file)
        args.temp_file_list.append(endpoint_args.high_cov_file)
        args.temp_file_list.append(args.args_file)
        args.temp_file_list.append(args.stat_file)
        args.temp_file_list.append(args.global_distribution_file)
        args.temp_file_list.append(args.node35_file)
        args.temp_file_list.append(args.node33_file)
        args.temp_file_list.append(args.node53_file)
        args.temp_file_list.append(args.node55_file)
        args.temp_file_list.append(args.node_cluster33_file)
        args.temp_file_list.append(args.node_cluster35_file)
        args.temp_file_list.append(args.node_cluster53_file)
        args.temp_file_list.append(args.node_cluster55_file)
        args.temp_file_list.append(args.quantified_bk_pair_file)
        args.temp_file_list.append(args.bk_cand_pair_file)
        args.temp_file_list.append(endpoint_args.bcd22_file)
        args.temp_file_list.append(args.weird_reads_file)
        args.temp_file_list.append(endpoint_args.low_mapq_bcd21_file)
        args.temp_file_list.append(endpoint_args.bcd21_file + '.split')
        args.temp_file_list.append(dbo_args.bcd11_file)
        args.temp_file_list.append(dbo_args.bcd12_file)


        for temp_file in args.temp_file_list:
            if os.path.exists(temp_file): os.remove(temp_file)
    
    return

def detect_decreased_barcode_overlap (args, dbo_args, endpoint_args):

    if args.is_wgs:
        win_size = 10000
    else:
        win_size = 40000
    
    ### calculating read depth | output file: args.read_depth_file
    task = 'calculating read depth'
    if args.run_from_begining == False and my_utils.check_file_exists(args.read_depth_file) == True:
        my_utils.myprint ('read depth file existed, skipped %s' % task)
    else:
        my_utils.myprint (task)
        cmd_args_list = [args.cal_read_depth_from_bcd21, endpoint_args.bcd21_file, args.read_depth_file, args.faidx_file, str(dbo_args.bin_size), str(args.min_mapq) ]
        my_utils.myprint('running command: %s' % (' '.join(cmd_args_list) ) ) 
        subprocess.call( cmd_args_list ) 
        my_utils.myprint ('finished %s' % task)
       
	### counting overlapping barcodes | output files dbo_args.bcd11_file
    task = 'counting overlapping barcodes between twin windows'
    if args.run_from_begining == False and my_utils.check_file_exists(dbo_args.bcd11_file) == True:
        my_utils.myprint ('bcd11 files existed, skipped %s' % task )
    else:
        my_utils.myprint (task)
        cmd_args_list = [args.cal_twin_win_bcd_cnt, endpoint_args.bcd21_file, dbo_args.bcd11_file, args.faidx_file, str(dbo_args.bin_size), str(win_size), str(args.min_mapq)]
        my_utils.myprint('running command: %s' % (' '.join(cmd_args_list) ) ) 
        subprocess.call( cmd_args_list )
        my_utils.myprint ('finished %s' % task)

    ### calculating centroid | output file: dbo_args.bcd12_file 
    task = 'calculating centroid'
    if args.run_from_begining == False and my_utils.check_file_exists(dbo_args.bcd12_file) == True:
        my_utils.myprint ('bcd12 files existed, skipped %s' % task )
    else:
        my_utils.myprint (task)
        cmd_args_list = [args.cal_centroid_from_read_depth, args.read_depth_file, dbo_args.bcd11_file, dbo_args.bcd12_file, args.faidx_file]
        my_utils.myprint('running command: %s' % (' '.join(cmd_args_list) ) ) 
        subprocess.call( cmd_args_list )
        my_utils.myprint ('finished %s' % task)

	### calculating expected overlap | output file: dbo_args.bcd13_file

    task = 'calculating barcode similarity and p-value'
    if args.run_from_begining == False and my_utils.check_file_exists(dbo_args.bcd13_file) == True:
        my_utils.myprint ('bcd12 files existed, skipped %s' % task )
    else:
        my_utils.myprint (task)
        if args.is_wgs:
            is_wgs = 1
        else:
            is_wgs = 0

        cmd_args_list = [args.cal_expected_overlap_value, dbo_args.bcd12_file, dbo_args.bcd13_file, str(is_wgs)]
        my_utils.myprint('running command: %s' % (' '.join(cmd_args_list) ) ) 
        subprocess.call( cmd_args_list ) 
        my_utils.myprint ('finished %s' % task)

    return 
    

def detect_increased_fragment_ends(args, dbo_args, endpoint_args):

    gc.enable()

    ### clustering reads | output file: bcd22 file
    task = 'clustering reads'

    if args.is_wgs: 
        is_wgs = 1
    else:
        is_wgs = 0

    if args.run_from_begining == False and my_utils.check_file_exists (endpoint_args.bcd22_file):
        my_utils.myprint('bcd22 file existed, skipped %s' % (task))
    else:
        my_utils.myprint(task)
    
        cmd = '%s %s %s %s %d %d %d' % (args.cluster_reads, endpoint_args.bcd21_file, endpoint_args.bcd22_file, args.weird_reads_file, is_wgs, args.user_defined_min_reads_in_fragment, args.min_mapq)
        my_utils.myprint(cmd)
        os.system(cmd)
    
    gc.collect()

    ### searching for extremely high coverage region  

    task = 'searching for extremely high coverage region'
    if args.run_from_begining == False and my_utils.check_file_exists (endpoint_args.barcode_cov_file):
        my_utils.myprint ('high coverage region file existed, skipped %s' % (task))
    else:
        my_utils.myprint(task)
        get_high_coverage_regions.get_high_coverage_regions(args, dbo_args, endpoint_args) 

    gc.collect()

    ### estimating distribution parameters

    if args.global_distribution_calculated == False: 
        global_distribution.estimate_global_distribution(args, dbo_args, endpoint_args, endpoint_args.bcd22_file)

    arguments.output_arguments2file(args, dbo_args, endpoint_args)

    gc.collect()

    ## find paired breakpoints ##

    task = 'searching for paired breakpoints'
    if args.run_from_begining == False and my_utils.check_file_exists(args.bk_cand_pair_file) == True:
        my_utils.myprint('paired breakpoint file existed, skipped %s' % (task))
    else:
        my_utils.myprint(task)
        find_paired_bk.find_paired_bk(args, dbo_args, endpoint_args)

    gc.collect()

    return

def extract_barcode_from_bam (args, endpoint_args):
    
    ## sort bam by barcode ##

    cmd = '%s %s | %s sort -l 1 -m 2G -@ %d -t BX -o %s -' % (args.output_bam_coreinfo, args.bam, args.samtools, args.n_thread, args.sortbx_bam)

    if my_utils.check_file_exists(args.sortbx_bam):
        my_utils.myprint('File: %s existed, skipped sorting bam by barcode' % args.sortbx_bam)
    else:
        my_utils.myprint('sorting bam file by barcode')
        my_utils.myprint('running command: %s' % cmd)
        os.system(cmd)

    ## extract barcode info ##
    n_compress_threads = args.n_thread - 1
    if n_compress_threads < 1: n_compress_threads = 1
    cmd = '%s %s __STDOUT__ %s | %s --fast --processes %d - > %s' % (args.extract_barcode, args.sortbx_bam, args.stat_file, args.pigz, n_compress_threads, endpoint_args.bcd21_file)

    if args.run_from_begining == False and my_utils.check_file_exists(args.stat_file) and my_utils.check_file_exists(endpoint_args.bcd21_file):
        my_utils.myprint('File: %s existed, skipped extracting barcode from bam' % endpoint_args.bcd21_file)
    else:
        my_utils.myprint('extracting barcode info from bam file')
        my_utils.myprint('running command: %s' % cmd)
        os.system(cmd)


    task = 'extracting low mapq bcd21'
    if args.run_from_begining == False and my_utils.check_file_exists(endpoint_args.low_mapq_bcd21_file) == True:
        my_utils.myprint('%s existed, skipped %s' % (endpoint_args.low_mapq_bcd21_file, task) )
    else:
        my_utils.myprint(task)
        get_low_mapq_bcd21_file(endpoint_args.bcd21_file, endpoint_args.low_mapq_bcd21_file, args.min_mapq) 

    if args.rm_temp_files and my_utils.check_file_exists(endpoint_args.bcd21_file):
        if os.path.exists(args.sortbx_bam): os.remove(args.sortbx_bam)

    return

def get_low_mapq_bcd21_file(bcd21_file, low_mapq_bcd21_file, min_mapq):

    if bcd21_file[-2:] == 'gz':
        bcd21_fp = gzip.open(bcd21_file, 'r')
    else:
        bcd21_fp = open(bcd21_file, 'r')

    low_mapq_bcd21_fp = gzip.open(low_mapq_bcd21_file, 'w')

    while 1:
        line = bcd21_fp.readline()
        if not line: break 
        if line[0] == '#':
            low_mapq_bcd21_fp.write(line)
            continue

        col_list = line.strip().split(tab)
        mapq = int(col_list[3])
        if mapq < min_mapq:
            low_mapq_bcd21_fp.write(line)

    low_mapq_bcd21_fp.close()
    bcd21_fp.close()

if __name__ == '__main__':
    main()
