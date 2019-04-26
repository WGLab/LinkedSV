#!/usr/bin/env python

import os
import sys
from fragment import *
from my_utils import *
import numpy as np


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

def main():


    args, dbo_args, endpoint_args = parse_user_arguments()

    get_high_coverage_regions(args, dbo_args, endpoint_args)

    return

def get_high_coverage_regions(args, dbo_args, endpoint_args):

    bcd22_file = endpoint_args.bcd22_file
    faidx_file = args.faidx_file

    bin_size = 100

    chr_len_list = get_chr_length(faidx_file)
    n_chr = len(chr_len_list)

    frm_cov_dict = dict()
    for i in range(0, n_chr):
        chr_len = chr_len_list[i]
        n_bin = int(chr_len / bin_size) + 1
        frm_cov_dict[i] = [0] * n_bin
    
    bcd22_fp = open(bcd22_file, 'r')
    while 1:
        line = bcd22_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        frm = Fragment(line)
        if frm.tid not in frm_cov_dict: continue
        start_bin_idx = int (frm.start/bin_size)
        end_bin_idx = int(frm.end/bin_size)
        for i in range(start_bin_idx, end_bin_idx+1):
            if i < len(frm_cov_dict[frm.tid]): frm_cov_dict[frm.tid][i] += 1
        
    bcd22_fp.close()

    alt_tid_set = args.alt_tid_set

    barcode_cov_fp = open(endpoint_args.barcode_cov_file, 'w')
    high_cov_temp_file = args.out_prefix + '.high_cov_temp_file.bed'
    high_cov_temp_fp = open(high_cov_temp_file, 'w')

    for tid in range(0, len(chr_len_list)): 
        if tid in alt_tid_set: continue
        cov_list = frm_cov_dict[tid]
        pos_cov_list = list()
        for cov in cov_list: 
            if cov > 0: pos_cov_list.append(cov)
        if len(pos_cov_list) == 0: continue

        median = np.median(pos_cov_list)
        std = np.percentile(pos_cov_list, 84.13) - median

        upper_limit = median + 3 * std 
        pos_cov_list2 = list() 
        for cov in pos_cov_list:
            if cov < upper_limit: pos_cov_list2.append(cov)

        if len(pos_cov_list2) == 0: continue
        median = np.median(pos_cov_list2)
        upper_limit = int(median * 4 + 0.5)
        barcode_cov_fp.write('#upper_limit\t%d\n' % upper_limit)
        barcode_cov_fp.write('#tid\tstart\tend\tbarcode_coverage\n')

        for bin_idx in range(0, len(cov_list)):
            barcode_cov_fp.write('%d\t%d\t%d\t%d\n' % (tid, bin_idx * bin_size, (bin_idx+1) * bin_size, cov_list[bin_idx]))
            if cov_list[bin_idx] > upper_limit:
                high_cov_temp_fp.write('%d\t%d\t%d\t%d\n' % (tid, bin_idx * bin_size, (bin_idx+1) * bin_size, cov_list[bin_idx]))

    barcode_cov_fp.close()
    high_cov_temp_fp.close()

    cmd = 'bedtools merge -i %s > %s' % (high_cov_temp_file, endpoint_args.high_cov_file)
    os.system(cmd)
    if os.path.exists(high_cov_temp_file): os.remove(high_cov_temp_file)

    return 

def read_barcode_cov_file (barcode_cov_file):

    high_cov_dict = dict()
    barcode_cov_fp = open(barcode_cov_file, 'r')
    while 1:
        line = barcode_cov_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        if line[0] == '#upper_limit': 
            upper_limit = int(line[1])
            continue

        if len(line) < 3: continue

        tid = int(line[0])
        pos = int(line[1])
        cov = int(line[2])

        if cov < upper_limit: continue

        if tid not in high_cov_dict:  
            high_cov_dict[tid] = dict()

        high_cov_dict[tid][pos] = cov 

    barcode_cov_fp.close()
    
    return high_cov_dict

def in_high_cov_region(tid, pos, high_cov_dict): 
    
    pos_idx = int(pos/100) * 100
    if tid in high_cov_dict and pos_idx in high_cov_dict[tid]:
        return True
    else:
        return False

def read_high_cov_region_file(high_cov_file):

    high_cov_dict = dict()
    high_cov_R_dict = dict()
    high_cov_L_dict = dict()

    high_cov_fp = open(high_cov_file, 'r')
    while 1:
        line = high_cov_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        tid   = int(line[0])
        start = int(line[1])
        end   = int(line[2])
        
        if tid not in high_cov_dict:
            high_cov_dict[tid]    = dict()
            high_cov_R_dict[tid] = dict()
            high_cov_L_dict[tid] = dict()

        for i in range(start, end, 100):
            index = int(i / 100)
            high_cov_dict[tid][index]    = 1
            high_cov_R_dict[tid][index] = 1
            high_cov_L_dict[tid][index] = 1
            
        for i in range(end, end + 10000, 100):
            index = int(i / 100)
            high_cov_R_dict[tid][index] = 1

        for i in range(start-10000, start, 100):
            index = int(i / 100)
            high_cov_L_dict[tid][index] = 1

    return high_cov_dict, high_cov_R_dict, high_cov_L_dict

if __name__ == '__main__':
    main()
