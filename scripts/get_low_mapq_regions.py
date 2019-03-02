#!/usr/bin/env python

import os
import sys
from fragment import *
from my_utils import *
import numpy as np


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in.bcd21> <faidx_file> <out.txt>' 
argc  = 3 

def main():
    if len(arg) < argc:
        print usage
        sys.exit()

    bcd21_file = arg.pop(0)
    faidx_file = arg.pop(0)
    out_file = arg.pop(0)
    
    get_low_mapq_regions(bcd21_file, faidx_file, out_file)
    return

def get_low_mapq_regions(bcd21_file, faidx_file, out_file):

    bin_size = 10

    tid2chrname_list, chrname2tid_dict = get_chrnames(faidx_file)
    chr_len_list = get_chr_length(faidx_file)
    n_chr = len(chr_len_list)

    cov_dict = dict()
    for i in range(0, n_chr):
        chr_len = chr_len_list[i]
        n_bin = int(chr_len / bin_size) + 1
        cov_dict[i] = [0] * n_bin
    
    print 'reading bcd21 file'
    bcd21_fp = open(bcd21_file, 'r')
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        tid = int(line[0])
        start = int(line[1])
        end = int(line[2])
        start_bin_idx = int(start/bin_size)
        end_bin_idx   = int(end/bin_size)

        for i in range(start_bin_idx, end_bin_idx+1):
            cov_dict[tid][i] += 1
        
    bcd21_fp.close()

    print 'writing to the output file'
    out_fp = open(out_file, 'w')

    for tid in range(0, len(chr_len_list)): 

        cov_list = cov_dict[tid]
        pos_cov_list = list()
        for cov in cov_list: 
            if cov > 0: pos_cov_list.append(cov)

        if len(pos_cov_list) == 0: 
            continue

        median = np.median(pos_cov_list)
        std = median - np.percentile(pos_cov_list, 100.0-84.13) 
        lower_limit = min(median - 4 * std, median / 5.0)

        out_fp.write('#lower_limit=%.4f,median=%.4f,std=%.4f\n' % (lower_limit, median, std))

        for bin_idx in range(0, len(cov_list)):
            if cov_list[bin_idx] < lower_limit: out_fp.write('%s\t%d\t%d\t%d\n' % (tid2chrname_list[tid], bin_idx * bin_size, (bin_idx+1) * bin_size, cov_list[bin_idx]))

    out_fp.close()

    return 

if __name__ == '__main__':
    main()
