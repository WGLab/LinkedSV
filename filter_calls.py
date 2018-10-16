#!/usr/bin/env python

import os
import sys
from my_utils import *

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in_bedpe> <black_region_bed_file> <gap_region_bed_file> <alt_ctg_file> <out_file>'
argc  = 5 

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    filter_calls(args, dbo_args, endpoint_args)

def filter_calls(args, dbo_args, endpoint_args):

    sv_bedpe_file            = args.merged_bedpe_file

    black_region_bed_file    = args.black_region_bed_file

    gap_region_bed_file      = args.gap_region_bed_file

    low_mapq_region_bed_file = args.low_mapq_region_bed_file

    alt_ctg_file = args.alt_ctg_file

    out_file     = args.filter_bedpe_file

    bin_size     = 100

    black_reg_dict = read_black_reg_bed_file (black_region_bed_file, bin_size)

    gap_left_region_dict, gap_right_region_dict = read_gap_region_file (gap_region_bed_file, bin_size)

    low_mapq_left_region_dict, low_mapq_right_region_dict = read_low_mapq_region_file (low_mapq_region_bed_file,  bin_size)

    alt_chr_name_set = read_alternative_contig_file (alt_ctg_file)

    sv_bedpe_fp = open(sv_bedpe_file, 'r')

    out_fp = open(out_file, 'w')
    n_alt_chr = 0
    n_black_reg = 0
    n_gap = 0
    while 1:
        line = sv_bedpe_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        chr1 = line[0]
        pos1 = int(line[1])
        chr2 = line[3]
        pos2 = int(line[4])
        endtype1 = line[10]
        endtype2 = line[11]
        index1 = int(pos1 / bin_size)
        index2 = int(pos2 / bin_size)
        if chr1 != chr2: 
            sv_length = int(1e10)
        else:
            sv_length = abs(pos2 - pos1)

        if chr1 in alt_chr_name_set: 
            n_alt_chr += 1
            continue
        if chr2 in alt_chr_name_set: 
            n_alt_chr += 1
            continue

        if chr1 in black_reg_dict and index1 in black_reg_dict[chr1]: 
            n_black_reg += 1
            continue
        if chr2 in black_reg_dict and index2 in black_reg_dict[chr2]:
            n_black_reg += 1
            continue

        if sv_length < 200000 and endtype1 == '3p_end' and chr1 in gap_left_region_dict and index1 in gap_left_region_dict[chr1]: 
            n_gap += 1
            continue
        if sv_length < 200000 and endtype1 == '5p_end' and chr1 in gap_right_region_dict and index1 in gap_right_region_dict[chr1]: 
            n_gap += 1
            continue

        if sv_length < 200000 and endtype2 == '3p_end' and chr2 in gap_left_region_dict and index2 in gap_left_region_dict[chr2]: 
            n_gap += 1
            continue

        if sv_length < 200000 and endtype2 == '5p_end' and chr2 in gap_right_region_dict and index2 in gap_right_region_dict[chr2]: 
            n_gap += 1
            continue


        if sv_length < 100000 and endtype1 == '3p_end' and chr1 in low_mapq_left_region_dict and index1 in low_mapq_left_region_dict[chr1] and endtype2 == '3p_end' and chr2 in low_mapq_left_region_dict and index2 in low_mapq_left_region_dict[chr2]: 
            n_gap += 1
            continue

        if sv_length < 100000 and endtype1 == '5p_end' and chr1 in low_mapq_right_region_dict and index1 in low_mapq_right_region_dict[chr1] and endtype2 == '3p_end' and chr2 in low_mapq_left_region_dict and index2 in low_mapq_left_region_dict[chr2]: 
            n_gap += 1
            continue

        if sv_length < 100000 and endtype1 == '3p_end' and chr1 in low_mapq_left_region_dict and index1 in low_mapq_left_region_dict[chr1] and endtype2 == '5p_end' and chr2 in low_mapq_right_region_dict and index2 in low_mapq_right_region_dict[chr2]: 
            n_gap += 1
            continue

        if sv_length < 100000 and endtype1 == '5p_end' and chr1 in low_mapq_right_region_dict and index1 in low_mapq_right_region_dict[chr1] and endtype2 == '5p_end' and chr2 in low_mapq_right_region_dict and index2 in low_mapq_right_region_dict[chr2]: 
            n_gap += 1
            continue


        out_fp.write(tab.join(line) + endl)

    sv_bedpe_fp.close()
    out_fp.close()

    return

def read_gap_region_file(gap_region_bed_file, bin_size):

    gap_left_region_dict = dict()
    gap_right_region_dict = dict()

    if not os.path.isfile(gap_region_bed_file): 
        return gap_left_region_dict, gap_right_region_dict

    in_gap_bed_fp = open(gap_region_bed_file, 'r')

    while 1:
        line = in_gap_bed_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        chrm  = line[0]
        start = int(line[1])
        end   = int(line[2])
        if chrm not in gap_left_region_dict:
            gap_left_region_dict[chrm] = dict()
        if chrm not in gap_right_region_dict:
            gap_right_region_dict[chrm] = dict()

        for i in range(start - 5000, start+1000, bin_size):
            index = int(i / bin_size)
            gap_left_region_dict[chrm][index] = 1

        for i in range(end-1000, end+5000, bin_size):
            index = int(i / bin_size)
            gap_right_region_dict[chrm][index] = 1

    in_gap_bed_fp.close()
    return gap_left_region_dict, gap_right_region_dict

def read_low_mapq_region_file (low_mapq_region_bed_file,  bin_size):

    low_mapq_left_region_dict = dict()

    low_mapq_right_region_dict = dict()

    if not os.path.isfile(low_mapq_region_bed_file):
        return low_mapq_left_region_dict, low_mapq_right_region_dict

    low_mapq_region_bed_fp = open(low_mapq_region_bed_file, 'r')

    while 1:
        line = low_mapq_region_bed_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        chrm  = line[0]
        start = int(line[1])
        end   = int(line[2])

        if chrm not in low_mapq_left_region_dict:
            low_mapq_left_region_dict[chrm] = dict()
        if chrm not in low_mapq_right_region_dict:
            low_mapq_right_region_dict[chrm] = dict()

        for i in range(start - 10000, start+1000, bin_size):
            index = int(i / bin_size)
            low_mapq_left_region_dict[chrm][index] = 1

        for i in range(end-1000, end+10000, bin_size):
            index = int(i / bin_size)
            low_mapq_right_region_dict[chrm][index] = 1

    low_mapq_region_bed_fp.close()
    return low_mapq_left_region_dict, low_mapq_right_region_dict

def read_black_reg_bed_file(black_region_bed_file, bin_size):
    black_reg_dict = dict()
    if not os.path.isfile(black_region_bed_file): 
        return black_reg_dict

    in_black_reg_bed_fp = open(black_region_bed_file, 'r')
    while 1:
        line = in_black_reg_bed_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        chrm   = line[0]
        start = int(line[1])
        end   = int(line[2])
     
        if chrm not in black_reg_dict:
            black_reg_dict[chrm]    = dict()

        for i in range(start, end, bin_size):
            index = int(i / bin_size)
            black_reg_dict[chrm][index]    = 1 
    
    in_black_reg_bed_fp.close()
    return black_reg_dict



if __name__ == '__main__':
    main()
