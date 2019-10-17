#!/usr/bin/env python

import os
import sys
import gzip
from scipy.spatial import * 
import bisect
import math
import numpy as np

try:
    from scripts import plot_read_depth
except ImportError:
    import plot_read_depth
try:
    from scripts import my_utils
except ImportError:
    import my_utils
try:
    from scripts import bedpe
except ImportError:
    import bedpe
try:
    from scripts import arguments
except ImportError:
    import arguments


tab = '\t'
endl = '\n'

class Bcd21Core:

    def __init__(self, attr_list):
        self.tid, self.start, self.end, self.mapq, self.bcd = attr_list[0:5]
        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)
        self.mapq = int(self.mapq)

    def key_start(self):
        return self.tid * my_utils.FIX_LENGTH + self.start

    def key_end(self):
        return self.tid * my_utils.FIX_LENGTH + self.end

class D2:
    def __init__(self, pos1, pos2_list):
        self.pos1 = pos1
        self.pos2_list = pos2_list


def main():

    args, dbo_args, endpoint_args = arguments.parse_user_arguments()

    filter_calls(args, dbo_args, endpoint_args)

def filter_calls(args, dbo_args, endpoint_args):

    my_utils.myprint('filtering SV calls')

    bin_size = 100
    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(args.faidx_file)
    alt_chr_name_set = my_utils.read_alternative_contig_file (args.alt_ctg_file)

    if os.path.exists(args.black_region_bed_file):
        my_utils.myprint ('reading black region bed file: %s' % args.black_region_bed_file)
        black_reg_dict = read_black_reg_bed_file (args.black_region_bed_file, bin_size)
    else:
        if args.ref_version == 'hg19' or args.ref_version == 'hg38' or args.ref_version == 'b37':
            my_utils.myprint ('ERROR! black list file is missing: %s' % black_region_bed_file)
        black_reg_dict = dict()

    if os.path.exists(args.gap_region_bed_file):
        gap_left_region_dict, gap_right_region_dict = read_gap_region_file (args.gap_region_bed_file, bin_size)
        my_utils.myprint ('reading gap region bed file: %s' % args.gap_region_bed_file)
    else:
        if args.ref_version == 'hg19' or args.ref_version == 'hg38' or args.ref_version == 'b37':
            my_utils.myprint ('ERROR! gap region file is missing: %s' % args.gap_region_bed_file)
        gap_left_region_dict = dict()
        gap_right_region_dict = dict()

    raw_svcall_list = my_utils.read_object_file(args.merged_bedpe_file, bedpe.QuantifiedBKCandCore)
    for i in range(0, len(raw_svcall_list)):
        raw_svcall_list[i].ft = '.'

    round1_retained_sv_list = filter_1d_blacklist(raw_svcall_list, black_reg_dict, alt_chr_name_set, gap_left_region_dict, gap_right_region_dict, bin_size )
  
    n_retained_sv = 0
    for svcall in round1_retained_sv_list:
        if svcall.ft == '.': n_retained_sv += 1
    my_utils.myprint('number of retained SVs: %d' % n_retained_sv)
        
    round2_retained_sv_list = filter_low_mapq_gaps(round1_retained_sv_list, endpoint_args, chrname2tid_dict) 
  
    n_retained_sv = 0
    for svcall in round2_retained_sv_list:
        if svcall.ft == '.': n_retained_sv += 1

    if args.ref_version == 'b37':
        remove_chr_prefix = True
    else:
        remove_chr_prefix = False

    round3_retained_sv_list = filter_calls_2d (round2_retained_sv_list, args.black_region_2d_file, args.filter_bedpe_file, remove_chr_prefix)

    n_retained_sv = 0
    for svcall in round3_retained_sv_list:
        if svcall.ft == '.': n_retained_sv += 1
   
    round4_retained_sv_list = filter_dbo_score (round3_retained_sv_list, args)
  
    n_retained_sv = 0
    for svcall in round4_retained_sv_list:
        if svcall.ft == '.': n_retained_sv += 1

    round5_retained_sv_list = filter_read_depth (round4_retained_sv_list, args)
  
    round6_retained_sv_list = filter_sv_length (round5_retained_sv_list, args)

    final_retained_sv_list = round6_retained_sv_list

    n_retained_sv = 0
    for svcall in final_retained_sv_list:
        if svcall.ft == '.': n_retained_sv += 1
    my_utils.myprint('number of retained SVs: %d' % n_retained_sv)

    header  = '#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\t'
    header += 'sv_type\tsv_id\tsv_length\tqual_score\tfilter\tinfo\n'

    out_file  = args.filter_bedpe_file
    out_fp    = open(out_file, 'w')
    out_fp.write(header)
    sv_id = 0
    n_svcall = len(final_retained_sv_list)
    n_digit = int(math.log10(n_svcall) + 2)

    for svcall in final_retained_sv_list:
        if svcall.ft == '.':
            svcall.ft = 'PASS'
            sv_id += 1
            sv_id_str = str(sv_id)
            sv_id_str = '0' * (n_digit - len(sv_id_str)) + sv_id_str
            svcall.sv_id = 'ID%s' % sv_id_str
            out_fp.write(svcall.output_core2() + endl)

    out_fp.close()

    return

def filter_sv_length (input_sv_list, args):
    for j in range(0, len(input_sv_list)):
        svcall = input_sv_list[j]
        if svcall.chrm1 != svcall.chrm2: continue
        if int(svcall.svlength)  < 5000: 
            input_sv_list[j].ft = 'LENGTH_FILTER'

    return input_sv_list

def filter_read_depth(input_sv_list, args):

    if args.is_wgs == False or args.germline_mode == False:
        return input_sv_list 

    chr_len_list = my_utils.get_chr_length(args.faidx_file)

    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(args.faidx_file)

    wg_high_mapq_depth_list, wg_total_depth_list, bin_size = plot_read_depth.get_wg_depth_list(args.read_depth_file, chr_len_list)

    wg_depth_mean, wg_depth_std = get_mean_std_depth(wg_high_mapq_depth_list)

    for j in range(0, len(input_sv_list)):
        svcall = input_sv_list[j]
        if svcall.svtype == 'INV' or svcall.svtype == 'TRA' or svcall.chrm1 != svcall.chrm2:
            continue

        bk1_pos = svcall.start1
        bk2_pos = svcall.end2
        tid = chrname2tid_dict[svcall.chrm1]

        sv_region_total_depth_list = get_region_depth_list(wg_total_depth_list, tid, bk1_pos, bk2_pos, bin_size)
        n_dots = len(sv_region_total_depth_list)
        if n_dots < 2: continue

        mean_reg_total_depth = np.median(sv_region_total_depth_list)
        q1_reg_total_depth = np.percentile(sv_region_total_depth_list, 0.25)
        q3_reg_total_depth = np.percentile(sv_region_total_depth_list, 0.75)

        if svcall.svtype == 'DEL' and (mean_reg_total_depth <= wg_depth_mean * 0.667 and q3_reg_total_depth < wg_depth_mean):
            continue
        elif svcall.svtype == 'DUP' and (mean_reg_total_depth >= wg_depth_mean * 1.25 and q1_reg_total_depth > wg_depth_mean):
            continue
        else:
            input_sv_list[j].ft = 'DEPTH_FILTER'

    return input_sv_list


def get_region_depth_list(wg_depth_list, tid, start, end, bin_size):

    start_idx = int(start / bin_size)
    end_idx = int(end / bin_size) + 1
    return wg_depth_list[tid][start_idx:end_idx]

def get_mean_std_depth(wg_depth_list):

    linear_list = list()
    for tid in range(0, len(wg_depth_list)):
        linear_list += wg_depth_list[tid]

    wg_depth_mean = np.mean(linear_list)
    wg_depth_std  = np.std(linear_list)

    return wg_depth_mean, wg_depth_std

def filter_dbo_score(input_sv_list, args):

    min_inv_dbo_score = 1
    min_tra_dbo_score = 0.5
    if args.is_wgs and args.germline_mode:
        min_inv_dbo_score = 2
        min_tra_dbo_score = 1

    for i in range(0, len(input_sv_list)):
        svcall = input_sv_list[i]
        if svcall.svtype == 'INV' and min(svcall.dbo_score1, svcall.dbo_score2) < min_inv_dbo_score:
            input_sv_list[i].ft = 'LOW_DBO_SCORE'
        elif svcall.chrm1 != svcall.chrm2 and min(svcall.dbo_score1, svcall.dbo_score2) < min_tra_dbo_score:
            input_sv_list[i].ft = 'LOW_DBO_SCORE'

    return input_sv_list

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

def filter_calls_2d (svcall_list,  black_list_file, out_file, remove_chr_prefix = False):

    black_list_2array_dict, bin_size = read_2d_blacklist_file(black_list_file, remove_chr_prefix)

    mean_fragment_length = 20000
    box_length = mean_fragment_length

    for i in range(0, len(svcall_list)):
        svcall = svcall_list[i]
        if svcall.ft != '.':  continue

        chr1 = svcall.chrm1
        pos1 = svcall.start1
        chr2 = svcall.chrm2 
        pos2 = svcall.start2 
        end_type1 = svcall.endtype1
        end_type2 = svcall.endtype2
        
        if end_type1 == 'R_end':
            start1 = pos1 - box_length
        elif end_type1 == 'L_end':
            start1 = pos1 

        if end_type2 == 'R_end':
            start2 = pos2 - box_length
        elif end_type2 == 'L_end':
            start2 = pos2

        end1 = start1 + box_length
        end2 = start2 + box_length

        key1 = two_chr_to_key(chr1, chr2)
        key2 = two_chr_to_key(chr2, chr1)

        if key1 in black_list_2array_dict:
            pos1_list, pos2_list_list = black_list_2array_dict[key1]
            number_of_points = get_number_of_points_from_black_list_file(start1, end1, start2, end2, pos1_list, pos2_list_list, bin_size) 
        elif key2 in black_list_2array_dict:
            my_utils.myprint('switch chr1 and chr2')
            pos1_list, pos2_list_list = black_list_2array_dict[key2]
            number_of_points = get_number_of_points_from_black_list_file(start2, end2, start1, end1, pos1_list, pos2_list_list, bin_size) 
        else:
            number_of_points = 0
      
        if number_of_points >= 20 : svcall_list[i].ft = '2D_BLACKLIST'
        
    return svcall_list


def get_number_of_points_from_black_list_file(start1, end1, start2, end2, pos1_list, pos2_list_list, bin_size):

    start1_bin = start1 - start1 % bin_size 
    end1_bin = end1 - end1 % bin_size

    start2_bin = start2 - start2 % bin_size 
    end2_bin = end2 - end2 % bin_size


    number_of_points = 0
    idx1 = max(0, bisect.bisect_left(pos1_list, start1_bin) - 1)


    while idx1 < len(pos1_list) and pos1_list[idx1] <= end1_bin: 
        if pos1_list[idx1] >= start1_bin:
            pos2_list = pos2_list_list[idx1]
            idx2 = max(0, bisect.bisect_left(pos2_list, start2_bin) - 1)

            while idx2 < len(pos2_list) and pos2_list[idx2] <= end2_bin: 
                if pos2_list[idx2] >= start2_bin:
                    number_of_points += 1
                idx2 += 1

        idx1 += 1

    return number_of_points

def read_2d_blacklist_file(black_list_file, remove_chr_prefix):

    black_list_2d_dict = dict()

    black_list_fp = gzip.open(black_list_file, 'rt')

    bin_size = 0
    
    while 1:
        line = black_list_fp.readline()
        if not line: break
        if line[0] == '#': 
            line = line.strip().split('=')
            bin_size = int(line[1])
            my_utils.myprint ('bin_size = %d' % bin_size)
            continue

        if line[0] == '>':
            chr_list = line[1:].strip().split(',')
            chr1 = chr_list[0]
            chr2 = chr_list[1]

            if remove_chr_prefix:
                chr1 = chr1[3:]
                chr2 = chr2[3:]

            key = two_chr_to_key(chr1, chr2) 
            black_list_2d_dict[key] = list()
            continue

        line = line.strip().split(tab)
        pos1 = int(line[0])
        pos2_list = line[1].split(',')
        for i in range(0, len(pos2_list)):
            pos2_list[i] = int(pos2_list[i])

        d2 = D2(pos1, pos2_list)
        black_list_2d_dict[key].append(d2)
        
    black_list_fp.close() 
    
    black_list_2array_dict = dict()
    for key in black_list_2d_dict:
        black_list_2d_dict[key].sort(key = lambda d2: d2.pos1) 
        black_list_2array_dict[key] = ( list(), list() )
        for d2 in black_list_2d_dict[key]:
            black_list_2array_dict[key][0].append(d2.pos1)
            black_list_2array_dict[key][1].append(d2.pos2_list)

    del black_list_2d_dict

    return black_list_2array_dict, bin_size

def two_chr_to_key(chr1, chr2):

    return '%s\t%s' % (chr1, chr2)

def filter_1d_blacklist(raw_svcall_list, black_reg_dict, alt_chr_name_set, gap_left_region_dict, gap_right_region_dict, bin_size ):

    for i in range(0, len(raw_svcall_list)):
        svcall = raw_svcall_list[i]
        if svcall.ft != '.': continue
        if svcall.chrm1 != svcall.chrm2: 
            sv_length = int(1e10)
        else:
            sv_length = abs(svcall.start2 - svcall.start1)

        if svcall.chrm1 in alt_chr_name_set or svcall.chrm2 in alt_chr_name_set: 
            raw_svcall_list[i].ft = 'NON_PRIMARY_CONTIG'
            continue

        index1 = int(svcall.start1 / bin_size)
        index2 = int(svcall.start2 / bin_size)

        if svcall.chrm1 in black_reg_dict and index1 in black_reg_dict[svcall.chrm1]: 
            raw_svcall_list[i].fl = '1D_BLACKLIST'
            continue

        if svcall.chrm2 in black_reg_dict and index2 in black_reg_dict[svcall.chrm2]: 
            raw_svcall_list[i].fl = '1D_BLACKLIST'
            continue

        if sv_length < 200000 and svcall.endtype1 == 'R_end' and svcall.chrm1 in gap_left_region_dict and index1 in gap_left_region_dict[svcall.chrm1]: 
            raw_svcall_list[i].fl = 'GAP_REGION'

        if sv_length < 200000 and svcall.endtype1 == 'L_end' and svcall.chrm1 in gap_right_region_dict and index1 in gap_right_region_dict[svcall.chrm1]:
            raw_svcall_list[i].fl = 'GAP_REGION'

        if sv_length < 200000 and svcall.endtype2 == 'R_end' and svcall.chrm2 in gap_left_region_dict and index2 in gap_left_region_dict[svcall.chrm2]: 
            raw_svcall_list[i].fl = 'GAP_REGION'

        if sv_length < 200000 and svcall.endtype2 == 'L_end' and svcall.chrm2 in gap_right_region_dict and index2 in gap_right_region_dict[svcall.chrm2]:
            raw_svcall_list[i].fl = 'GAP_REGION'

    return raw_svcall_list


def filter_low_mapq_gaps(input_sv_list, endpoint_args, chrname2tid_dict):

    all_supp_barcode_dict = dict()

    for j in range(0, len(input_sv_list)):
        svcall = input_sv_list[j]
        if svcall.ft != '.': continue

        support_barcode_list = svcall.support_barcodes.rstrip(',').split(',')
        for bcd in support_barcode_list:
            all_supp_barcode_dict[bcd] = list()

    if os.path.exists(endpoint_args.low_mapq_bcd21_file) == False:
        my_utils.myprint('WARNING! low mapq bcd21 file does not exist.')
        return input_sv_list

    my_utils.myprint('reading low mapq bcd21 file: %s'  % endpoint_args.low_mapq_bcd21_file)
    low_mapq_bcd21_fp = my_utils.gzopen(endpoint_args.low_mapq_bcd21_file, 'rt')
    i = 0 
    while 1:
        line = low_mapq_bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        i += 1

        attr_list = line.strip().split(tab)
        bcd21 = Bcd21Core(attr_list)
        if bcd21.bcd in all_supp_barcode_dict:
            all_supp_barcode_dict[bcd21.bcd].append(bcd21)
        if i % 10000000 == 0:
            my_utils.myprint('processed %d reads' % i)

    low_mapq_bcd21_fp.close()

    for bcd in all_supp_barcode_dict:
        all_supp_barcode_dict[bcd].sort(key = lambda bcd21: bcd21.key_start() )

    my_utils.myprint('finished reading low mapq bcd21 file: %s'  % endpoint_args.low_mapq_bcd21_file)

    region_size = 10 * 1000
    # for deletion, the region is the deletion region, for other type of svs, the region is 10 kb of either breakpoint

    for j in range(0, len(input_sv_list)):
        svcall = input_sv_list[j]
        if svcall.ft != '.': continue

        support_barcode_list = svcall.support_barcodes.rstrip(',').split(',')
        n_low_mapq_bcd = 0
        region_key_start1   = -1
        region_key_end1     = -1
        region_key_start2   = -1
        region_key_end2     = -1

        tid1 = chrname2tid_dict[svcall.chrm1]
        tid2 = chrname2tid_dict[svcall.chrm2]

        if svcall.endtype1 == 'L_end':
            region_key_start1 = tid1 * my_utils.FIX_LENGTH + svcall.start1 - region_size
        else:
            region_key_start1 = tid1 * my_utils.FIX_LENGTH + svcall.start1 

        if svcall.endtype2 == 'L_end':
            region_key_start2 = tid2 * my_utils.FIX_LENGTH + svcall.start2 - region_size
        else:
            region_key_start2 = tid2 * my_utils.FIX_LENGTH + svcall.start2 

        region_key_end1 = region_key_start1 + region_size
        region_key_end2 = region_key_start2 + region_size

        if svcall.svtype == 'DEL':
            region_key_start1 = tid1 * my_utils.FIX_LENGTH + svcall.start1
            region_key_end1   = tid2 * my_utils.FIX_LENGTH + svcall.start2
            region_key_start2 = region_key_start1
            region_key_end2 = region_key_end1

        n_low_mapq_bcd = 0

        for bcd in support_barcode_list:
            if bcd not in all_supp_barcode_dict: continue
            for bcd21 in all_supp_barcode_dict[bcd]:
                if (bcd21.key_start() > region_key_start1 and bcd21.key_end() < region_key_end1) or (bcd21.key_start() > region_key_start2 and bcd21.key_end() < region_key_end2) : 
                    n_low_mapq_bcd += 1
                    break

        n_supp_bcd = svcall.num_fragment_support 
        ratio_low_mapq_bcd = float(n_low_mapq_bcd) / float(n_supp_bcd)
        if (not (ratio_low_mapq_bcd < 0.2 and svcall.score * (1-ratio_low_mapq_bcd) > 20)):
            input_sv_list[j].ft = 'LOW_MAPQ_BETWEEN_BK'

    return input_sv_list





if __name__ == '__main__':
    main()
