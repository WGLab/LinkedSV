#!/usr/bin/env python

import numpy as np
from my_utils import *
from fragment import *

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()
    cal_endpoint_density (args, dbo_args, endpoint_args)

    return

def cal_endpoint_density (args, dbo_args, endpoint_args):

    out_prefix = args.out_prefix
    faidx_file = args.faidx_file

    bcd22_file = endpoint_args.bcd22_file
    out1_file = endpoint_args.endpoints_density_file 

    bin_len = endpoint_args.bin_len_for_calculate_endpoint_density
    total_win_length = args.gap_distance950 

    if total_win_length > args.median_fragment_length * 0.2:
        total_win_length = args.median_fragment_length * 0.2
    
    win_len = int(round(float(total_win_length)/bin_len))

    myprint('window length for calculating enriched endpoints: %d' % (win_len * bin_len))

    in_fp     = open(bcd22_file, 'r')
    fai_fp    = open(faidx_file, 'r')
    out1_fp   = open(out1_file, 'w')
    myprint('reading fai file: ' + faidx_file)
    chr_len = list()
    while 1:
        line = fai_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        chr_len.append(int(line[1]))
    fai_fp.close()

    myprint('init depth list')
    frag_depth_list = range(0, len(chr_len))
    frag_start_list = range(0, len(chr_len))
    frag_end_list   = range(0, len(chr_len))
    for i in range(0, len(chr_len)):
        frag_depth_list[i] = [0] * (int(chr_len[i]/bin_len)+1)
        frag_start_list[i] = [0] * (int(chr_len[i]/bin_len)+1)
        frag_end_list[i]   = [0] * (int(chr_len[i]/bin_len)+1)

    frag_window_depth = range(0, len(chr_len))
    frag_window_start = range(0, len(chr_len))
    frag_window_end   = range(0, len(chr_len))
    start_density     = range(0, len(chr_len))
    end_density       = range(0, len(chr_len))

    for i in range(0, len(chr_len)):
        frag_window_depth[i] = [0] * (int(chr_len[i]/bin_len)+1)
        frag_window_start[i] = [0] * (int(chr_len[i]/bin_len)+1)
        frag_window_end[i]   = [0] * (int(chr_len[i]/bin_len)+1)
        start_density[i]     = [0] * (int(chr_len[i]/bin_len)+1)
        end_density[i]       = [0] * (int(chr_len[i]/bin_len)+1)


    myprint('reading bcd22 file: ' + bcd22_file)
    line_c = 0
    while 1:
        line = in_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)

        bcd22_frm = Fragment(line) 
        frag_tid = bcd22_frm.tid
        frag_len = bcd22_frm.length
        if frag_len < endpoint_args.min_frag_length: continue
        frag_start = bcd22_frm.start
        frag_end = bcd22_frm.end

        bin_start = int(frag_start / bin_len)
        bin_end = int(frag_end / bin_len)

        frag_start_list[frag_tid][bin_start] += 1
        frag_end_list[frag_tid][bin_end] += 1
        for j in range(bin_start, bin_end+1):
            frag_depth_list[frag_tid][j] += 1

        line_c += 1
        if line_c % 1e6 == 0:
            myprint('processed ' + str(line_c) + ' fragments')

    myprint('outputing results into file: ' + out1_file)
    start_density_list = list()
    end_density_list = list()
    for i in range(0, len(chr_len)):
        j = 0;
        frag_window_depth[i][j] = 0
        frag_window_start[i][j] = 0
        frag_window_end[i][j] = 0
        for k in range(0, win_len):
            frag_window_depth[i][j] += frag_depth_list[i][j+k]
            frag_window_start[i][j] += frag_start_list[i][j+k]
            frag_window_end[i][j]   += frag_end_list[i][j+k]

        if frag_window_depth[i][j] > 0:
            start_density[i][j] = float(frag_window_start[i][j]) / float(frag_window_depth[i][j])
            end_density[i][j]   = float(frag_window_end[i][j]) / float(frag_window_depth[i][j])
            start_density_list.append(start_density[i][j])
            end_density_list.append(end_density[i][j])
        else:
            start_density[i][j] = 0
            end_density[i][j] = 0

        out1_fp.write(str(i) + tab + str(j*bin_len) + tab + str((j+win_len)*bin_len-1) + tab + str(frag_window_depth[i][j]) + tab + str(frag_window_start[i][j]) + tab + str(frag_window_end[i][j]) + tab + str(start_density[i][j]) + tab + str(end_density[i][j]) + endl)

        for j in range(1, int(chr_len[i]/bin_len)-win_len+1):
            frag_window_depth[i][j] = frag_window_depth[i][j-1] - frag_depth_list[i][j-1] + frag_depth_list[i][j+win_len-1] 
            frag_window_start[i][j] = frag_window_start[i][j-1] - frag_start_list[i][j-1] + frag_start_list[i][j+win_len-1]
            frag_window_end[i][j]   = frag_window_end[i][j-1]   - frag_end_list[i][j-1]   + frag_end_list[i][j+win_len-1]
            if frag_window_depth[i][j] > 0:
                start_density[i][j] = float(frag_window_start[i][j]) / float(frag_window_depth[i][j])
                end_density[i][j]   = float(frag_window_end[i][j]) / float(frag_window_depth[i][j])
                start_density_list.append(start_density[i][j])
                end_density_list.append(end_density[i][j])
            else:
                start_density[i][j] = 0
                end_density[i][j] = 0
            out1_fp.write(str(i) + tab + str(j*bin_len) + tab + str((j+win_len)*bin_len-1) + tab + str(frag_window_depth[i][j]) + tab + str(frag_window_start[i][j]) + tab + str(frag_window_end[i][j]) + tab + str(start_density[i][j]) + tab + str(end_density[i][j]) + endl)
    
    out1_fp.close()

    myprint('finished calculating endpoint density')

    return 

if __name__ == '__main__':
    main()
