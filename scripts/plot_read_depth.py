#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from scripts import my_utils
except ImportError:
    import my_utils

tab  = '\t'
endl = '\n'


def plot_read_depth_for1region(chrom, tid, bk_pos1, bk_pos2, out_file, figure_title, wg_high_mapq_depth_list, wg_total_depth_list, chr_len_list, bin_size, wg_avg_depth):
    
    plt.figure(figsize=(10, 5)) 

    if bk_pos2 < bk_pos1:
        temp = bk_pos2
        bk_pos2 = bk_pos1
        bk_pos1 = temp

    sv_len = bk_pos2 - bk_pos1 

    win_start = max(0, bk_pos1 - sv_len)
    win_end = min(chr_len_list[tid], bk_pos2 + sv_len)
    
    win_start_idx = int(win_start/bin_size)
    win_end_idx = int(win_end/bin_size) + 1
   
    x = range(win_start_idx * bin_size, win_end_idx * bin_size, bin_size)
    y1 = wg_high_mapq_depth_list[tid][win_start_idx:win_end_idx]
    y2 = wg_total_depth_list[tid][win_start_idx:win_end_idx]

    ymean = np.mean(y2)

    ymax = ymean * 3

    if ymax < wg_avg_depth * 2:
        ymax = wg_avg_depth * 2

    plt.title(figure_title)

    plt.xlabel('%s position' % chrom)
    plt.ylabel('Read depth')

    plt.plot(x, y2, '-', color = 'grey')
    plt.plot(x, y1, '-', color = 'black')

    plt.axis([win_start, win_end, 0, ymax])
    plt.axvline(x=bk_pos1, color='r', linestyle = '--')
    plt.axvline(x=bk_pos2, color='r', linestyle = '--')
    plt.axhline(y=wg_avg_depth, color='b', linestyle = '--')
    plt.ticklabel_format(axis='both', style='plain')
    plt.xticks(np.arange(min(x), max(x)+1, sv_len))
    plt.rcParams.update({'font.size': 16})
    plt.show()
    plt.savefig(out_file, dpi=200)
    plt.close('all')
    my_utils.myprint('saved figure: %s' % out_file)

    return
 

def get_wg_depth_list (in_depth_file, chr_len_list):

    in_depth_fp = my_utils.gzopen(in_depth_file, 'r') 

    bin_size = 1
    while 1:
        line = in_depth_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        start_pos = int(line[1])
        end_pos   = int(line[2])
        bin_size = end_pos - start_pos 
        if bin_size >= 10:
            break
        else:
            my_utils.myprint('ERROR! bin_size < 1 in depth file: %s ' % in_depth_file)
            sys.exit()


    n_chr = len(chr_len_list)
    wg_high_mapq_depth_list = [0] * n_chr 
    wg_total_depth_list     = [0] * n_chr 
    for tid in range(0, n_chr):
        wg_high_mapq_depth_list[tid] = list()
        wg_total_depth_list[tid] = list()

    in_depth_fp.seek(0, 0)
    while 1:
        line = in_depth_fp.readline()
        if not line: break
        if line[0] == '#': continue

        line = line.strip().split(tab) 
        tid = int(line[0])
        wg_high_mapq_depth_list[tid].append(float(line[3]))
        wg_total_depth_list[tid].append(float(line[4]))

    in_depth_fp.close()

    my_utils.myprint('finished reading file: %s' % in_depth_file)
    return wg_high_mapq_depth_list, wg_total_depth_list, bin_size

if __name__ == '__main__':
    main()
