#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from my_utils import *
from bed import *

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + ' <in.depth.txt> <in.bed> <faidx_file> <out_dir> <type_of_depth (read or bcd)> '
argc  = 5 

def main():
    if len(arg) < argc:
        print (usage)
        sys.exit()

    in_depth_file         = arg.pop(0)
    in_bed_file           = arg.pop(0)
    faidx_file            = arg.pop(0)
    out_dir               = arg.pop(0)
    depth_type            = arg.pop(0)


    depth_type = depth_type.lower()
    if depth_type == 'barcode': depth_type = 'bcd'

    if depth_type != 'read' and depth_type != 'bcd':
        print('ERROR! type_of_depth can only be \'read\' or \'bcd\'')
        print(usage)
        sys.exit()
        
    if os.path.exists(out_dir) == False:
        os.makedirs(out_dir)

    chr_len_list = get_chr_length(faidx_file)
    tid2chrname_list, chrname2tid_dict = get_chrnames(faidx_file)
    wg_high_mapq_depth_list, wg_total_depth_list, bin_size = get_wg_depth_list(in_depth_file, chr_len_list)
    
    wg_total_depth = 0
    wg_n_bin = 0
    for tid in range(0, len(wg_high_mapq_depth_list)):
        for depth in wg_high_mapq_depth_list[tid]:
            wg_total_depth += depth
            wg_n_bin += 1

    wg_avg_depth = float(wg_total_depth) / wg_n_bin 
    input_bed_list = read_bedsvcall_file(in_bed_file, chrname2tid_dict)

    for bed in input_bed_list:
        plot_read_depth_for1bed(bed, out_dir, wg_high_mapq_depth_list, wg_total_depth_list, chr_len_list, bin_size, depth_type, wg_avg_depth)

    return

def plot_read_depth_for1bed(bed, out_dir, wg_high_mapq_depth_list, wg_total_depth_list, chr_len_list, bin_size, depth_type, wg_avg_depth):
    
    plt.figure(figsize=(10, 5)) 
    tid = bed.tid
    if depth_type == 'read':
        out_figure_name = os.path.join(out_dir, '%s.%s_%d_%d.read_depth.png' % (bed.svtype, bed.chrm, bed.start, bed.end))
    else:
        out_figure_name = os.path.join(out_dir, '%s.%s_%d_%d.barcode_depth.png' % (bed.svtype, bed.chrm, bed.start, bed.end))

    bed_len = bed.end - bed.start
    win_start = max(0, bed.start - bed_len)
    win_end = min(chr_len_list[tid], bed.end + bed_len)
    
    win_start_idx = int(win_start/bin_size)
    win_end_idx = int(win_end/bin_size) + 1
   
    x = range(win_start_idx * bin_size, win_end_idx * bin_size, bin_size)
    y1 = wg_high_mapq_depth_list[tid][win_start_idx:win_end_idx]
    y2 = wg_total_depth_list[tid][win_start_idx:win_end_idx]

    ymean = np.mean(y2)

    ymax = ymean * 2

    if ymax < wg_avg_depth * 2:
        ymax = wg_avg_depth * 2

    if depth_type == 'read':
        plt.title('%s:%d-%d %s Read Depth' % (bed.chrm, bed.start, bed.end, bed.svtype))
    else:
        plt.title('%s:%d-%d %s Barcode Depth' % (bed.chrm, bed.start, bed.end, bed.svtype))

    plt.xlabel('Position')
    if depth_type == 'read':
        plt.ylabel('Read depth')
    else:
        plt.ylabel('Barcode per bp')

    plt.plot(x, y2, '-', color = 'grey')
    plt.plot(x, y1, '-', color = 'black')

    plt.axis([win_start, win_end, 0, ymax])
    plt.axvline(x=bed.start, color='r', linestyle = '--')
    plt.axvline(x=bed.end, color='r', linestyle = '--')
    plt.axhline(y=wg_avg_depth, color='b', linestyle = '--')
    plt.ticklabel_format(axis='both', style='plain')
    plt.xticks(np.arange(min(x), max(x)+1, bed_len))
    plt.rcParams.update({'font.size': 16})
    plt.show()
    plt.savefig(out_figure_name)
    plt.close('all')
    return
 

def get_wg_depth_list (in_depth_file, chr_len_list):

    if (in_depth_file[-2:] == 'gz'):
        in_depth_fp = gzip.open(in_depth_file, 'r') 
    else:
        in_depth_fp = open(in_depth_file, 'r') 

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
            myprint('ERROR! bin_size < 1 in depth file: %s ' % in_depth_file)
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

    myprint('finished reading file: %s' % in_depth_file)
    return wg_high_mapq_depth_list, wg_total_depth_list, bin_size

if __name__ == '__main__':
    main()
