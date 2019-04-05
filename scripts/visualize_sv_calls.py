#!/usr/bin/env python

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import my_utils
import bedpe
import bed
import plot_read_depth
import plot_2d_barcodes


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in.svcalls.bedpe> <bcd13_file> <bcd21_file> <faidx_file> <cal_2d_overlapping_barcodes_binary> <cal_read_depth_from_bcd21> <out_dir> <out_prefix>'
argc  = 8 


def main():
    if len(arg) < argc:
        print (usage)
        sys.exit()

    in_svcalls_bedpe_file = arg.pop(0)
    bcd13_file = arg.pop(0)
    bcd21_file = arg.pop(0)
    faidx_file = arg.pop(0)
    cal_2d_overlapping_barcodes_binary = arg.pop(0)
    cal_read_depth_from_bcd21_binary = arg.pop(0)
    out_dir = arg.pop(0)
    out_prefix = arg.pop(0)

    visualize_sv_calls(in_svcalls_bedpe_file, bcd13_file, bcd21_file, faidx_file, cal_2d_overlapping_barcodes_binary, cal_read_depth_from_bcd21_binary, out_dir, out_prefix)

    return

def visualize_sv_calls(in_svcalls_bedpe_file, bcd13_file, bcd21_file, faidx_file, cal_2d_overlapping_barcodes_binary, cal_read_depth_from_bcd21_binary, out_dir, out_prefix):

    out_dir = os.path.abspath(out_dir)
    my_utils.make_dir(out_dir)

    cal_2d_overlapping_barcodes_binary = os.path.abspath(cal_2d_overlapping_barcodes_binary)
    cal_read_depth_from_bcd21_binary = os.path.abspath(cal_read_depth_from_bcd21_binary)

    chr_len_list = my_utils.get_chr_length(faidx_file)
    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(faidx_file)
    in_svcall_list = bedpe.read_svcall_bedpe_file(in_svcalls_bedpe_file, chrname2tid_dict)
    for i in range(0, len(in_svcall_list)):
        in_svcall_list[i].format() 

    plot_depth(cal_read_depth_from_bcd21_binary, bcd21_file, in_svcall_list, faidx_file, chr_len_list, tid2chrname_list, chrname2tid_dict, out_dir, out_prefix)

    plot_twin_window_barcode_similarity(in_svcall_list, bcd13_file, faidx_file, out_dir, chr_len_list, tid2chrname_list, chrname2tid_dict, out_prefix)

    flank_dist = 100 * 1000 # set flank distance to be 100 kb
    plot_heatmap (in_svcall_list, bcd21_file, faidx_file, out_dir, flank_dist, chr_len_list, tid2chrname_list, chrname2tid_dict, cal_2d_overlapping_barcodes_binary, out_prefix)

    return


def generate_target_region_bedpe_list(in_svcall_list, chr_len_list, flank_dist, chrname2tid_dict):

    target_region_bedpe_list = list()
    for i in range(0, len(in_svcall_list)):
        svcall = in_svcall_list[i]
        tid1 = svcall.tid1
        start1 = svcall.start1 - flank_dist
        end1 = svcall.start1 + flank_dist
        tid2 = svcall.tid2
        start2 = svcall.start2 - flank_dist
        end2 = svcall.start2 + flank_dist

        if start1 < 0: start1 = 0 
        if start2 < 0: start2 = 0

        if end1 > chr_len_list[tid1]: end1 = chr_len_list[tid1]
        if end2 > chr_len_list[tid2]: end2 = chr_len_list[tid2]

        bedpe1 = bedpe.BedpeSVCall(svcall.chrm1, start1, end1, svcall.chrm2, start2, end2, svcall.svtype, svcall.sv_id, chrname2tid_dict)
        target_region_bedpe_list.append(bedpe1)

    return target_region_bedpe_list


def plot_heatmap(in_svcall_list, bcd21_file, faidx_file, out_dir, flank_dist, chr_len_list, tid2chrname_list, chrname2tid_dict, cal_2d_overlapping_barcodes_binary, out_prefix):

    if os.path.exists(cal_2d_overlapping_barcodes_binary) == False:
        my_utils.myprint('ERROR! The binary file doesn\'t exist: %s' % cal_2d_overlapping_barcodes_binary)
        my_utils.myprint('Skipped plotting the heat maps')
        return

    out_dir = os.path.join(out_dir, '2D_heatmap')
    my_utils.make_dir(out_dir)

    my_utils.myprint('plotting heat maps of overlapping barcodes')

    target_region_bedpe_list = generate_target_region_bedpe_list(in_svcall_list, chr_len_list, flank_dist, chrname2tid_dict)
    target_region_bedpe_file = os.path.join(out_dir, 'target_region.bedpe') 
    target_region_bedpe_fp = my_utils.gzopen(target_region_bedpe_file, 'w')  
    for bedpe1 in target_region_bedpe_list: 
        target_region_bedpe_fp.write(bedpe1.output_svcall() + endl)
    target_region_bedpe_fp.close()

    target_region_2d_ovl_with_low_mapq_file = os.path.join(out_dir, '%s.2d_heatmap.with_low_mapq_reads.txt' % out_prefix)

    bin_size = 1000
    max_ovl_num = 100

    cmd_args_list1 = [cal_2d_overlapping_barcodes_binary, bcd21_file, target_region_bedpe_file, target_region_2d_ovl_with_low_mapq_file, faidx_file, str(bin_size), '1']
    subprocess.call(cmd_args_list1)
    plot_2d_barcodes.plot_2d_overlapping_barcodes(target_region_2d_ovl_with_low_mapq_file, target_region_bedpe_list, bin_size, max_ovl_num, out_dir, out_prefix)

    return

def plot_depth(cal_read_depth_from_bcd21_binary, bcd21_file, in_svcalls_list, faidx_file, chr_len_list, tid2chrname_list, chrname2tid_dict, out_dir, out_prefix):

    if os.path.exists(cal_read_depth_from_bcd21_binary) == False:
        my_utils.myprint('ERROR! The binary file doesn\'t exist:%s\nSkipped plotting read depth' % cal_read_depth_from_bcd21_binary)
        return

    out_dir = os.path.join(out_dir, 'read_depth')
    my_utils.make_dir(out_dir)

    bin_size = 1000
    read_depth_file = os.path.join(out_dir, '%s.read_depth.txt' % out_prefix)
    cmd_args_list = [cal_read_depth_from_bcd21_binary, bcd21_file, read_depth_file, faidx_file, str(bin_size), '20']
    my_utils.myprint('calculating read depth from file: %s' % bcd21_file)
    subprocess.call(cmd_args_list)

    my_utils.myprint('plotting read depth')
    wg_high_mapq_depth_list, wg_total_depth_list, bin_size = plot_read_depth.get_wg_depth_list(read_depth_file, chr_len_list)
    
    wg_total_depth = 0 
    wg_n_bin = 0 
    for tid in range(0, len(wg_high_mapq_depth_list)):
        for depth in wg_high_mapq_depth_list[tid]:
            wg_total_depth += depth
            wg_n_bin += 1

    wg_avg_depth = float(wg_total_depth) / wg_n_bin 

    for svcall in in_svcalls_list:
        if svcall.chrm1 != svcall.chrm2: continue
        out_file = os.path.join(out_dir, '%s.%s.read_depth.png' % (out_prefix, svcall.sv_id))
        figure_title = 'Read depth (%s, %d bp %s)' % (svcall.sv_id, svcall.end2 - svcall.start1, svcall.svtype)
        plot_read_depth.plot_read_depth_for1region(svcall.chrm1, svcall.tid1, svcall.start1, svcall.end2, out_file, figure_title, wg_high_mapq_depth_list, wg_total_depth_list, chr_len_list, bin_size, wg_avg_depth)

    return

def plot_twin_window_barcode_similarity(in_svcall_list, bcd13_file, faidx_file, out_dir, chr_len_list, tid2chrname_list, chrname2tid_dict, out_prefix):

    out_dir = os.path.join(out_dir, 'twin_window_barcode_similarity')
    my_utils.make_dir(out_dir)

    wg_pvalue_list, bin_size = get_wg_pvalue_list_from_bcd13_file(bcd13_file, chr_len_list)

    max_length_in_one_figure = 500 * 1000
    for svcall in in_svcall_list:
        if svcall.svtype == 'DEL' or svcall.svtype == 'DUP': continue
        flank_dist = abs(svcall.key2 - svcall.key1)
        if flank_dist > 50 * 1000: flank_dist = 50 * 1000
        if flank_dist < 10 * 1000: flank_dist = 10 * 1000
        if abs(svcall.key2 - svcall.key1) > max_length_in_one_figure:
            out_file1 = os.path.join(out_dir, '%s.%s.breakpoint1.twin_window_barcode_similarity.png' % (out_prefix, svcall.sv_id )) 
            reg_start = svcall.start1 - flank_dist 
            reg_end = svcall.end1 + flank_dist 
            if reg_start < 0: reg_start = 0
            if reg_end > chr_len_list[svcall.tid1]: reg_end = chr_len_list[svcall.tid1]
            bk_pos1 = svcall.start1
            bk_pos2 = bk_pos1 
            figure_title = '%s, %s (breakpoint 1)' % (svcall.sv_id, svcall.svtype)
            plot_twin_window_barcode_similarity_for1region(svcall.chrm1, svcall.tid1, reg_start, reg_end, bk_pos1, bk_pos2, out_file1, figure_title, wg_pvalue_list, bin_size)

            out_file2 = os.path.join(out_dir, '%s.%s.breakpoint2.twin_window_barcode_similarity.png' % (out_prefix, svcall.sv_id )) 
            reg_start = svcall.start2 - flank_dist 
            reg_end = svcall.end2 + flank_dist 
            bk_pos1 = svcall.start2
            bk_pos2 = bk_pos1 
            if reg_start < 0: reg_start = 0
            if reg_end > chr_len_list[svcall.tid2]: reg_end = chr_len_list[svcall.tid2]
            figure_title = '%s, %s (breakpoint 2)' % (svcall.sv_id, svcall.svtype)
            plot_twin_window_barcode_similarity_for1region(svcall.chrm2, svcall.tid2, reg_start, reg_end, bk_pos1, bk_pos2, out_file2, figure_title, wg_pvalue_list, bin_size) 
        else:
            out_file  = os.path.join(out_dir, '%s.%s.both_breakpoints.twin_window_barcode_similarity.png' % (out_prefix, svcall.sv_id)) 
            reg_start = svcall.start1 - flank_dist 
            reg_end = svcall.end2 + flank_dist 
            if reg_start < 0: reg_start = 0
            if reg_end > chr_len_list[svcall.tid1]: reg_end = chr_len_list[svcall.tid1]
            bk_pos1 = svcall.start1
            bk_pos2 = svcall.start2
            figure_title = '%s, %s (breakpoint 1 and 2)' % (svcall.sv_id, svcall.svtype)
            plot_twin_window_barcode_similarity_for1region(svcall.chrm1, svcall.tid1, reg_start, reg_end, bk_pos1, bk_pos2, out_file, figure_title, wg_pvalue_list, bin_size) 

    return

def plot_twin_window_barcode_similarity_for1region(chrom, tid, reg_start, reg_end, bk_pos1, bk_pos2, out_file, figure_title, wg_pvalue_list, bin_size):

    min_x_idx = int(reg_start / bin_size)
    max_x_idx = int(reg_end / bin_size) + 1

    x_list = list()
    y_list = list()

    for idx in range(min_x_idx, max_x_idx):
        x = idx * bin_size
        x_list.append(x)
        if idx < len(wg_pvalue_list[tid]):
            y = wg_pvalue_list[tid][idx]
            y_list.append(y)
        else:
            break

    xmin = min(x_list)
    xmax = max(x_list)
    x_range = xmax - xmin
    ymax = max(y_list)
    

    plt.figure(figsize=(10, 5))
    plt.title(figure_title)
    plt.xlabel('%s position' % chrom)
    plt.ylabel('-log10(P-value)')
    plt.plot(x_list, y_list, '-', color = 'black')
    plt.axis([xmin, xmax, 0, ymax])
    plt.axvline(x=bk_pos1, color='r', linestyle = '--')
    plt.axvline(x=bk_pos2, color='r', linestyle = '--')
    plt.xticks(np.arange(xmin, xmax+1, x_range/4))
    plt.show()
    plt.savefig(out_file, dpi=200)
    plt.close('all')
    my_utils.myprint('saved figure: %s' % out_file)

    return

def get_wg_pvalue_list_from_bcd13_file(bcd13_file, chr_len_list):

    bin_size = get_bin_size_from_bcd13_file(bcd13_file)
    n_chr = len(chr_len_list)
    wg_pvalue_list = [0] * n_chr
    for i in range(0, n_chr):
        chr_len = chr_len_list[i]
        n_bin = int(chr_len / bin_size + 2)
        wg_pvalue_list[i] = [0] * n_bin

    bcd13_fp = my_utils.gzopen(bcd13_file, 'r')
    while 1:
        line = bcd13_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        tid = int(line[0])
        pos = int(line[1])
        idx = int(pos / bin_size)
        pvalue = float(line[8])
        if idx < len(wg_pvalue_list[tid]):
            wg_pvalue_list[tid][idx] = pvalue

    bcd13_fp.close()

    return wg_pvalue_list, bin_size

def get_bin_size_from_bcd13_file(bcd13_file):

    bcd13_fp = my_utils.gzopen(bcd13_file, 'r')
    pos_list = list()

    while 1:
        line = bcd13_fp.readline()
        if not line: break
        if line[0] == '#': continue

        line = line.strip().split(tab)
        pos = int(line[1])
        pos_list.append(pos)
        if len(pos_list) > 1000000: break
    bcd13_fp.close()
    
    interval_list = list()
    for i in range(1, len(pos_list)):
        interval = pos_list[i] - pos_list[i-1]
        if interval > 0: interval_list.append(interval)

    if len(interval_list) < 1:
        my_utils.myprint('Failed to get bin size from file: %s' % bcd13_file)
        sys.exit(1)

    bin_size = int(np.median(interval_list))

    del pos_list
    del interval_list

    return bin_size
            








if __name__ == '__main__':
    main()
