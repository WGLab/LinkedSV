#!/usr/bin/env python

from __future__ import division
from my_utils import *
from collections import deque

def main():
    
    args, dbo_args, endpoint_args = parse_user_arguments()
    count_barcode_overlap(args, dbo_args, endpoint_args)

    return

def count_barcode_overlap(args, dbo_args, endpoint_args):
    
    skip_seperate_bcd_chr = True
    for chr_bcd_file in dbo_args.chr_bcd_file_list:
        if os.path.exists(chr_bcd_file) == False:
            skip_seperate_bcd_chr = False
            break

    if args.run_from_begining == False and skip_seperate_bcd_chr:
        myprint('bcd chr files existed, skipped seperating bcd chr files')
    else:
        myprint('splitting bcd file')
        seperate_bcd_chr(args.copy(), dbo_args.copy())
        myprint('finished splitting bcd file')

    for tid in range(0, len(dbo_args.chr_bcd_file_list)):
        count_barcode_overlap_for_1_chr(tid, args.copy(), dbo_args.copy())

    return

def seperate_bcd_chr(args, dbo_args):

    # in_file: bcd file of all chrs)
    in_file = args.bcd_file
    chr_bcd_file_list = dbo_args.chr_bcd_file_list

    for tid in range(0, len(chr_bcd_file_list)):
        out_file = chr_bcd_file_list[tid]
        out_fp = open(out_file, 'w')
        out_fp.close()

    in_fp = open(in_file, 'r')
    prev_tid = -10
    file_idx = 0
    out_file = 0
    out_fp = ''
    while 1:
        line = in_fp.readline()
        if not line:
            break
        line2 = line.strip().split(tab)
        tid = int(line2[0])
        if tid != prev_tid:
            if out_file: 
                out_fp.close()
            if tid >= len(chr_bcd_file_list) or tid < 0:
                continue
                
            out_file = chr_bcd_file_list[tid] 
            out_fp = open(out_file, 'w')
            myprint('current chr_bcd file: %s' %out_file)
            prev_tid = tid

        out_fp.write(line)

    out_fp.close()
    in_fp.close()    
    return 

def cal_weight_center(rd_cnt_dq, window_side): 
    total_rd_cnt = 0
    total_weight_dist = 0
    max_dist = len(rd_cnt_dq)
    i = 0
    for rd_cnt in rd_cnt_dq: 
        total_weight_dist += (max_dist-i) * rd_cnt 
        total_rd_cnt += rd_cnt 
        i += 1

    if total_rd_cnt == 0:
        center1 = max_dist/2
    else:
        center1 = float(total_weight_dist)/total_rd_cnt 
    
    if window_side == 1:
        return center1
    else:
        return float(max_dist)-center1

def merge_window_all_bcdset(window_bcdset_dq):

    window_all_bcd_set = set()
    for bcdset in window_bcdset_dq:
        window_all_bcd_set = window_all_bcd_set.union(bcdset)

    return window_all_bcd_set

def count_barcode_overlap_for_1_chr(tid, args, dbo_args):

    myprint('counting barcode overlap for tid: ' + str(tid))
    min_map_qual = args.min_mapq
    out_prefix = args.out_prefix 
    bin_size = dbo_args.bin_size
    window_size = dbo_args.win_size
    chr_bcd_file_list = dbo_args.chr_bcd_file_list
    bcd11_file_list = dbo_args.bcd11_file_list 
    bcd12_file_list = dbo_args.bcd12_file_list 


    myprint('generating bcd11 file')
    bcd_chr_file = chr_bcd_file_list[tid]
    in_fp = open(bcd_chr_file, 'r')

    out_file = bcd11_file_list[tid] 
    out_fp = open(out_file, 'w')

    n_bin = 0
    barcode_dict = dict()
    prev_bin_id = -1
    while 1:
        line = in_fp.readline()
        if not line:
            break

        line = line.strip().split(tab)
        pos = int(line[1]) 
        mapq = int(line[3])
        if mapq < min_map_qual:
            continue

        barcode = line[4]
        hap_type = int(line[5])
        bin_id = int(pos/bin_size)

        if (bin_id > prev_bin_id):
            if (bin_id > prev_bin_id + 1):
                for i in range(prev_bin_id + 1, bin_id):
                    out_fp.write(str(i) + '\n')
            out_fp.write(str(bin_id) + '\t')
            for item in barcode_dict:
                out_fp.write(item + ':' + str(barcode_dict[item]) + '\t')
            out_fp.write('\n')
            barcode_dict = dict() 

        if barcode in barcode_dict: 
            barcode_dict[barcode] += 1
        else:
            barcode_dict[barcode] = 1
        prev_bin_id = bin_id
    
    in_fp.close()
    out_fp.close()
    # finished writing out_file1


    # start reading out_file1
    myprint('generating bcd12 file')
    in_file2 = out_file
    in_fp2 = open(in_file2, 'r')
    out_file2 = bcd12_file_list[tid] 
    out_fp2 = open(out_file2, 'w')


    # window1
    window_bcdset_dq1 = deque()
    rd_cnt_dq1 = deque()
    window1_all_bcd_set = set()
    for i in range(0, window_size):
        line = in_fp2.readline().strip().split(tab)
        bin_bcdset = set()
        rd_cnt = 0
        for j in range(1, len(line)):
            [barcode, count] = line[j].split(':')
            bin_bcdset.add(barcode)
            rd_cnt += int(count)

        rd_cnt_dq1.append(rd_cnt)
        window_bcdset_dq1.append(bin_bcdset)
        window1_all_bcd_set = window1_all_bcd_set.union(bin_bcdset)
    center1 = cal_weight_center(rd_cnt_dq1, 1)


        
    # window2
    window_bcdset_dq2 = deque()
    rd_cnt_dq2 = deque()
    window2_all_bcd_set = set()
    for i in range(0, window_size):
        line = in_fp2.readline().strip().split(tab)
        bin_bcdset = set()
        rd_cnt = 0
        for j in range(1, len(line)):
            [barcode, count] = line[j].split(':')
            bin_bcdset.add(barcode)
            rd_cnt += int(count)

        rd_cnt_dq2.append(rd_cnt)
        window_bcdset_dq2.append(bin_bcdset)
        window2_all_bcd_set = window2_all_bcd_set.union(bin_bcdset)
    center2 = cal_weight_center(rd_cnt_dq2, 2)

    window1_bc_cnt = len(window1_all_bcd_set)
    window2_bc_cnt = len(window2_all_bcd_set)
    overlapped_bc_cnt = len(window1_all_bcd_set.intersection(window2_all_bcd_set))
    middle_pos = (window_size-1) * bin_size

    out_fp2.write(str(middle_pos) + tab + str(window1_bc_cnt) + tab + str(window2_bc_cnt) + tab + str(overlapped_bc_cnt) + tab + '%.2f'%center1 + tab + '%.2f'%center2 + endl)

    line = in_fp2.readline().strip().split('\t')
    middle_pos += bin_size
    while line[0]:
        rd_cnt_dq1.popleft()
        rd_cnt_dq1.append(rd_cnt_dq2.popleft())

        window_bcdset_dq1.popleft()
        window_bcdset_dq1.append(window_bcdset_dq2.popleft())
        
        bin_bcdset = set()
        rd_cnt = 0
        for j in range (1, len(line)):
            [barcode, count] = line[j].split(':')
            bin_bcdset.add(barcode)
            rd_cnt += int(count)

        rd_cnt_dq2.append(rd_cnt)
        window_bcdset_dq2.append(bin_bcdset)

        center1 = cal_weight_center(rd_cnt_dq1, 1)
        center2 = cal_weight_center(rd_cnt_dq2, 2)
        
        window1_all_bcd_set = merge_window_all_bcdset(window_bcdset_dq1)
        window2_all_bcd_set = merge_window_all_bcdset(window_bcdset_dq2)

        window1_bc_cnt = len(window1_all_bcd_set)
        window2_bc_cnt = len(window2_all_bcd_set)
        overlapped_bc_cnt = len(window1_all_bcd_set.intersection(window2_all_bcd_set))

        out_fp2.write(str(middle_pos) + tab + str(window1_bc_cnt) + tab + str(window2_bc_cnt) + tab + str(overlapped_bc_cnt) + tab + '%.2f'%center1 + tab + '%.2f'%center2 + endl)

        line = in_fp2.readline().strip().split(tab)
        middle_pos += bin_size


    in_fp2.close()
    out_fp2.close()

    return
    

def cal_twin_empirical_dist(bcd12_file_list, out_file):

    twin_empirical_dist = list()
    out_fp = open(out_file, 'w')
    for bcd12_file in bcd12_file_list:
        bcd12_file = bcd12_file.strip()
        empirical_list1 =  get1twin_window_empirical_list(bcd12_file)
        twin_empirical_dist += empirical_list1

    for stat in twin_empirical_dist:
        out_fp.write(str(stat) + endl)

    out_fp.close()
        
    return twin_empirical_dist

def get1twin_window_empirical_list(bcd12_file):
    empirical_list1 = list()
    bcd12_fp = open(bcd12_file, 'r')
    while 1:
        line = bcd12_fp.readline()
        if not line:
            break
        line = line.strip().split(tab)
        if len(line) < 4:
            print "ERROR!" , line, bcd12_file
        m1 = float(line[1])
        m2 = float(line[2])
        m0 = float(line[3])
        if m1 == 0 or m2 == 0 or m0 == 0:
            continue
        stat = (m1 * m2) / m0

        empirical_list1.append(stat) 

    bcd12_fp.close()
    return empirical_list1
        

if __name__ == '__main__':
    main()
