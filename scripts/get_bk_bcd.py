#!/usr/bin/env python

from my_utils import *

arg = sys.argv[1:]
arg.reverse()

usage = 'python ' + __file__ + ' ' + '<bcd11.list> <bk.bed> <out.peak> <win_size> <bin_size>'
argc  = 5

class bk_candidate:
    def __init__(self, tid, start, end, info):
        self.tid = int(tid)
        self.start = int(start)
        self.end = int(end)
        self.info = info
        self.left_bcd_set = set()
        self.right_bcd_set = set()

def main():
    if len(arg) < argc:
        print usage
        sys.exit()
    
    bcd11_list_file = arg.pop()
    bk_bed_file     = arg.pop()
    out_peak_file   = arg.pop()
    win_size        = int(arg.pop())
    bin_size        = int(arg.pop())

    bcd11_list_fp = open(bcd11_list_file, 'r') 
    lines = list(bcd11_list_fp)
    bcd11_list = list()
    for line in lines:
        bcd11_list.append(line.strip())
        
    get_bk_bcd(bcd11_list, bk_bed_file, out_peak_file, win_size, bin_size) 

    return

def get_bk_candidates(bcd11_list, bk_bed_file):

    bk_bed_fp = open(bk_bed_file, 'r')
    lines = list(bk_bed_fp)
    bk_bed_fp.close()
    bk_cand_db = dict()

    for tid in range(0, len(bcd11_list)):
        bk_cand_db[tid] = list()

    for line in lines:
        line = line.strip().split(tab)
        tid = int(line[0])
        bk_cand = bk_candidate(line[0], line[1], line[2], tab.join(line[3:]))
        bk_cand_db[tid].append(bk_cand)

    return bk_cand_db

def get_bk_bcd(bcd11_list, bk_bed_file, out_peak_file, win_size, bin_size):

    out_peak_fp = open(out_peak_file, 'w')
    
    myprint ('reading bk file: ' + bk_bed_file)
    bk_cand_db = get_bk_candidates(bcd11_list, bk_bed_file)
        
    myprint ('extracting barcodes for breakpoint candidates...')
    for tid in range(0, len(bcd11_list)):
        myprint ('processing ' + str(tid) + '-th chr')
        bcd11_file = bcd11_list[tid]
        bcd11_lines = list()
        bcd11_pos_list = list()
        bcd11_fp = open(bcd11_file, 'r')
        while 1:
            line = bcd11_fp.readline()
            if not line:
                break
            bcd11_lines.append(line)
            line1 = line.strip().split(tab)
            pos = int(line1[0]) * bin_size
            bcd11_pos_list.append(pos)
        bcd11_fp.close()


        for i in range(0, len(bk_cand_db[tid])):
            bk_cand = bk_cand_db[tid][i]
            left_end_pos = bk_cand.end
            left_end_index = bisect.bisect_left(bcd11_pos_list, left_end_pos)
            left_start_index = left_end_index - win_size

            right_start_index = left_end_index
            right_end_index = right_start_index + win_size

            left_bcd_set = set()
            right_bcd_set = set()

            for line in bcd11_lines[left_start_index:left_end_index]:
                line = line.strip().split(tab)
                for item in line[1:]:
                    left_bcd_set.add(item.split(':')[0])

            for line in bcd11_lines[right_start_index:right_end_index]:
                line = line.strip().split(tab)
                for item in line[1:]:
                    right_bcd_set.add(item.split(':')[0])

            bk_cand_db[tid][i].left_bcd_set = left_bcd_set        
            bk_cand_db[tid][i].right_bcd_set = right_bcd_set        
            intersect_bcd_set = bk_cand.left_bcd_set.intersection(bk_cand.right_bcd_set) 

            info2 = 'left_bcd_cnt=' + str(len(bk_cand.left_bcd_set)) + ';right_bcd_cnt=' + str(len(bk_cand.right_bcd_set)) + ';overlapped_bcd_cnt=' +  str(len(intersect_bcd_set))
            out_peak_fp.write(str(bk_cand.tid) + tab + str(bk_cand.start) + tab + str(bk_cand.end) + tab + bk_cand.info + tab + info2 + tab )

            for bcd in bk_cand.left_bcd_set:
                out_peak_fp.write(bcd + ',')
            out_peak_fp.write(tab)

            for bcd in bk_cand.right_bcd_set:
                out_peak_fp.write(bcd + ',')
            out_peak_fp.write(endl)
           
    return     


if __name__ == '__main__':
    main()
