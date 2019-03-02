#!/usr/bin/env python

import os
import sys
import bisect
from my_utils import *
from fragment import *

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]
arg.reverse()

usage = 'python ' + __file__ + ' ' + '<in.bcd22> <in.bedpe> <faidx_file> <output_file>'
argc  =  4

MAX_FRM_LENGTH = int(2e5)
MIN_FRM_LENGTH = int(2000)

'''
class Fragment:
    def __init__(self, bcd22_line):

        line = bcd22_line.strip().split(tab)
        self.tid       = int(line[0])  # chr number , start from 0
        self.start     = int(line[1])
        self.end       = int(line[2])
        self.length    = int(line[3])
        self.bcd       = line[4]
        self.frag_id   = int(line[5])
        self.num_reads = int(line[6])
        self.map_pos   = line[7]
        self.key       = self.tid * FIX_LENGTH + self.start

    def output(self, chrname_list):
        chrm = chrname_list[self.tid]
        return '%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s\n' % (chrm, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads, self.map_pos)
        
'''

class FragmentwithKey(Fragment):
    def __init__(self, attr_list):
        Fragment.__init__(self, attr_list)
        self.key = self.tid * FIX_LENGTH + self.start

class Bedpe:
    def __init__(self, bed_line):

        line = bed_line.strip().split(tab)
        self.chrm1  = line[0]
        self.start1 = int(line[1])
        self.end1   = int(line[2]) 
        self.chrm2  = line[3]
        self.start2 = int(line[4])
        self.end2   = int(line[5]) 
        if self.end1 - self.start1 < 20000:
            mean1 = (self.start1 + self.end1) / 2
            self.start1 = mean1 - 10000
            self.end1 = mean1 + 10000
        if self.end2 - self.start2 < 20000:
            mean2 = (self.start2 + self.end2) / 2
            self.start2 = mean2 - 10000
            self.end2 = mean2 + 10000


def main():
    
    if len(arg) < argc:
        print usage
        sys.exit()

    in_bcd22_file = os.path.abspath(arg.pop())
    in_bedpe_file = os.path.abspath(arg.pop())
    faidx_file    = os.path.abspath(arg.pop())
    out_file      = os.path.abspath(arg.pop())

    chrname_list, chrname_dict = get_chrnames(faidx_file)
    frm_db, frag_key_db = get_frm_db_from_bcd22_file(in_bcd22_file) 
    bedpe_list = read_bedpe_file(in_bedpe_file)


    out_fp = open(out_file, 'w')
    myprint ('start searching shared barcodes')
    i = 0 
    all_output1 = ''
    all_output2 = ''
    for bedpe in bedpe_list:
        output1, output2 = search_for_shared_barcodes(frm_db, frag_key_db, bedpe, chrname_dict, chrname_list)
        if len(output1) > 1:
            all_output1 += output1
            all_output2 += output2
        i += 1
        if i % 1e3 == 0: myprint('searched %d bed pairs' % i)
    
    out_fp.write(all_output1)
    out_fp.write('\n\n\n\n\n')
    out_fp.write(all_output2)
    out_fp.close()
    myprint ('finished searching shared barcodes')
    return

def mim_frm_gap(support_frm_list):
    if len(support_frm_list) < 2:
        return 0, 0

    start_gap_list = list()
    start_sorted_frm_list = sorted(support_frm_list, key = lambda frm: frm.start)
    for i in range(1, len(start_sorted_frm_list)): 
        gap = start_sorted_frm_list[i].start -start_sorted_frm_list[i-1].start
        start_gap_list.append(gap)


    end_gap_list = list()
    end_sorted_frm_list = sorted(support_frm_list, key = lambda frm: frm.end)
    for i in range(1, len(end_sorted_frm_list)): 
        gap = end_sorted_frm_list[i].end -end_sorted_frm_list[i-1].end
        end_gap_list.append(gap)

    min_start_gap = min(start_gap_list)
    min_end_gap = min(end_gap_list)
    return min_start_gap, min_end_gap, start_gap_list, end_gap_list

def search_for_shared_barcodes(frm_db, frag_key_db, bedpe, chrname_dict, chrname_list):

    frm_list1 = get_interval_frm_list(frm_db, frag_key_db, bedpe.chrm1, bedpe.start1, bedpe.end1, chrname_dict)
    frm_list2 = get_interval_frm_list(frm_db, frag_key_db, bedpe.chrm2, bedpe.start2, bedpe.end2, chrname_dict)
    
    frm_bcd_set1 = set()
    frm_bcd_set2 = set()
    for frm in frm_list1:
        frm_bcd_set1.add(frm.bcd)
    for frm in frm_list2:
        frm_bcd_set2.add(frm.bcd)

    shared_bcd_set = frm_bcd_set1.intersection(frm_bcd_set2)
    if len(shared_bcd_set) < 2:
        return '' , ''

    support_frm_list1 = list()
    support_frm_list2 = list()

    for frm in frm_list1:
        if frm.bcd in shared_bcd_set and frm.length:
            support_frm_list1.append(frm)
    for frm in frm_list2:
        if frm.bcd in shared_bcd_set and frm.length:
            support_frm_list2.append(frm)

    min_start_gap1, min_end_gap1, start_gap_list1, end_gap_list1 = mim_frm_gap(support_frm_list1)
    min_start_gap2, min_end_gap2, start_gap_list2, end_gap_list2 = mim_frm_gap(support_frm_list1)


    output2 = '\n\nsupport fragments for breakpoint %s:%d-%d;%s:%d-%d (bk1)\n\n' % (bedpe.chrm1, bedpe.start1, bedpe.end1, bedpe.chrm2, bedpe.start2, bedpe.end2)
    for frm in frm_list1:
        if frm.bcd in shared_bcd_set and frm.length:
            output2 += frm.output()

    output2 += '\n\nsupport fragments for breakpoint %s:%d-%d;%s:%d-%d (bk2)\n\n' % (bedpe.chrm1, bedpe.start1, bedpe.end1, bedpe.chrm2, bedpe.start2, bedpe.end2)

    for frm in frm_list2:
        if frm.bcd in shared_bcd_set and frm.length:
            output2 += frm.output()

    output2 += '\n\n' 

    output1 = '%s:%d-%d\t%s:%d-%d\t%d\tsupport_info\t' % (bedpe.chrm1, bedpe.start1, bedpe.end1, bedpe.chrm2, bedpe.start2, bedpe.end2, len(shared_bcd_set))

    total_num_reads1 = 0
    total_num_reads2 = 0
    average_frm_length1 = 0
    average_frm_length2 = 0
    for frm in support_frm_list1:
        total_num_reads1 += frm.num_reads
        average_frm_length1 += frm.length
    for frm in support_frm_list2:
        total_num_reads2 += frm.num_reads
        average_frm_length2 += frm.length

    average_frm_length1 = float(average_frm_length1)/len(support_frm_list1)
    average_frm_length2 = float(average_frm_length2)/len(support_frm_list2)

    output1 += '%d\t%d\t%.f\t%.f\t%d\t%d\t%d\t%d\t' % (total_num_reads1, total_num_reads2, average_frm_length1, average_frm_length2, min_start_gap1, min_end_gap1, min_start_gap2, min_end_gap2)

    for frm in support_frm_list1:
        output1 += str(frm.num_reads) + ','
    output1 += tab
    for frm in support_frm_list2:
        output1 += str(frm.num_reads) + ','
    output1 += tab

    for frm in support_frm_list1:
        output1 += str(frm.length) + ','
    output1 += tab

    for frm in support_frm_list2:
        output1 += str(frm.length) + ','

    
    output1 += endl
    return output1, output2

def get_interval_frm_list(frm_db, frag_key_db, chrm, start, end, chrname_dict):

    frm_list = list()
    tid = chrname_dict[chrm] 
    search_start = start - MAX_FRM_LENGTH
    if search_start < 0:
        search_start = 0

    start_key = tid * FIX_LENGTH + search_start 
    end_key = tid * FIX_LENGTH + end

    start_index = bisect.bisect_left(frag_key_db[tid],start_key)
    end_index = bisect.bisect_right(frag_key_db[tid], end_key)
    for i in range(start_index, end_index):
        frm = frm_db[tid][i]
        if frm.end > start and frm.start < end and frm.length > MIN_FRM_LENGTH:
            frm_list.append(frm) 
    
    return frm_list

def read_bedpe_file(bedpe_file):
    bedpe_list = list()
    bedpe_fp = open(bedpe_file, 'r')
    while 1:
        line = bedpe_fp.readline()
        if not line:
            break
        bedpe = Bedpe(line)
        bedpe_list.append(bedpe)

    bedpe_fp.close()
    return bedpe_list

def get_frm_db_from_bcd22_file(in_bcd22_file):

    myprint ('reading bcd22 file')
    frag_db = dict()
    bcd22_fp = open(in_bcd22_file, 'r')
    line = bcd22_fp.readline() # skip the header line
    while 1:
        line = bcd22_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        frm = FragmentwithKey(line) 
        if frm.tid not in frag_db:
            frag_db[frm.tid] = list()

        frag_db[frm.tid].append(frm)
    
    bcd22_fp.close() 
    frag_key_db = dict()
    
    myprint('sorting frm list...')
    for tid in frag_db:
        frag_key_db[tid] = list()
        frag_db[tid].sort(key=lambda x: x.key)
        for frm in frag_db[tid]: 
            frag_key_db[tid].append(frm.key)

    myprint('finished sorting frm list...')
    return frag_db, frag_key_db
            
if __name__ == '__main__':
    main()
