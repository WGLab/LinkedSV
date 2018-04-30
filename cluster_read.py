#!/usr/bin/env python

import numpy as np
from my_utils import *
from fragment import *
import global_distribution
import math

tab = '\t'
endl = '\n'
def main():

    args, dbo_args, endpoint_args = parse_user_arguments()
    cluster_reads(args, dbo_args, endpoint_args)

    return

class Bcd21:
    def __init__(self, attr_list):
        self.tid, self.start, self.end, self.mapq, self.bcd, self.hptype, self.prev_read_pos, self.next_read_pos, self.prev_distance, self.next_distance = attr_list[0:10]
        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)
        self.mapq = int(self.mapq)
        self.hptype = int(self.hptype)
        self.prev_distance = int(self.prev_distance)
        self.next_distance = int(self.next_distance) 

def bcd21_to_bcd22_file(bcd21_file, bcd22_file, length_cut):

    bcd21_fp = open(bcd21_file, 'r')
    bcd22_fp = open(bcd22_file, 'w')
    bcd22_fp.write('#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos, gap_distance_cut_off=%d\n' % (length_cut))
    fragment_bcd21_list = list()
    line = bcd21_fp.readline().strip().split(tab)
    curr_bcd21_term = Bcd21(line)
    fragment_bcd21_list.append(curr_bcd21_term) 
    frag_id = 0
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        new_bcd21_term = Bcd21(line)
        if new_bcd21_term.bcd == curr_bcd21_term.bcd and new_bcd21_term.prev_distance < length_cut and new_bcd21_term.prev_read_pos != '0:-1':
            fragment_bcd21_list.append(new_bcd21_term)
        else:
            bcd22 = convert_bcd21list_to_bcd22 (fragment_bcd21_list, frag_id)
            frag_id += 1
            bcd22_fp.write(bcd22.output() + endl)
            if frag_id % 1000000 == 0: myprint('grouped %d fragments' % frag_id)
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)
            curr_bcd21_term = new_bcd21_term

    if len(fragment_bcd21_list) > 0:
        bcd22 = convert_bcd21list_to_bcd22 (fragment_bcd21_list, frag_id)
        bcd22_fp.write(bcd22.output() + endl)

    bcd21_fp.close()
    bcd22_fp.close()

    return

def cluster_reads(args, dbo_args, endpoint_args):

    ## first round ##
    length_cut = 100 * 1000 # first round length cut is 100k
    myprint('first round clustering reads, length cut is %d, output file is: %s' % (length_cut, endpoint_args.tmpbcd22_file))
    bcd21_to_bcd22_file(endpoint_args.bcd21_file, endpoint_args.tmpbcd22_file, length_cut)
    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args)

    myprint('second round clustering reads, length cut is %d, output file is: %s' % (args.gap_distance_cutoff, endpoint_args.bcd22_file))
    bcd21_to_bcd22_file(endpoint_args.bcd21_file, endpoint_args.bcd22_file, args.gap_distance_cutoff)
    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args)
    
    return

def convert_bcd21list_to_bcd22(fragment_bcd21_list, frag_id):
    
    frm_tid = fragment_bcd21_list[0].tid 
    frm_start = fragment_bcd21_list[0].start
    frm_end = fragment_bcd21_list[-1].end
    frm_length = frm_end - frm_start
    frm_bcd = fragment_bcd21_list[0].bcd
    num_reads = len(fragment_bcd21_list)
    hptype = [0] * 3
    map_pos = ''
    for bcd21 in fragment_bcd21_list:
        hptype[bcd21.hptype] += 1
        map_pos += '%d,%d;' % (bcd21.start, bcd21.end)
    
    map_pos.strip(';') 

    attr_list = [frm_tid, frm_start, frm_end, frm_length, frm_bcd, frag_id, num_reads, hptype[0], hptype[1], hptype[2], map_pos]  

    return Fragment(attr_list)
    
if __name__ == '__main__':
    main()
