#!/usr/bin/env python

import numpy as np
from my_utils import *
from fragment import *
from cluster_weird_reads import *
import global_distribution
import math

tab = '\t'
endl = '\n'

class Bcd21:

    def __init__(self, attr_list):
        #barcode    tid start_post  end_pos map_qual    hap_type    ReadID  flag    n_left_clip n_right_clip    insert_size mate_tid    mate_pos
        self.tid, self.start, self.end, self.mapq, self.bcd, self.hptype, self.read_id, self.flag, self.n_left_clip, self.n_right_clip, self.insert_size, self.mate_tid, self.mate_pos = attr_list[0:13]
        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)
        self.mapq = int(self.mapq)
        self.hptype = int(self.hptype)
        self.flag = int(self.flag)
        self.n_left_clip = int(self.n_left_clip)
        self.n_right_clip = int(self.n_right_clip)
        self.insert_size = int(self.insert_size)
        self.mate_tid = int(self.mate_tid)
        self.mate_pos = int(self.mate_pos)

    def key_start(self):
        return self.tid * FIX_LENGTH + self.start

    def key_end(self):
        return self.tid * FIX_LENGTH + self.end

    def output_read_info (self):
        output = '[%s]%d]%d]%d]' % (self.read_id, self.tid, self.start, self.end) 
        output += '%d]%d]%d]%d]' % (self.mapq, self.hptype, self.flag, self.insert_size)
        output += '%d]%d]' % (self.mate_tid, self.mate_pos)
        return output 

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()
    cluster_reads(args, dbo_args, endpoint_args)

    return


def cluster_reads(args, dbo_args, endpoint_args):

    weired_readname_dict = dict() 

    ## first round ##
    length_cut = 200 * 1000 # first round length cut is 200k
    if args.run_from_begining == False and check_file_exists(endpoint_args.tmpbcd22_file):
        myprint('tmpbcd22_file existed, skipped first round clustering')
    else:
        myprint('first round clustering reads, length cut is %d, output file is: %s' % (length_cut, endpoint_args.tmpbcd22_file))
        bcd21_to_bcd22_file (endpoint_args.bcd21_file, endpoint_args.tmpbcd22_file, length_cut, weired_readname_dict)

    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args, endpoint_args.tmpbcd22_file, is_fast_mode = True)

    # get weird reads #

    myprint('getting weird read names')

    tid2chrname_list, chrname2tid_dict = get_chrnames(args.faidx_file)

    weired_readname_dict = get_weired_readname_dict(chrname2tid_dict, args.weired_reads_file)

    myprint('finished getting weired read names')

    ## second round ## 

    if args.run_from_begining == False and check_file_exists(endpoint_args.bcd22_file):
        myprint('bcd22_file existed, skipped second round clustering')
    else:
        myprint('second round clustering reads, length cut is %d, output file is: %s' % (args.gap_distance_cutoff, endpoint_args.bcd22_file))
        bcd21_to_bcd22_file (endpoint_args.bcd21_file, endpoint_args.bcd22_file, args.gap_distance_cutoff, weired_readname_dict)

    del weired_readname_dict

    gc.collect()

    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args, endpoint_args.bcd22_file, is_fast_mode = False)

    return


def bcd21_to_bcd22_file(bcd21_file, bcd22_file, length_cut, weired_readname_dict, is_fast_mode = False):

    bcd21_fp = open(bcd21_file, 'r')
    bcd22_fp = open(bcd22_file, 'w')

    bcd22_fp.write('#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos\tnum_left_weird_reads\tnum_right_weird_reads\tleft_weird_reads_info\tright_weird_reads_info\tother_weird_reads_info; gap_distance_cut_off=%d\n' % (length_cut))

    fragment_bcd21_list = list()

    frag_id = 0
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        new_bcd21_term = Bcd21(line)
        if len(fragment_bcd21_list) == 0:
            fragment_bcd21_list.append(new_bcd21_term)
            continue
        if new_bcd21_term.bcd == fragment_bcd21_list[-1].bcd and new_bcd21_term.key_start() - fragment_bcd21_list[-1].key_end() < length_cut:
            fragment_bcd21_list.append(new_bcd21_term)
        else:
            bcd22 = convert_bcd21list_to_bcd22 (fragment_bcd21_list, frag_id, is_fast_mode, weired_readname_dict)
            bcd22_fp.write(bcd22.output() + endl)
            frag_id += 1
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)
            if frag_id % 1000000 == 0: myprint('grouped %d fragments' % frag_id)

    if len(fragment_bcd21_list) > 0:
        bcd22 = convert_bcd21list_to_bcd22 (fragment_bcd21_list, frag_id, is_fast_mode, weired_readname_dict)
        bcd22_fp.write(bcd22.output() + endl)

    bcd21_fp.close()
    bcd22_fp.close()

    return


def convert_bcd21list_to_bcd22(fragment_bcd21_list, frag_id, is_fast_mode, weired_readname_dict):
    
    frm_tid   = fragment_bcd21_list[0].tid 
    frm_start = fragment_bcd21_list[0].start
    frm_end   = fragment_bcd21_list[-1].end
    frm_length = frm_end - frm_start
    frm_bcd = fragment_bcd21_list[0].bcd
    num_reads = len(fragment_bcd21_list)
    hptype = [0] * 3
    map_pos = ''

    for bcd21 in fragment_bcd21_list:
        hptype[bcd21.hptype] += 1
        map_pos += '%d,%d;' % (bcd21.start, bcd21.end)
    map_pos = map_pos.strip(';') 

    if is_fast_mode:
        attr_list = [frm_tid, frm_start, frm_end, frm_length, frm_bcd, frag_id, num_reads, hptype[0], hptype[1], hptype[2], map_pos]  
        attr_list += [0, 0, '.', '.', '.']
        return Fragment(attr_list)
        
    n_left_weird_reads = 0
    n_right_weird_reads = 0
    weird_reads_list = list()

    for bcd21 in fragment_bcd21_list:
        if bcd21.read_id in weired_readname_dict:
            weird_reads_list.append(bcd21)

    left_side_weird_reads_list = list()
    right_side_weird_reads_list = list()
    other_weird_reads_list = list()

    for bcd21 in weird_reads_list: 
        if bcd21.end - frm_start < 400 and bcd21.flag & 0x10: 
            n_left_weird_reads +=1
            left_side_weird_reads_list.append(bcd21)
        elif frm_end - bcd21.start < 400 and (not bcd21.flag & 0x10): 
            n_right_weird_reads +=1
            right_side_weird_reads_list.append(bcd21)
        else:
            other_weird_reads_list.append(bcd21)

    left_weird_reads_output = ''
    right_weird_reads_output = ''
    other_weird_reads_output = ''

    for bcd21 in left_side_weird_reads_list:
        left_weird_reads_output  += bcd21.output_read_info()

    for bcd21 in right_side_weird_reads_list:
        right_weird_reads_output += bcd21.output_read_info()

    for bcd21 in other_weird_reads_list:
        other_weird_reads_output += bcd21.output_read_info()

    if left_weird_reads_output == '': left_weird_reads_output = '.'
    if right_weird_reads_output == '': right_weird_reads_output = '.'
    if other_weird_reads_output == '': other_weird_reads_output = '.'

    attr_list = [frm_tid, frm_start, frm_end, frm_length, frm_bcd, frag_id, num_reads, hptype[0], hptype[1], hptype[2], map_pos ]  
    attr_list += [n_left_weird_reads, n_right_weird_reads, left_weird_reads_output, right_weird_reads_output, other_weird_reads_output]

    return Fragment(attr_list)

if __name__ == '__main__':
    main()
