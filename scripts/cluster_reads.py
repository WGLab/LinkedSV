#!/usr/bin/env python

import numpy as np
import math
import gzip

try:
    from scripts.my_utils import *
except ImportError:
    from my_utils import *

try:
    from scripts.fragment import *
except ImportError:
    from fragment import *
try:
    from scripts.cluster_weird_reads import *
except ImportError:
    from cluster_weird_reads import *
try:
    import scripts.global_distribution
except ImportError:
    import global_distribution


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


    ## first round ##
    length_cut = 50 * 1000 # first round length cut is 50k
    if args.run_from_begining == False and check_file_exists(endpoint_args.tmpbcd22_file):
        myprint('tmpbcd22_file existed, skipped first round clustering')
    else:
        myprint('first round clustering reads, length cut is %d, output file is: %s' % (length_cut, endpoint_args.tmpbcd22_file))
        bcd21_to_bcd22_file (args, endpoint_args.bcd21_file, endpoint_args.tmpbcd22_file, length_cut, set(), is_fast_mode = True, with_weird_reads = False)

    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args, endpoint_args.tmpbcd22_file, is_fast_mode = True)

    if args.run_from_begining == False and check_file_exists(endpoint_args.bcd22_file):
        myprint ('bcd22_file existed, skipped second round clustering')
        global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args, endpoint_args.bcd22_file, is_fast_mode = False)

        if args.rm_temp_files:
            if check_file_exists(endpoint_args.tmpbcd22_file): os.remove(endpoint_args.tmpbcd22_file)

        return

    ## split bcd21 file ##
    
    bcd21_split_file = endpoint_args.bcd21_file + '.split'
    if args.run_from_begining == False and check_file_exists(args.num_split_bcd21_file) and check_file_exists(bcd21_split_file):
        myprint('bcd21_split_file existed, skipped splitting')
    else:
        myprint('splitting bcd21 file')
        cmd = 'python ' + args.split_weird_reads_program + ' ' + endpoint_args.bcd21_file + ' ' + args.weird_reads_file + ' ' + args.num_split_bcd21_file
        myprint ('running command: %s' % cmd)
        os.system(cmd)

    num_split_bcd21_fp = open(args.num_split_bcd21_file, 'r')
    n_split_file = int(num_split_bcd21_fp.readline().strip())
    num_split_bcd21_fp.close()

    args.temp_file_list.append(args.num_split_bcd21_file)

    if n_split_file <= 0: 
        myprint('ERROR!')
        sys.exit()

    split_barcode_set = set()
    bcd21_split_fp = open(bcd21_split_file, 'r')
    lines = list(bcd21_split_fp)
    for line in lines:
        bcd = line.strip().split(tab)[0]
        split_barcode_set.add(bcd)
    bcd21_split_fp.close()

    args.temp_file_list.append(bcd21_split_file)

    ## second round ## 
    if args.run_from_begining == False and check_file_exists(endpoint_args.bcd22_file):
        myprint('bcd22_file existed, skipped second round clustering')
    else:
        myprint('second round clustering reads, length cut is %d, output file is: %s' % (args.gap_distance_cutoff, endpoint_args.bcd22_file))

        bcd21_to_bcd22_file (args, endpoint_args.bcd21_file, endpoint_args.bcd22_file, args.gap_distance_cutoff, split_barcode_set, is_fast_mode = False, with_weird_reads=True)

    global_distribution.estimate_global_distribution (args, dbo_args, endpoint_args, endpoint_args.bcd22_file, is_fast_mode = False)

    if args.rm_temp_files and check_file_exists(endpoint_args.bcd22_file):
        if check_file_exists(endpoint_args.tmpbcd22_file): os.remove(endpoint_args.tmpbcd22_file)

        for temp_file in args.temp_file_list:
            if check_file_exists(temp_file): os.remove(temp_file) 
        args.temp_file_list = list()

    return


def bcd21_to_bcd22_file(args, bcd21_file, bcd22_file, length_cut, split_barcode_set, is_fast_mode = True, with_weird_reads = False):

    if args.is_wgs:
        min_num_good_reads = 6 
    else:
        min_num_good_reads = 3 

    if args.user_defined_min_reads_in_fragment > 0:
        min_num_good_reads = args.user_defined_min_reads_in_fragment

    weird_readname_dict = dict()

    if with_weird_reads == True:
        tid2chrname_list, chrname2tid_dict = get_chrnames(args.faidx_file)
        fid = 1
        weird_reads_file = args.weird_reads_file + '.split%d' % (fid) 
        args.temp_file_list.append(weird_reads_file)
        myprint('getting weird read names from file: %s' % weird_reads_file)
        weird_readname_dict = get_weird_readname_dict (chrname2tid_dict, weird_reads_file)
        myprint('finished getting weird read names')

    bcd21_fp = my_utils.gzopen(bcd21_file, 'r')

    bcd22_fp = open(bcd22_file, 'w')
    bcd22_header = '#tid\tfrag_start\tfrag_end\tfrag_length\tfrag_barcode\tfrag_ID\tnum_reads\thptype0\thptype1\thptype2\tmap_pos\tmap_qual\tnum_left_weird_reads\tnum_right_weird_reads\tleft_weird_reads_info\tright_weird_reads_info\tother_weird_reads_info; gap_distance_cut_off=%d\n' % (length_cut)
    bcd22_fp.write(bcd22_header)

    fragment_bcd21_list = list()

    frag_id = 0
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        new_bcd21_term = Bcd21(line)
        bcd = new_bcd21_term.bcd
        read_id = new_bcd21_term.read_id

        if with_weird_reads == True and bcd in split_barcode_set:

            bcd22 = convert_bcd21list_to_bcd22(args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
            if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)
            frag_id += 1
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)

            fid += 1
            split_barcode_set.remove(bcd)
            weird_reads_file = args.weird_reads_file + '.split%d' % (fid) 
            args.temp_file_list.append(weird_reads_file) 
            del weird_readname_dict
            weird_readname_dict = get_weird_readname_dict (chrname2tid_dict, weird_reads_file)
            continue

        if len(fragment_bcd21_list) == 0 or ( new_bcd21_term.bcd == fragment_bcd21_list[-1].bcd and new_bcd21_term.key_start() - fragment_bcd21_list[-1].key_end() < length_cut):
            fragment_bcd21_list.append(new_bcd21_term)
        else:
            bcd22 = convert_bcd21list_to_bcd22 (args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
            if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)
            frag_id += 1
            fragment_bcd21_list = list()
            fragment_bcd21_list.append(new_bcd21_term)

            if frag_id % 1000000 == 0: myprint('grouped %d fragments' % frag_id)

    if len(fragment_bcd21_list) > 0:
        bcd22 = convert_bcd21list_to_bcd22 (args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict)
        if n_good_reads(bcd22) >= min_num_good_reads: bcd22_fp.write(bcd22.output() + endl)

    bcd21_fp.close()
    bcd22_fp.close()

    if with_weird_reads == True: del weird_readname_dict

    gc.collect()

    return

def n_good_reads(frm):

    num_good_reads = 0
    for i in range(0, len(frm.map_qual)):
        if frm.map_qual[i] == '1': num_good_reads += 1
    return num_good_reads


def convert_bcd21list_to_bcd22(args, fragment_bcd21_list, frag_id, is_fast_mode, weird_readname_dict):
    
    frm_tid   = fragment_bcd21_list[0].tid 
    frm_start = fragment_bcd21_list[0].start
    frm_end   = fragment_bcd21_list[-1].end
    frm_length = frm_end - frm_start
    frm_bcd = fragment_bcd21_list[0].bcd
    num_reads = len(fragment_bcd21_list)
    hptype = [0] * 3
    map_pos = ''
    map_qual = ''

    bad_flag = 256 + 1024 + 2048
    for bcd21 in fragment_bcd21_list:
        hptype[bcd21.hptype] += 1
        map_pos += '%d,%d;' % (bcd21.start, bcd21.end)
        if bcd21.mapq < args.min_mapq or bcd21.flag & bad_flag: 
            q = '0'
        else:
            q = '1'
        map_qual += q

    map_pos = map_pos.strip(';') 

    if is_fast_mode:
        attr_list = [frm_tid, frm_start, frm_end, frm_length, frm_bcd, frag_id, num_reads, hptype[0], hptype[1], hptype[2], map_pos, map_qual]
        attr_list += [0, 0, '.', '.', '.']
        return Fragment(attr_list)
        
    n_left_weird_reads = 0
    n_right_weird_reads = 0
    weird_reads_list = list()

    for bcd21 in fragment_bcd21_list:
        if bcd21.read_id in weird_readname_dict:
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

    attr_list = [frm_tid, frm_start, frm_end, frm_length, frm_bcd, frag_id, num_reads, hptype[0], hptype[1], hptype[2], map_pos, map_qual]  
    attr_list += [n_left_weird_reads, n_right_weird_reads, left_weird_reads_output, right_weird_reads_output, other_weird_reads_output]

    frm = Fragment(attr_list)
    return frm

def calculate_num_reads_from_bcd21_file(in_bcd21_file, min_mapq):

    n_reads = 0
    in_bcd21_fp = gzopen(in_bcd21_file, 'r')
    while 1:
        line = in_bcd21_fp.readline()
        if not line: break
        line = line.strip().split(tab) 
        bcd21 = Bcd21(line)
        if bcd21.flag & 256 or bcd21.flag & 1024 or bcd21.flag & 2048: continue
        if bcd21.mapq < min_mapq: continue
        n_reads += 1 

    in_bcd21_fp.close()

    return n_reads


if __name__ == '__main__':
    main()
