#!/usr/bin/env python

import os
import sys

try:
    from scripts.my_utils import *
except ImportError:
    from my_utils import *

class Fragment:
    def __init__(self, attr_list):

        self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads, self.hp0, self.hp1, self.hp2, self.map_pos, self.map_qual = attr_list[0:12]
        self.n_left_weird_reads, self.n_right_weird_reads = attr_list[12:14]
        self.left_weird_reads_output, self.right_weird_reads_output, self.other_weird_reads_output = attr_list[14:17]

        self.tid       = int(self.tid)  
        self.start     = int(self.start)
        self.end       = int(self.end)
        self.length    = int(self.length )
        self.frag_id   = int(self.frag_id)
        self.num_reads = int(self.num_reads)
        self.hp0       = int(self.hp0)
        self.hp1       = int(self.hp1)
        self.hp2       = int(self.hp2)

        self.n_left_weird_reads  = int(self.n_left_weird_reads)
        self.n_right_weird_reads = int(self.n_right_weird_reads)
        self.key       = self.tid * FIX_LENGTH + self.start

    def output(self):
        outstring = '%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t' % (self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads, self.hp0, self.hp1, self.hp2, self.map_pos, self.map_qual)
        outstring += '%d\t%d\t%s\t%s\t%s' % (self.n_left_weird_reads, self.n_right_weird_reads, self.left_weird_reads_output, self.right_weird_reads_output, self.other_weird_reads_output)

        return outstring

    def key_start(self):
        return self.tid * FIX_LENGTH + self.start

    def key_end(self):
        return self.tid * FIX_LENGTH + self.end

    def gap_distances(self):
        map_pos = self.map_pos.strip(';').split(';')
        gap_distance_list = list()
        start_list = list()
        end_list = list()
        for pos in map_pos:
            start, end = pos.split(',')
            start_list.append(int(start))
            end_list.append(int(end))
        for i in range(1, len(end_list)):
            gap_distance_list.append(end_list[i] - start_list[i-1])

        return gap_distance_list

class CoreFragment:
    def __init__(self, attr_list):
        self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads = attr_list[0:7]
        self.tid       = int(self.tid)  
        self.start     = int(self.start)
        self.end       = int(self.end)
        self.length    = int(self.length )
        self.frag_id   = int(self.frag_id)
        self.num_reads = int(self.num_reads)

    def output(self):
        return '%d\t%d\t%d\t%d\t%s\t%d\t%d' % (self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads)

    def key_start(self):
        return self.tid * FIX_LENGTH + self.start

    def key_end(self):
        return self.tid * FIX_LENGTH + self.end

class ReadInfo:
    def __init__(self, attr_list):
        self.read_id, self.tid, self.start, self.end = attr_list[0:4]
        self.mapq, self.hptype, self.flag, self.insert_size = attr_list[4:8]
        self.mate_tid, self.mate_pos = attr_list[8:10]

        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)
        self.mapq = int(self.mapq)
        self.hptype = int(self.hptype)
        self.flag = int(self.flag)
        self.insert_size = int(self.insert_size)
        self.mate_tid = int(self.mate_tid)
        self.mate_pos = int(self.mate_pos)

    def output (self):
        outstring  = '[%s]%d]%d]%d]' % (self.read_id, self.tid, self.start, self.end)
        outstring += '%d]%d]%d]%d]' % (self.mapq, self.hptype, self.flag, self.insert_size)
        outstring += '%d]%d]' % (self.mate_tid, self.mate_pos)
        return outstring

def get_weird_reads_info_from_fragments(frm):

    left_weird_reads_info_list = list()
    right_weird_reads_info_list = list()
    other_weird_reads_output = list()

    if frm.left_weird_reads_output != '.':
        read_list = frm.left_weird_reads_output.lstrip('[').split('[')
        for read in read_list:
            read = read.rstrip(']').split(']')
            readinfo = ReadInfo(read)
            left_weird_reads_info_list.append(readinfo)

    if frm.right_weird_reads_output != '.':
        read_list = frm.right_weird_reads_output.lstrip('[').split('[')
        for read in read_list:
            read = read.rstrip(']').split(']')
            readinfo = ReadInfo(read)
            right_weird_reads_info_list.append(readinfo)

    if frm.other_weird_reads_output != '.':
        read_list = frm.other_weird_reads_output.lstrip('[').split('[')
        for read in read_list:
            read = read.rstrip(']').split(']')
            readinfo = ReadInfo(read)
            other_weird_reads_output.append(readinfo)

    return left_weird_reads_info_list, right_weird_reads_info_list, other_weird_reads_output


def convert_reads_info_list_to_string(reads_info_list):
    if len(reads_info_list) == 0:  
        output = '.' 
    else:
        output = ''
    for readinfo in reads_info_list:
        output += readinfo.output()
    
    return output 

def exist_read_pair_support(frm1, frm2, endtype1, endtype2): 

    if endtype1 == 'R_end' and frm1.n_right_weird_reads == 0: return False
    if endtype1 == 'L_end' and frm1.n_left_weird_reads == 0: return False
    if endtype2 == 'R_end' and frm2.n_right_weird_reads == 0: return False
    if endtype2 == 'L_end' and frm2.n_left_weird_reads == 0: return False 

    frm1_left_weird_reads_info_list, frm1_right_weird_reads_info_list, frm1_other_weird_reads_info_list = get_weird_reads_info_from_fragments(frm1)
    frm2_left_weird_reads_info_list, frm2_right_weird_reads_info_list, frm2_other_weird_reads_info_list = get_weird_reads_info_from_fragments(frm2)

    if endtype1 == 'R_end' and endtype2 == 'R_end':
        for read_info1 in frm1_right_weird_reads_info_list:
            for read_info2 in frm2_right_weird_reads_info_list:
                if read_info1.read_id == read_info2.read_id and read_info1.flag & 0x10 == 0 and read_info2.flag & 0x10 == 0:  return True


    if endtype1 == 'L_end' and endtype2 == 'L_end':
        for read_info1 in frm1_left_weird_reads_info_list:
            for read_info2 in frm2_left_weird_reads_info_list:
                if read_info1.read_id == read_info2.read_id and read_info1.flag & 0x10 and read_info2.flag & 0x10:  return True

    if endtype1 == 'L_end' and endtype2 == 'R_end':
        for read_info1 in frm1_left_weird_reads_info_list:
            for read_info2 in frm2_right_weird_reads_info_list:
                if read_info1.read_id == read_info2.read_id and read_info1.flag & 0x10 and read_info2.flag & 0x10 == 0: return True

    if endtype1 == 'R_end' and endtype2 == 'L_end':
        for read_info1 in frm1_right_weird_reads_info_list:
            for read_info2 in frm2_left_weird_reads_info_list:
                if read_info1.read_id == read_info2.read_id and read_info1.flag & 0x10 == 0 and read_info2.flag & 0x10: return True

    return False

def read_bcd22_file(bcd22_file, min_frm_length = 0):
    frm_list = list()
    bcd22_fp = open(bcd22_file, 'r')
    while 1:
        line = bcd22_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        frm = Fragment(line)
        if frm.length < min_frm_length: continue
        frm_list.append(frm)

    bcd22_fp.close()
    return frm_list

def extract_frm_from_bcd22_file(bcd22_file, frm_id_set):

    frm_list = list()
    bcd22_fp = open(bcd22_file, 'r')
    while 1:
        line = bcd22_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        frm = Fragment(line)
        if frm.frag_id in frm_id_set: frm_list.append(frm)

    bcd22_fp.close()

    return frm_list

def read_bcd22_file_core(bcd22_file, min_frm_length = 0):
    frm_list = list()
    bcd22_fp = open(bcd22_file, 'r')
    while 1:
        line = bcd22_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line = line.strip().split(tab)
        frm = CoreFragment(line)
        if frm.length < min_frm_length: continue
        frm_list.append(frm)

    bcd22_fp.close()
    return frm_list


