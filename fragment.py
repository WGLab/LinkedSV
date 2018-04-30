#!/usr/bin/env python

import os
import sys
from my_utils import *

class Fragment:
    def __init__(self, attr_list):
        self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads, self.hp0, self.hp1, self.hp2, self.map_pos = attr_list[0:11]
        self.tid       = int(self.tid)  
        self.start     = int(self.start)
        self.end       = int(self.end)
        self.length    = int(self.length )
        self.frag_id   = int(self.frag_id)
        self.num_reads = int(self.num_reads)
        self.hp0       = int(self.hp0)
        self.hp1       = int(self.hp1)
        self.hp2       = int(self.hp2)
        self.key       = self.tid * FIX_LENGTH + self.start

    def output(self):
        return '%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s' % (self.tid, self.start, self.end, self.length, self.bcd, self.frag_id, self.num_reads, self.hp0, self.hp1, self.hp2, self.map_pos)

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


if __name__ == '__main__':
    main()
