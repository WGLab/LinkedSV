#!/usr/bin/env python

import os 
import sys
import bisect
import psutil
import gc
import math
import gzip
from datetime import datetime
from arguments import *

tab = '\t'
endl = '\n'
endll = '\n\n'

ChrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'] 
ChrSet = set(ChrList)
FIX_LENGTH = int(1e10)

##### system tools #####

TimeFormat = '%m/%d/%Y %H:%M:%S'

def curr_time():

    return  '[' + datetime.now().strftime(TimeFormat) + '] '

def myprint(string):

    process = psutil.Process(os.getpid())
    mem_usage = round(float(process.memory_info().rss))
    if mem_usage > 1e9:
        mem_usage = mem_usage/1e9
        mem_usage_str = ' (%3.3f GB)' % mem_usage
    elif mem_usage > 1e6:
        mem_usage = mem_usage/1e6
        mem_usage_str = ' (%3.3f MB)' % mem_usage
    else:
        mem_usage = mem_usage/1e3
        mem_usage_str = ' (%3.3f KB)' % mem_usage

    print >> sys.stderr, '[' + datetime.now().strftime(TimeFormat) +  mem_usage_str + ']', string

    return

##### list tools #####

def init_2d_list(n_row):

    a = [0] * n_row 
    for i in range(0, n_row):
        a[i] = list()

    return a

def int_list2string_list(int_list):

    string_list = list()
    for int_number in int_list:
        string_list.append(str(int_number))
    return string_list


##### file tools #####

def get_file_list(list_file):
    file_list = list()
    list_fp = open(list_file, 'r')
    lines = list(list_fp)
    for line in lines:
        line = line.strip()
        file_list.append(line)

    list_fp.close()
    return file_list

def check_file_exists(input_file):

    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
        return True
    else:
        return False

def line_count(in_file):
    n = 0
    in_file_fp = gzopen(in_file, 'r')
    while 1:
        line = in_file_fp.readline()
        if not line: break
        n += 1

    in_file_fp.close()
    return n

def gzopen(in_file, mode):

    if in_file[-2:] == 'gz':
        in_fp = gzip.open(in_file, mode)
    else:
        in_fp = open(in_file, mode)

    return in_fp

##### reference tools #####

def get_chrnames(faidx_file):
    fai_fp = open(faidx_file)
    tid2chrname_list = list()
    chrname2tid_dict = dict()
    tid = 0
    while 1:
        line = fai_fp.readline()
        if not line:
            break
        line = line.strip().split(tab)
        chrname = line[0]
        tid2chrname_list.append(chrname)
        chrname2tid_dict[chrname] = tid
        tid += 1

    return tid2chrname_list, chrname2tid_dict


def get_chr_length(faidx_file):
    chr_len_list = list()
    faidx_fp = open(faidx_file, 'r')
    line = faidx_fp.readline().strip().split('\t')
    tid = 0
    while len(line) > 1:
        chr_len_list.append(int(line[1]))
        line = faidx_fp.readline().strip().split('\t')
        tid += 1
    faidx_fp.close()
    return chr_len_list

def calculate_genome_length(faidx_file):

    genome_length = 0
    faidx_fp = open(faidx_file, 'r')
    while 1:
        line = faidx_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        chr_len = int(line[1])
        genome_length += chr_len

    return genome_length

def linear2normal(linear_pos_list, chr_len_list):
    n_chr = len(chr_len_list)
    cum_len_list = []
    cum_len_list.append(0)
    for i in range(1, n_chr):
        cum_len_list.append(chr_len_list[i-1] + cum_len_list[i-1])

    cum_len_list.append (sum(chr_len_list))  # cum_len_list has n_chr + 1 elements
        
    normal_pos_list = dict()   # normal_pos_list: dict(list())
    for i in range(0, n_chr):
        normal_pos_list[i] = list()

    for i in range(0, len(linear_pos_list)):
        for tid in range(0, n_chr): 
            if (linear_pos_list[i] >= cum_len_list[tid] and linear_pos_list[i] < cum_len_list[tid+1]):
                normal_pos = linear_pos_list[i] - cum_len_list[tid]
                normal_pos_list[tid].append(normal_pos) 

    return normal_pos_list

def normal2linear1pos (tid, pos, chr_len_list):
    for i in range(0, tid):
        pos += chr_len_list[i]
    return pos


def get_bcd_set(bcd11_file, pos_list, bin_size, win_size):
    '''
    this function returns the bcd set of (pos, pos + win_size * bin_size)
    '''
    # bcd11_file: string, pos_list: list

    bcd11_fp = open(bcd11_file, 'r')
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

    
    bcd_set = list()  # bcd_set: list(set())
    for i in range(0, len(pos_list)):  # for i-th position 
        start_pos   = pos_list[i] - pos_list[i] % bin_size
        start_index = bisect.bisect_left(bcd11_pos_list, start_pos)
        end_index = start_index + win_size
    
        bcd1set = set()

        for line in bcd11_lines[start_index:end_index]:
            line = line.strip().split(tab)
            for item in line[1:]:
                bcd1set.add(item.split(':')[0])

        bcd_set.append(bcd1set)

    bcd11_fp.close() 
    return bcd_set


def read_object_file(object_file, Classname, options = None):
    
    object_list = list()
    object_fp = open(object_file, 'r')
    while 1:
        line = object_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        if options == None: 
            obj = Classname(line)
        else:
            obj = Classname(line, options)

        object_list.append(obj)

    object_fp.close()

    return object_list

def read_alternative_contig_file(alternative_contig_file):
    alt_chr_name_set = set()
    in_fp = open(alternative_contig_file, 'r')
    while 1:
        line = in_fp.readline()
        if not line: break
        line = line.strip()
        alt_chr_name_set.add(line)

    in_fp.close()
    return alt_chr_name_set

def get_alternative_tid_set(alternative_contig_file, faidx_file):
    
    alt_chr_name_set = read_alternative_contig_file(alternative_contig_file)
    tid2chrname_list, chrname2tid_dict = get_chrnames(faidx_file) 
    alt_tid_set = chrname_set_2_tid_set(alt_chr_name_set, chrname2tid_dict)

    return alt_tid_set

def chrname_set_2_tid_set(chrname_set, chrname2tid_dict): 
    tid_set = set()
    for chrname in chrname_set:
        if chrname in chrname2tid_dict:
            tid = chrname2tid_dict[chrname]
            tid_set.add(tid)

    return tid_set

def make_dir(out_dir):

    try:
        os.makedirs(out_dir)
    except:
        if not os.path.isdir(out_dir): raise

    return
