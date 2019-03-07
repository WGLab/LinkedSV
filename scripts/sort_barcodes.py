#!/usr/bin/env python

import gc
from my_utils import *

tab = '\t'
endl = '\n'

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()
    sort_barcodes(args, dbo_args, endpoint_args)

    return

class BcdInfo:
    def __init__(self, tid, start, end, mapq, bcd, other_info):

        self.tid = int(tid)
        self.start = int(start)
        self.end = int(end)
        self.mapq = int(mapq)
        self.bcd = bcd
        self.other_info = other_info

    def output(self):
        outstring = '%d\t%d\t%d\t%d\t%s\t%s' % (self.tid, self.start, self.end, self.mapq, self.bcd, self.other_info)  
        return outstring


def sort_barcodes(args, dbo_args, endpoint_args):

    ## first round ##
    bcd_file_prefix = args.out_prefix + '.bcd'
    sort_bcd_file = endpoint_args.bcd21_file
    sort_bcd_fp = open(sort_bcd_file, 'w')
    for i in range(0, 1024):
        bcd_file  = bcd_file_prefix + '.part%d' % i
        sort_barcodes_for1file(args, dbo_args, endpoint_args, bcd_file, sort_bcd_fp, i)
        gc.collect()
        #os.remove(bcd_file)

    sort_bcd_fp.close()
    return

def sort_barcodes_for1file(args, dbo_args, endpoint_args, bcd_file, sort_bcd_fp, fid):

    myprint('start sorting barcodes for file: %s' % bcd_file) 
    headerline, bcdinfo_list = read_bcd_file(bcd_file)
    if fid == 0: sort_bcd_fp.write(headerline)

    sorted_bcdinfo_list = sorted(bcdinfo_list, key = lambda bcdinfo:(bcdinfo.bcd, bcdinfo.tid, bcdinfo.start))

    for bcdinfo in sorted_bcdinfo_list:
        sort_bcd_fp.write(bcdinfo.output() + endl)

    myprint('finished sorting barcodes for file: %s' % bcd_file) 
    return

def read_bcd_file(bcd_file):

    bcd_fp = open(bcd_file, 'r')
    headerline = ''
    bcdinfo_list = list()
    while 1:
        line = bcd_fp.readline()
        if not line: break
        if line[0] == '#': 
            headerline += line
            continue
        line = line.strip().split(tab)
        tid = line[0]
        start = line[1]
        end = line[2]
        mapq = line[3]
        bcd = line[4]
        other_info = tab.join(line[5:])
        bcdinfo = BcdInfo(tid, start, end, mapq, bcd, other_info)
        bcdinfo_list.append(bcdinfo)

    bcd_fp.close()
    return headerline, bcdinfo_list 

if __name__ == '__main__':
    main()
