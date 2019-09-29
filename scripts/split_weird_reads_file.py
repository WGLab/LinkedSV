#!/usr/bin/env python 

import os
import sys
import gzip

from my_utils import *

def main():

    if len(sys.argv) < 4: 
        print ('at least 3 arguments')
        exit(1)

    bcd21_file = sys.argv[1]
    weird_reads_file = sys.argv[2]
    num_split_bcd21_file = sys.argv[3]

    n_bcd_per_file = 50000

    prefix = os.path.splitext(bcd21_file)[0]

    myprint ('extracting weird reads names')
    readname_fileid_dict = get_weird_reads_names(weird_reads_file)
    myprint ('number of weird reads: %d' % len(readname_fileid_dict) )

    myprint ('processing bcd21 file')
    n_existed_bcd21_split_file = split_bcd21_file(bcd21_file, n_bcd_per_file, readname_fileid_dict)

    myprint ('processing weird reads file')
    split_weird_reads_file (weird_reads_file, readname_fileid_dict, n_existed_bcd21_split_file)

    num_split_bcd21_fp = open (num_split_bcd21_file, 'w')
    num_split_bcd21_fp.write ('%d\n' % n_existed_bcd21_split_file)
    num_split_bcd21_fp.close ()

    del readname_fileid_dict
    gc.collect()
    myprint ('finished processing weird reads file')

    return

def get_weird_reads_names(weird_reads_file):
    readname_fileid_dict = dict()
    weird_reads_fp = open(weird_reads_file, 'r')
    while 1:
        line = weird_reads_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line2 = line.strip().split(tab)
        readname = line2[8].split('@')[0]
        readname_fileid_dict[readname] = 'N'

    weird_reads_fp.close()
    return readname_fileid_dict

def split_weird_reads_file (weird_reads_file, readname_fileid_dict, n_existed_bcd21_split_file):

    out_fp_list = [0] * (n_existed_bcd21_split_file +1)

    for i in range(1, len(out_fp_list)):
        split_weird_reads_file = weird_reads_file + '.split%d' % (i)
        out_fp_list[i] = open(split_weird_reads_file, 'w')

    weird_reads_fp = open(weird_reads_file, 'r')
    while 1:
        line = weird_reads_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line2 = line.strip().split(tab)
        readname = line2[8].split('@')[0]
        if readname not in readname_fileid_dict: 
            continue

        file_id = readname_fileid_dict[readname]
        if file_id == 'N': 
            continue

        out_fp_list[file_id].write(line)
    
    for i in range(1, len(out_fp_list)):
        out_fp_list[i].close()

    weird_reads_fp.close()

    return    

def split_bcd21_file(bcd21_file, n_bcd_per_file, readname_fileid_dict):

    file_id = 1 
    bcd_id = 0
    curr_bcd = '' 

    split_barcode_list = list()

    if bcd21_file[-2:] == 'gz':
        bcd21_fp = gzip.open(bcd21_file, 'rt')
    else:
        bcd21_fp = open(bcd21_file, 'r')

    bcd21_split_file = bcd21_file + '.split'
    bcd21_split_fp = open(bcd21_split_file, 'w')
    while 1:
        line = bcd21_fp.readline()
        if not line: break
        if line[0] == '#': continue
        line2 = line.strip().split(tab)
        bcd = line2[4]
        readname = line2[6]
        if bcd != curr_bcd:
            curr_bcd = bcd 
            bcd_id += 1

        if bcd_id <= n_bcd_per_file:
            if readname in readname_fileid_dict: readname_fileid_dict[readname] = file_id
        else:
            file_id += 1
            bcd21_split_fp.write('%s\t%d\n' % (bcd, file_id))
            if readname in readname_fileid_dict: readname_fileid_dict[readname] = file_id
            bcd_id = 1
            
    bcd21_split_fp.close()
    bcd21_fp.close()

    return file_id

if __name__ == '__main__':
    main()
