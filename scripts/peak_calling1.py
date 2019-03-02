#!/usr/bin/env python

import numpy as np
import peakutils
from scipy.stats import norm
from my_utils import *


def main():
    
    args, dbo_args, endpoint_args = parse_user_arguments()
    peak_calling(bcd13_list, min_distance, output_bk_candidate_file)

    return

def peak_calling1(args, dbo_args, endpoint_args):

    bcd13_list = dbo_args.bcd13_file_list
    min_distance = dbo_args.min_distance_for_peak_calling
    output_bk_candidate_file = dbo_args.bk_file
    
    output_bk_candidate_fp = open(output_bk_candidate_file, 'w')
    for tid in range(0, len(bcd13_list)):
        bcd13_file = bcd13_list[tid]
        peak_calling_bcd13_file(tid, bcd13_file, min_distance, output_bk_candidate_fp)
    output_bk_candidate_fp.close()
    return

def read_bcd13_file(bcd13_file):
    
    myprint ('start reading bcd13 files')
    bcd13_line_list = list()
    xlist = list()
    ylist = list()
    effective_ylist = list() # y values where m1 > 50 and m2 > 50

    bcd13_fp = open(bcd13_file, 'r')
    while 1:
        bcd13_line = bcd13_fp.readline()
        if not bcd13_line:
            break
        bcd13_line_list.append(bcd13_line)
        line = bcd13_line.strip().split(tab)
        m1 = int(line[1])
        m2 = int(line[2])
        m0 = int(line[3])
        xlist.append(int(line[0]))
        ylist.append(float(line[7]))
        if m1 > 50 and m2 > 50:
            effective_ylist.append(float(line[7]))

    bcd13_fp.close()
    myprint ('finished reading bcd13 files')
    return bcd13_line_list, xlist, ylist, effective_ylist

def output(np_log_effective_y):
    outfile = 'logy_values.txt'
    outfp = open(outfile, 'w')
    for logy in np_log_effective_y:
        outfp.write(str(logy) + endl)
    outfp.close()

def peak_calling_bcd13_file(tid, bcd13_file, min_distance, output_bk_candidate_fp):
    myprint('peak calling for: ' + bcd13_file)
    bcd13_line_list, xlist, ylist, effective_ylist = read_bcd13_file(bcd13_file)
    npy = np.array(ylist) 
    max_y_avlue = max(npy)
    
    np_log_effective_y = np.log(effective_ylist)


    logymean = np.mean(np_log_effective_y)
    logystd  = np.std(np_log_effective_y)
    myprint ('number of effective y values: %d' %(len(np_log_effective_y)))

    abs_threshold1 = np.percentile(npy, 95)
    abs_threshold2 = 1.2 
    abs_threshold  = max(abs_threshold1, abs_threshold2)

    for q in range(0, 101, 5):
        myprint ('%d-th percentile:\t%.2f' % (q, np.percentile(npy, q))) 

    rel_threshold = float(abs_threshold)/ max_y_avlue
    myprint ('abs threshold: %.2f,  max y value: %.2f' % (abs_threshold, max_y_avlue))
    myprint ('relative threshold: %.2f' % (rel_threshold))

    index_list = peakutils.indexes(npy, thres=rel_threshold, min_dist=min_distance)
    validated_index_list = list()
    for index in index_list:
        if validate_peak_index(index, ylist) == True:
            validated_index_list.append(index)

    myprint ('number of raw peaks: ' + str(len(index_list)) + ', number of validated peak: ' + str(len(validated_index_list)))

    max_mlogpvalue = 10.0
    for i in validated_index_list:
        logy = math.log(ylist[i]) 
        pvalue = 1.0 - norm.cdf(logy, logymean, logystd)
        if pvalue > 0:
            mlogpvalue = -math.log10(pvalue)
        else:
            mlogpvalue = max_mlogpvalue

        if mlogpvalue > max_mlogpvalue:
            mlogpvalue = max_mlogpvalue

        line = bcd13_line_list[i]
        pos, m1, m2, m0, d1, d2, predict_m0, ratio = line.strip().split(tab)
        if int(m1) < 100 or int(m2) < 100:
            continue
        info = 'm1=' + m1 + ';m2=' + m2 + ';m0=' + m0 + ';expected_observed_ratio=' + ratio
        output_bk_candidate_fp.write(str(tid) + tab + str(xlist[i]) + tab + str(xlist[i+1]) + tab + '%.3f'%mlogpvalue + tab + info + endl)

    return

def validate_peak_index(index, y):

    win_size = 100

    if win_size > len(y):
        win_size = len(y)

    if index < win_size:
        start = 0
        end = start + win_size
    elif index + win_size >= len(y):
        end = len(y)
        start = end - win_size
    else:
        start = index - win_size/2
        end = start + win_size

    array = y[start:end]
    mean = np.mean(array)
    std = np.std(array)
    
    if y[index] > mean: 
        return True
    else:
        return False

if __name__ == '__main__':
    main()

