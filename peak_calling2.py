#!/usr/bin/env python

import numpy as np
import peakutils
import bisect
from my_utils import *

class Interval:
    def __init__(self, tid, start, end, end5p_cnt, end3p_cnt):
        self.tid = int(tid)  # chr number , start from 0
        self.start = int(start) 
        self.end = int(end)
        self.end5p_cnt = int(end5p_cnt)
        self.end3p_cnt = int(end3p_cnt)

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    peak_calling2(args, dbo_args, endpoint_args)

    return


def peak_calling2(args, dbo_args, endpoint_args):

    is_wgs = args.is_wgs

    endpoints_density_file = endpoint_args.endpoints_density_file 

    #min_distance = endpoint_args.min_distance_for_peak_calling 
    min_distance = 20

    output_bk_candidate_file = endpoint_args.bk_file 

    min_frag_length = endpoint_args.min_frag_length 

    myprint ('reading endpoint density file: ' + endpoints_density_file)

    interval_list, end5p_density_list, end3p_density_list = get_endpoint_density_from_file(endpoints_density_file)

    myprint ('peak calling2...')

    myprint ('sorting endpoint density list')

    non_zero_end5p_density_list = list()
    for item in end5p_density_list:
        if item > 0: non_zero_end5p_density_list.append(item)
            
    non_zero_end3p_density_list = list()
    for item in end3p_density_list:
        if item > 0: non_zero_end3p_density_list.append(item)

    sorted_end5p_density_list = sorted(non_zero_end5p_density_list)
    sorted_end3p_density_list = sorted(non_zero_end3p_density_list)
    
    max_end5p_density = sorted_end5p_density_list[-1]
    max_end3p_density = sorted_end3p_density_list[-1]

    abs_end5_density_threshold = sorted_end5p_density_list[int(float(len(sorted_end5p_density_list)) * 0.95)]
    abs_end3_density_threshold = sorted_end3p_density_list[int(float(len(sorted_end3p_density_list)) * 0.95)]

    rel_end5_density_threshold = float(abs_end5_density_threshold) / max_end5p_density
    rel_end3_density_threshold = float(abs_end3_density_threshold) / max_end3p_density

    myprint ('abs_end5_density_threshold=%f, abs_end3_density_threshold=%f' % (abs_end5_density_threshold, abs_end3_density_threshold))
    myprint ('rel_end5_density_threshold=%f, rel_end3_density_threshold=%f' % (rel_end5_density_threshold, rel_end3_density_threshold))

    np_end5p_density_list = np.array(end5p_density_list)
    np_end3p_density_list = np.array(end3p_density_list)

    myprint ('peak calling for np_end5p_density_list')
    end5p_index_nparray = peakutils.indexes(np_end5p_density_list, thres=rel_end5_density_threshold, min_dist=min_distance)

    myprint ('peak calling for np_end3p_density_list')
    end3p_index_nparray = peakutils.indexes(np_end3p_density_list, thres=rel_end3_density_threshold, min_dist=min_distance)

    validated_end5p_index_list = end5p_index_nparray 
    validated_end3p_index_list = end3p_index_nparray

    '''
    myprint ('validating peaks...')
    for index in end5p_index_nparray: 
        if validate_peak_index(index, np_end5p_density_list, abs_end5_density_threshold, is_wgs) == True:
            validated_end5p_index_list.append(index)

    for index in end3p_index_nparray: 
        if validate_peak_index(index,np_end3p_density_list, abs_end3_density_threshold, is_wgs) == True:
            validated_end3p_index_list.append(index)

    myprint ('raw 5p peaks: ' + str(len(end5p_index_nparray)) + ', validated peaks: ' + str(len(validated_end5p_index_list)))
    myprint ('raw 3p peaks: ' + str(len(end3p_index_nparray)) + ', validated peaks: ' + str(len(validated_end3p_index_list)))
    '''

    myprint ('outputing results to file: ' + output_bk_candidate_file)
    validated_end5p_index_set = set(validated_end5p_index_list)
    validated_end3p_index_set = set(validated_end3p_index_list)
    all_validated_index_set = validated_end5p_index_set.union(validated_end3p_index_set)
    
    sorted_end5p_density_list_length = float(len(sorted_end5p_density_list))
    sorted_end3p_density_list_length = float(len(sorted_end3p_density_list))

    output_bk_candidate_fp = open(output_bk_candidate_file, 'w') 
    for i in range(0, len(interval_list)):
        if i not in all_validated_index_set: continue
        interval1 = interval_list[i]

        end5p_density = end5p_density_list[i]
        end3p_density = end3p_density_list[i]

        quantile5p_index = bisect.bisect(sorted_end5p_density_list, end5p_density)-1 
        quantile5p = float(sorted_end5p_density_list_length - quantile5p_index) / sorted_end5p_density_list_length
        if quantile5p < 1e-100: quantile5p = 1e-100
        mlogquantile5p = -math.log(quantile5p, 10)

        quantile3p_index = bisect.bisect(sorted_end3p_density_list, end3p_density)-1 
        quantile3p = float(sorted_end3p_density_list_length - quantile3p_index) / sorted_end3p_density_list_length
        if quantile3p < 1e-100: quantile3p = 1e-100
        mlogquantile3p = -math.log(quantile3p, 10)

        if i in validated_end5p_index_set:
            endtype = '5p_end'
            bkstart = int(round(interval1.start + interval1.end + 1)/2.0) - 50
            bkend   = bkstart + 100
            mlogquantile = mlogquantile5p
        else:
            endtype = '3p_end'
            bkstart = int(round(interval1.start + interval1.end + 1)/2.0) - 50
            bkend   = bkstart + 100
            mlogquantile = mlogquantile3p

        info = 'endtype=' + endtype + ';window_start=' + str(interval1.start) + ';window_end=' + str(interval1.end) + ';5p_ends_cnt=' + str(interval1.end5p_cnt) + ';3p_ends_cnt=' + str(interval1.end3p_cnt) + ';5p_ends_score=' + '%.3f'%mlogquantile5p + ';3p_ends_score=' + '%.3f'%mlogquantile3p

        output_bk_candidate_fp.write(str(interval1.tid) + tab + str(bkstart) + tab + str(bkend) + tab + '%.3f'%mlogquantile + tab + info + endl)

    myprint ('finished peak calling for endpoint density file: ' + endpoints_density_file)
    
    return

def validate_peak_index(index, np_density_list, abs_threshold, is_wgs):

    return True
    '''
    win_size = 500
    np_density_list_length = len(np_density_list)    
    if win_size > np_density_list_length:
        win_size = np_density_list_length

    if index < win_size:
        start = 0
        end = start + win_size
    elif index + win_size >= np_density_list_length:
        end = np_density_list_length
        start = end - win_size
    else:
        start = index - win_size/2
        end = start + win_size

    array = np_density_list[start:end]
    mean = np.mean(array)
    #std = np.std(array)

    del array
    
    if np_density_list[index] > abs_threshold and np_density_list[index] > mean:
        return True
    else:
        return False
    '''

def get_endpoint_density_from_file(endpoints_density_file):
    
    interval_list = list()
    end5p_density_list = list()
    end3p_density_list = list()
    endpoints_density_fp = open(endpoints_density_file, 'r')

    while 1:
        line = endpoints_density_fp.readline()
        if not line:    
            break
        line = line.strip().split(tab)
        interval1 = Interval(line[0], line[1], line[2], line[4], line[5])
        interval_list.append(interval1)
        end5p_density_list.append(float(line[4])) # end5p_cnt
        end3p_density_list.append(float(line[5])) # end3p_cnt
    
    endpoints_density_fp.close()

    return interval_list, end5p_density_list, end3p_density_list

if __name__ == '__main__':
    main()
