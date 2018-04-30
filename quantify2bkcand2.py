#!/usr/bin/env python

import math
import numpy as np
from my_utils import *
from fragment import *
from bed import *
from bedpe import *
import bisect
from scipy import spatial
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import poisson
from scipy.stats import gamma 
from scipy.stats import norm 
from scipy.stats import expon
from scipy import stats
import global_distribution

def main():
    
    args, dbo_args, endpoint_args = parse_user_arguments()
    quantify2bkcand(args, dbo_args, endpoint_args)
    return 


def quantify2bkcand(args, dbo_args, endpoint_args):

    myprint('estimation global parameters')
    if args.global_distribution_calculated == False:
        global_distribution.estimate_global_distribution(args, dbo_args, endpoint_args)

    myprint ('reading paired breakpoint candidate file: %s' % args.bk_cand_pair_file)
    paired_bk_cand_list = read_paired_bk_cand_file(args.bk_cand_pair_file)

    myprint ('reading bcd22 file: %s' % endpoint_args.bcd22_file)
    if args.is_wgs:
        min_frm_length = endpoint_args.min_frag_length/2
    else:
        min_frm_length = 0 

    bcd22_frm_list = read_bcd22_file(endpoint_args.bcd22_file, min_frm_length)

    if os.path.exists(args.node33_file):
        args.n_node33 = line_count(args.node33_file)
    else:
        myprint ('ERROR! node33 file does not exist!')
        sys.exit()

    myprint ('sorting fragments by start position')
    start_sorted_frm_list = sorted(bcd22_frm_list, key = lambda frm: frm.key_start())
    start_key_list = list()
    for frm in start_sorted_frm_list: start_key_list.append(frm.key_start())

    myprint ('sorting fragments by end position')
    end_sorted_frm_list = sorted(bcd22_frm_list, key = lambda frm: frm.key_end())
    end_key_list = list()
    for frm in end_sorted_frm_list: end_key_list.append(frm.key_end())

    myprint ('quantifying candidate calls')

    if args.is_wgs == False: # for targeted sequencing
        target_region_bed_db = dict()
        target_region_bed_list = list()
        target_region_bed_list = read_bed_file(args.target_region_bed, args.chrname2tid)
        for bed in target_region_bed_list:
            if bed.tid not in target_region_bed_db: target_region_bed_db[bed.tid] = list() 
            target_region_bed_db[bed.tid].append(bed)
        for tid in target_region_bed_db:
            target_region_bed_db[tid].sort(key = lambda bed: bed.start) 

        target_region_bed_startpos_db = dict() 
        for tid in target_region_bed_db:
            target_region_bed_startpos_db[tid] = list()
            for bed in target_region_bed_db[tid]: target_region_bed_startpos_db[tid].append(bed.start)

        args.target_region_bed_db = target_region_bed_db
        args.target_region_bed_startpos_db = target_region_bed_startpos_db


    quantified_bk_cand_list = list()
    i = 0
    for i in range(0, len(paired_bk_cand_list)):
        paired_bk_cand = paired_bk_cand_list[i]
        quantified_bk_cand = quantify1paired_bk_cand(args, dbo_args, endpoint_args, paired_bk_cand, bcd22_frm_list, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list)
        if (i+1) % 1000 == 0: myprint('quantified %d candidate calls' % (i+1))
        if quantified_bk_cand == None: continue
        quantified_bk_cand_list.append(quantified_bk_cand)


    myprint('finished quantifying candidate calls, outputing results to file: %s' % args.quantified_bk_pair_file)
    out_file = args.quantified_bk_pair_file
    out_fp = open(out_file, 'w')

    for quantified_bk_cand in quantified_bk_cand_list: out_fp.write(quantified_bk_cand.output() + endl) 
    out_fp.close()

    del bcd22_frm_list, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list
    return 

def split_same_fragment(args, same_fragment_list):

    new_gap_distance_cutoff = args.gap_distance950
    split_same_fragment_list1 = list()
    split_same_fragment_list2 = list()
    non_split_same_fragment_list = list()

    for frm in same_fragment_list:
        if frm.num_reads < 3:
            non_split_same_fragment_list.append(frm)
            continue
        map_pos = frm.map_pos 
        map_pos = map_pos.strip(';').split(';')
        read_start_list = list()
        read_end_list = list()
        for item in map_pos:
            read_start = int(item.split(',')[0])
            read_start_list.append(read_start)
            read_end = int(item.split(',')[1])
            read_end_list.append(read_end)
        max_gap_distance = 0
        max_gap_index = -1
        for i in range(1, len(read_start_list)):
            if max_gap_distance < read_start_list[i] - read_start_list[i-1]:
                max_gap_distance = read_start_list[i] - read_start_list[i-1]
                max_gap_index = i

        if max_gap_distance > new_gap_distance_cutoff:

            new_start1 = read_start_list[0]
            new_end1 = read_end_list[max_gap_index-1]
            new_start2 = read_start_list[max_gap_index]
            new_end2 = read_end_list[-1]
            new_num_reads1 = max_gap_index
            new_num_reads2 = frm.num_reads - new_num_reads1
            new_map_pos1 = ';'.join(map_pos[0:max_gap_index]) + ';'
            new_map_pos2 = ';'.join(map_pos[max_gap_index:]) + ';'
            new_frag_id1 = frm.frag_id + int(1e12)
            new_frag_id2 = frm.frag_id + int(2e12)

            new_attr_list1 = [frm.tid, new_start1, new_end1, new_end1-new_start1, frm.bcd, new_frag_id1, new_num_reads1, frm.hp0, frm.hp1, frm.hp2, new_map_pos1]
            new_frm1 = Fragment(new_attr_list1)

            new_attr_list2 = [frm.tid, new_start2, new_end2, new_end2-new_start2, frm.bcd, new_frag_id2, new_num_reads2, frm.hp0, frm.hp1, frm.hp2, new_map_pos2]
            new_frm2 = Fragment(new_attr_list2)
            split_same_fragment_list1.append(new_frm1)
            split_same_fragment_list2.append(new_frm2)
        else:
            non_split_same_fragment_list.append(frm)

    return split_same_fragment_list1, split_same_fragment_list2, non_split_same_fragment_list 


def get_region_frm_list(tid, start, end, endtype, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list):

    frm_list = list()
    ## get all fragments of which the endpoints in the region
    if endtype == '5p_end':
        searchstart_index = bisect.bisect_left(start_key_list, tid * FIX_LENGTH + start) 
        searchend_index   = bisect.bisect_right(start_key_list, tid * FIX_LENGTH + end)
        for i in range(searchstart_index, searchend_index):
            frm_list.append(start_sorted_frm_list[i])
    elif endtype == '3p_end': 
        searchstart_index = bisect.bisect_left(end_key_list, tid * FIX_LENGTH + start)
        searchend_index   = bisect.bisect_right(end_key_list,  tid * FIX_LENGTH + end)
        for i in range(searchstart_index, searchend_index):
            frm_list.append(end_sorted_frm_list[i])
    else:
        myprint ('ERROR! unknown endtype: %s' % endtype)
        sys.exit()

    return frm_list

def quantify1paired_bk_cand(args, dbo_args, endpoint_args, paired_bk_cand, bcd22_frm_list, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list):
    

    paired_bk_cand.format_self(args.chrname2tid)
    tid1 = paired_bk_cand.tid1(args.chrname2tid)
    tid2 = paired_bk_cand.tid2(args.chrname2tid)

    search_range = args.gap_distance950
    if search_range == None: search_range = paired_bk_cand.end1 - paired_bk_cand.start1

    if paired_bk_cand.endtype1 == '5p_end':
        start1 = paired_bk_cand.start1 
        end1 = start1 + search_range 
    else:
        end1 = paired_bk_cand.end1
        start1 = end1 - search_range

    if paired_bk_cand.endtype2 == '5p_end':
        start2 = paired_bk_cand.start2 
        end2 = start2 + search_range 
    else:
        end2 = paired_bk_cand.end2
        start2 = end2 - search_range

    frm_list1 = get_region_frm_list(tid1, paired_bk_cand.start1, paired_bk_cand.end1, paired_bk_cand.endtype1, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list)
    frm_list2 = get_region_frm_list(tid2, paired_bk_cand.start2, paired_bk_cand.end2, paired_bk_cand.endtype2, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list)

    shared_fragment_list1, shared_fragment_list2, same_fragment_list = get_shared_frag_list(frm_list1, frm_list2)

    num_same_frm = len(same_fragment_list)
    num_shared_bcd_2mol = len(shared_fragment_list1)
    num_shared_bcd = num_same_frm + num_shared_bcd_2mol


    split_same_fragment_list1, split_same_fragment_list2, non_split_same_fragment_list = split_same_fragment(args, same_fragment_list)

    num_split_mol = len(split_same_fragment_list1)

    if num_shared_bcd_2mol + num_split_mol < args.min_support_fragments: return None

    shared_fragment_list1 += split_same_fragment_list1
    shared_fragment_list2 += split_same_fragment_list2

    shared_fragment_ids1 = ''
    shared_fragment_ids2 = ''
    support_barcodes = ''
    for frm in shared_fragment_list1:
        shared_fragment_ids1 += '%d,' % frm.frag_id
        support_barcodes += frm.bcd + ','
    for frm in shared_fragment_list2:
        shared_fragment_ids2 += '%d,' % frm.frag_id
    
    shared_fragment_ids1.strip(',')
    shared_fragment_ids2.strip(',')
    support_barcodes.strip(',')

    endtype1 = paired_bk_cand.endtype1
    endtype2 = paired_bk_cand.endtype2

    total_logp_gap_ontarget, total_logp_gap_offtarget, total_logp_frm_length = logp_nosv_one_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2)
    logp_barcode = logp_nosv_two_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2)

    llr_barcode_overlapping = min(max(total_logp_gap_ontarget, total_logp_frm_length), logp_barcode)
    llg1 = total_logp_gap_ontarget 
    llg2 = total_logp_gap_offtarget
    llg3 = total_logp_frm_length 
    llg4 = logp_barcode

    endtype1_logp, endtype2_logp, predicted_endtype1, predicted_endtype2, predicted_svtype, start1_logp, end1_logp, start2_logp, end2_logp, bk1_pos, bk2_pos = quantification_svtype(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, paired_bk_cand)

    svtype_logp = min(endtype1_logp, endtype2_logp)

    tid1 = shared_fragment_list1[0].tid
    tid2 = shared_fragment_list2[0].tid
    chrm1 = args.tid2chrname[tid1]
    chrm2 = args.tid2chrname[tid2]

    if tid1 == tid2:
        svlength = str(abs(bk1_pos - bk2_pos))
    else:
        svlength = 'N.A.'

    attr_list = [chrm1, bk1_pos, bk1_pos+1, chrm2, bk2_pos, bk2_pos+1, predicted_svtype, svlength, len(shared_fragment_list1), predicted_endtype1, predicted_endtype2, llr_barcode_overlapping, llg1, llg2, llg3, llg4, svtype_logp, endtype1_logp, endtype2_logp, start1_logp, end1_logp, start2_logp, end2_logp, shared_fragment_ids1, shared_fragment_ids2, support_barcodes, tid1, tid2]

    quantified_bk_cand = QuantifiedBKCand(attr_list)
    
    return quantified_bk_cand

    
def quantification_svtype(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, paired_bk_cand):
    
    endtype1 = paired_bk_cand.endtype1
    endtype2 = paired_bk_cand.endtype2

    start_pos_list1 = list() 
    start_pos_list2 = list()
    end_pos_list1 = list()
    end_pos_list2 = list()

    for i in range(0, len(shared_fragment_list1)):
        start_pos_list1.append(shared_fragment_list1[i].start)
        end_pos_list1.append(shared_fragment_list1[i].end) 

        start_pos_list2.append(shared_fragment_list2[i].start)
        end_pos_list2.append(shared_fragment_list2[i].end) 

    start1_logp = estimate_breakpoint_logp(args, start_pos_list1)
    end1_logp   = estimate_breakpoint_logp(args, end_pos_list1)

    start2_logp = estimate_breakpoint_logp(args, start_pos_list2)
    end2_logp   = estimate_breakpoint_logp(args, end_pos_list2)

    if endtype1 == '5p_end': 
        endtype1_logp = start1_logp - end1_logp   
    else:
        endtype1_logp = end1_logp - start1_logp  

    if endtype2 == '5p_end': 
        endtype2_logp = start2_logp - end2_logp   
    else:
        endtype2_logp = end2_logp - start2_logp  

    predicted_endtype1 = endtype1
    predicted_endtype2 = endtype2

    if predicted_endtype1 == '5p_end' and predicted_endtype2 == '5p_end':
        predicted_svtype = 'INV' 
    elif predicted_endtype1 == '3p_end' and predicted_endtype2 == '3p_end':
        predicted_svtype = 'INV' 
    elif predicted_endtype1 == '5p_end' and predicted_endtype2 == '3p_end':
        if np.median(start_pos_list1) > np.median(end_pos_list2):
            predicted_svtype = 'DEL'
        else: 
            predicted_svtype = 'DUP'
    elif predicted_endtype1 == '3p_end' and predicted_endtype2 == '5p_end':
        if np.median(end_pos_list1) > np.median(start_pos_list2):
            predicted_svtype = 'DUP'
        else: 
            predicted_svtype = 'DEL'
    else:
        predicted_svtype = 'UNK'

    if paired_bk_cand.chrm1 != paired_bk_cand.chrm2: predicted_svtype = 'TRA'

    bk1_pos, bk2_pos = predict_bk_pos(start_pos_list1, start_pos_list2, end_pos_list1, end_pos_list2, predicted_endtype1, predicted_endtype2, args.gap_distance500)

    return endtype1_logp, endtype2_logp, predicted_endtype1, predicted_endtype2, predicted_svtype, start1_logp, end1_logp, start2_logp, end2_logp, bk1_pos, bk2_pos 


def predict_bk_pos(start_pos_list1, start_pos_list2, end_pos_list1, end_pos_list2, predicted_endtype1, predicted_endtype2, median_gap_distance):

    n_supp = len(start_pos_list1)
    bk_shift = int(median_gap_distance)

    if predicted_endtype1 == '5p_end':
        bk1_pos = int(np.median(start_pos_list1)) - bk_shift
    else:
        bk1_pos = int(np.median(end_pos_list1)) + bk_shift

    if predicted_endtype2 == '5p_end': 
        bk2_pos = int(np.median(start_pos_list2)) - bk_shift
    else:
        bk2_pos = int(np.median(end_pos_list2)) + bk_shift

    return bk1_pos, bk2_pos


def estimate_breakpoint_logp(args, position_list):

    p_endpoint_nonbk = float(args.gap_distance750) / float(args.median_fragment_length)
    p_endpoint_bk = 0.75

    median_pos = np.median(position_list)     
    half_gap_distance750 = args.gap_distance750 / 2.0
    num_within_half_gap_distance750 = 0
    for i in range(0, len(position_list)):
        if abs(position_list[i]-median_pos) <= half_gap_distance750:
            num_within_half_gap_distance750 += 1
   
    num_outside_half_gap_distance750 = len(position_list) - num_within_half_gap_distance750

    if p_endpoint_bk < 1e-100: p_endpoint_bk = 1e-100
    if p_endpoint_nonbk < 1e-100: p_endpoint_nonbk = 1e-100
    if p_endpoint_nonbk > p_endpoint_bk: 
        p_endpoint_nonbk = p_endpoint_bk
        myprint('WARNING: median fragment length is too small or gap distance is too large')

    logp_bk = num_within_half_gap_distance750 * math.log(p_endpoint_bk, 10) + num_outside_half_gap_distance750 * math.log(1.0-p_endpoint_bk, 10) 
    logp_nonbk = num_within_half_gap_distance750 * math.log(p_endpoint_nonbk, 10) + num_outside_half_gap_distance750 * math.log(1.0-p_endpoint_nonbk, 10) 

    logp = logp_bk - logp_nonbk

    return logp

def add2logvalue(a, b):
    if a < b: 
        tmp = a
        a = b
        b = tmp

    return a + math.log(1 + pow(10, (b-a)))

def logp_sv_two_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2):
    logp = 0
    for i in range(0, len(shared_fragment_list1)):
        frm1 = shared_fragment_list1[i]
        frm2 = shared_fragment_list2[i]
        pbarcode  = float(args.gap_distance_cutoff) * int(args.median_num_fragment_per_bcd+1) / args.genome_length
        if pbarcode < 1e-100: pbarcode = 1e-100
        logp_barcode = math.log(pbarcode, 10)
        lop_fragment_length1 = logp_fragment_length(args, frm1.length) 
        lop_fragment_length2 = logp_fragment_length(args, frm2.length) 
        logp += logp_barcode + lop_fragment_length1 + lop_fragment_length2
    
    return logp

def logp_sv_one_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2):
    logp = 0
    for i in range(0, len(shared_fragment_list1)):
        frm1 = shared_fragment_list1[i]
        frm2 = shared_fragment_list2[i]
        logp_barcode = 0
        total_frm_length = frm1.length + frm2.length 
        logp_total_frm_length = logp_fragment_length (args, total_frm_length)
        logp += logp_barcode + logp_total_frm_length

    return logp

def logp_nosv_two_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2):

    area_fold = float(args.genome_length) / float(args.gap_distance_cutoff) 
    area_fold = area_fold * area_fold
    lmda = args.n_node33 / area_fold
    n_supp = len(shared_fragment_list1)
    logpmf = poisson.logpmf(n_supp, lmda) / math.log(10)
    logp = logpmf + math.log(area_fold, 10) 
    return -logp 
       
def logp_nosv_one_mol(args, dbo_args, endpoint_args, shared_fragment_list1, shared_fragment_list2, endtype1, endtype2):

    read_per_bp_genome = args.read_per_bp_genome 
    assumed_frm_length_list = list()
    p_gap_ontarget_list = list()
    p_gap_offtarget_list = list()
    for i in range(0, len(shared_fragment_list1)):
        frm1 = shared_fragment_list1[i]
        frm2 = shared_fragment_list2[i]
        gap_length, gap_tid, gap_start, gap_end = fragment_min_distance(frm1, frm2)
        if gap_length <= 0:
            logp_gap_ontarget = 0
            logp_gap_offtarget = 0
        else:
            if args.is_wgs:
                logp_gap_ontarget = (gap_length + 1) * math.log( (1.0 - args.read_per_bp_genome), 10)
                logp_gap_offtarget = 0
            else:
                on_target_length = get_ontarget_length(args, gap_tid, gap_start, gap_end)
                off_target_length = gap_length - on_target_length
                logp_gap_ontarget = (on_target_length + 1) * math.log((1.0 - args.read_per_bp_ontarget), 10)
                logp_gap_offtarget = (off_target_length + 1) * math.log((1.0 - args.read_per_bp_offtarget), 10)

            p_gap_ontarget_list.append (math.exp(logp_gap_ontarget))
            p_gap_offtarget_list.append (math.exp(logp_gap_offtarget))

        assumed_frm_length = frm1.length + frm2.length + gap_length 
        assumed_frm_length_list.append (assumed_frm_length)

    total_logp_frm_length    = cal_logp_fragment_length (args, assumed_frm_length_list)
    total_logp_gap_ontarget  = cal_logp_gap_distance (args, p_gap_ontarget_list) 
    total_logp_gap_offtarget = cal_logp_gap_distance (args, p_gap_offtarget_list) 

    return -total_logp_gap_ontarget, -total_logp_gap_offtarget, -total_logp_frm_length

def cal_logp_gap_distance(args, pvalue_list):

    if len(pvalue_list) == 0: return 0
    n_supp = len(pvalue_list)
    print pvalue_list
    sem = math.sqrt(1.0/ (12.0 * n_supp))
    pvalue_mean = float(np.mean(pvalue_list))
    z = (pvalue_mean - 0.5) / sem 
    logp = norm.logcdf(z) / math.log(10)
    if logp == -float('inf'): logp = -10000000
    n_compare_times = args.num_reads_genome / n_supp
    logp += math.log(n_compare_times, 10)
    print n_supp, sem, pvalue_mean, z, logp, n_compare_times
    return logp

def cal_logp_fragment_length(args, assumed_frm_length_list):

    n_supp = len(assumed_frm_length_list)
    mean_supp_fragments = np.mean(assumed_frm_length_list)
    mean_population = 1.0 / args.fragment_length_lmda
    logp = gamma.logsf(mean_supp_fragments * n_supp, a=n_supp, scale=mean_population) / math.log(10.0)
    if logp == -float('inf'): logp = -10000000
    n_compare_times = args.total_num_fragment / n_supp
    logp += math.log(n_compare_times, 10)
    return logp 

def cal_logp_fragment_length2(args, assumed_frm_length_list):
    n_supp = len(assumed_frm_length_list)
    mean_supp_fragments = np.mean(assumed_frm_length_list)
    population_mean = 1.0 / args.fragment_length_lmda
    p = expon.sf(mean_supp_fragments, scale=population_mean)
    lmda = n_supp * p
    log_pmf = poisson.logpmf(n_supp, lmda)  / math.log(10)
    if log_pmf == -float('inf'): log_pmf = -10000000
    log_expected = log_pmf + math.log(args.total_num_fragment / n_supp, 10)
    return log_expected
     
def get_ontarget_length(args, gap_tid, gap_start, gap_end):

    if gap_tid not in args.target_region_bed_db:
        return  0
    on_target_length = 0
    
    target_region_bed_db = args.target_region_bed_db
    index = bisect.bisect(args.target_region_bed_startpos_db[gap_tid], gap_start)
    index = index - 2
    if index < 0: index = 0
    for i in range(index, len(target_region_bed_db[gap_tid])):
        bed = target_region_bed_db[gap_tid][i]
        overlap_length = min(gap_end, bed.end) - max(gap_start, bed.start) 
        if overlap_length > 0: on_target_length += overlap_length
        if bed.start > gap_end: break 
    
    return on_target_length

def fragment_min_distance(frm1, frm2): 
    if frm1.tid != frm2.tid:
        return FIX_LENGTH, None, None, None

    maxstart = max(frm1.start, frm2.start)
    minend = min(frm1.end, frm2.end)
    distance = maxstart - minend

    return distance, frm1.tid, minend, maxstart 

def get_shared_frag_list(fragment_list1, fragment_list2):

    frm_id_set1 = set()
    frm_id_set2 = set() 
    frm_bcd_set1 = set()
    frm_bcd_set2 = set() 

    for frm in fragment_list1:
        frm_id_set1.add(frm.frag_id)
        frm_bcd_set1.add(frm.bcd)

    for frm in fragment_list2:
        frm_id_set2.add(frm.frag_id)
        frm_bcd_set2.add(frm.bcd)
    
    shared_frm_id_set = frm_id_set1.intersection(frm_id_set2)
    shared_bcd_set = frm_bcd_set1.intersection(frm_bcd_set2)

    same_fragment_list = list()
    for frm in fragment_list1:
        if frm.frag_id in shared_frm_id_set:
            same_fragment_list.append(frm) 

    shared_fragment_dict1 = dict()
    shared_fragment_dict2 = dict()
    
    for frm in fragment_list1:
        if frm.bcd in shared_bcd_set and frm.frag_id not in shared_frm_id_set:
            if frm.bcd in shared_fragment_dict1:
                if frm.num_reads > shared_fragment_dict1[frm.bcd]:
                    shared_fragment_dict1[frm.bcd] = frm
            else:
                shared_fragment_dict1[frm.bcd] = frm

    for frm in fragment_list2:
        if frm.bcd in shared_bcd_set and frm.frag_id not in shared_frm_id_set:
            if frm.bcd in shared_fragment_dict2:
                if frm.num_reads > shared_fragment_dict2[frm.bcd]:
                    shared_fragment_dict2[frm.bcd] = frm
            else:
                shared_fragment_dict2[frm.bcd] = frm
    
    shared_fragment_list1 = list()
    shared_fragment_list2 = list()

    for bcd in shared_fragment_dict1:
        shared_fragment_list1.append(shared_fragment_dict1[bcd])
        shared_fragment_list2.append(shared_fragment_dict2[bcd])

    return shared_fragment_list1, shared_fragment_list2, same_fragment_list


if __name__ == '__main__':
    main()

