#!/usr/bin/env python

import os
import sys
from bedpe import *
from my_utils import *
from scipy import spatial
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import numpy as np
import random
from scipy.stats import binom

tab  = '\t'
endl = '\n'


class AlignedRead:

    def __init__(self, sam_line):

        line = sam_line.strip().split(tab)
        self.readname = line[0]
        self.flag = int(line[1])
        self.chrm = line[2]
        self.left_pos = int(line[3])
        self.right_pos = None
        self.mapq = int(line[4])
        self.cigar_string = line[5]
        self.right_pos = None
        self.cigar_operation_list = None
        self.cigar_operation_length_list = None
        self.n_cigar = None
        self.end3p_pos = None
        self.cigar_analysis()
        self.get_right_ref_pos() 
        self.get_end3p_pos()

        self.bcd = ''    
        self.hp = 0
        for i in range(11, len(line)):
            a = line[i].split(':')
            if a[0] == 'BX': self.bcd = a[2]
            if a[0] == 'HP': self.hp = int(a[2])

    def output_core(self):
        outstring = '%s\t%d\t%s\t%d\t%s\t%d' % (self.readname, self.flag, self.chrm, self.left_pos, self.cigar_string, self.n_cigar)
        return outstring
        
    def map_orientation(self):

        if self.flag & 0x10:
            return -1
        else:
            return 1 

    def cigar_analysis(self):
        p1 = p2 = 0 
        number_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        n_cigar = 0 
        self.cigar_operation_list = list()
        self.cigar_operation_length_list = list()
        for i in range(0, len(self.cigar_string)):
            c = self.cigar_string[i]
            if c in number_list: continue
            p2 = i 
            cigar_operation = self.cigar_string[p2]
            self.cigar_operation_list.append(cigar_operation)
            cigar_operation_length = int(self.cigar_string[p1:p2])
            self.cigar_operation_length_list.append(cigar_operation_length)
            p1 = p2 + 1 
            p2 = p2 + 1 
    
        self.n_cigar = len(self.cigar_operation_list)


    def get_right_ref_pos(self):
        if self.cigar_operation_list == None:
            self.cigar_analysis()

        ref_align_length = 0
        affect_ref_align_cigar = ['M', 'D', 'N', '=', 'X']
        for i in range(0, self.n_cigar):
            if self.cigar_operation_list[i] in affect_ref_align_cigar:
                ref_align_length += self.cigar_operation_length_list[i]
        
        self.right_pos = self.left_pos + ref_align_length - 1 
        
    def has_left_clip(self):

        if self.n_cigar == 0 or self.n_cigar == None: 
            return False

        if self.cigar_operation_list[0] == 'S' or self.cigar_operation_list[0] == 'H':
            return True
        else:
            return False

    def has_right_clip(self): 
        if self.n_cigar == 0 or self.n_cigar == None: 
            return False

        if self.cigar_operation_list[-1] == 'S' or self.cigar_operation_list[-1] == 'H':
            return True
        else:
            return False

    def has_clip(self):
        return (self.has_left_clip() or self.has_right_clip())

    def get_end3p_pos(self):
        if self.flag & 0x10:
            self.end3p_pos = self.left_pos 
        else:
            self.end3p_pos = self.right_pos

    def major_clip_side(self):

        if self.has_left_clip() == True and self.has_right_clip() == False: 
            return 'left'
        elif self.has_left_clip() == False and self.has_right_clip() == True: 
            return 'right'
        else:
            left_clip_length = self.cigar_operation_length_list[0]
            right_clip_length = self.cigar_operation_length_list[-1]
            if left_clip_length > right_clip_length:
                return 'left'
            elif left_clip_length < right_clip_length:
                return 'right'
            else:
                randnumber = random.randint(1,2)
                if randnumber == 1:
                    return 'left'
                else:
                    return 'right'

class PairedEndSupport:

    def __init__(self, aligned_read1 = None, aligned_read2 = None):

        self.aligned_read1 = None
        self.aligned_read2 = None
        self.endtype1 = None
        self.endtype2 = None

        if aligned_read1 != None:
            self.aligned_read1 = aligned_read1
            self.aligned_read2 = aligned_read2

            if self.aligned_read1.map_orientation() == 1:
                self.endtype1 = '3p_end'
            else: 
                self.endtype1 = '5p_end'

            if self.aligned_read2.map_orientation() == 1:
                self.endtype2 = '3p_end'
            else: 
                self.endtype2 = '5p_end'

class SplitReadSupport:

    def __init__(self, aligned_read1 = None, aligned_read2 = None): 

        self.aligned_read1 = None
        self.aligned_read2 = None
        self.endtype1 = None
        self.endtype2 = None
        self.major_clip_pos1 = None
        self.major_clip_pos2 = None

        if aligned_read1 != None:
            self.aligned_read1 = aligned_read1
            self.aligned_read2 = aligned_read2

            if self.aligned_read1.major_clip_side() == 'left': 
                self.endtype1 = '5p_end'
            else:
                self.endtype1 = '3p_end'

            if self.aligned_read2.major_clip_side() == 'left': 
                self.endtype2 = '5p_end'
            else:
                self.endtype2 = '3p_end'

            if self.aligned_read1.major_clip_side() == 'left':
                self.major_clip_pos1 = self.aligned_read1.left_pos
            else:
                self.major_clip_pos1 = self.aligned_read1.right_pos

            if self.aligned_read2.major_clip_side() == 'left':
                self.major_clip_pos2 = self.aligned_read2.left_pos
            else:
                self.major_clip_pos2 = self.aligned_read2.right_pos

            
class PESupportGroup:

    def __init__(self, pe_support_list = None, endtype1 = None, endtype2 = None, median_readpair_gap_distance = None):
        self.pe_support_list = list()
        self.n_pe_support    = 0
        self.endtype1 = None
        self.endtype2 = None
        self.bk1_pos = -1
        self.bk2_pos = -1
        self.median_end3p_pos1 = -1
        self.median_end3p_pos2 = -1
        self.map_orientation1  = 0
        self.map_orientation2  = 0

        if pe_support_list != None and median_readpair_gap_distance != None:
            self.pe_support_list = pe_support_list # list of objects of class PairedEndSupport 
            self.n_pe_support = len(self.pe_support_list)
            self.endtype1 = endtype1
            self.endtype2 = endtype1
            self.map_orientation1 = self.pe_support_list[0].aligned_read1.map_orientation() 
            self.map_orientation2 = self.pe_support_list[0].aligned_read2.map_orientation() 
            self.predict_bk_pos(median_readpair_gap_distance) 

    def predict_bk_pos(self, median_readpair_gap_distance):
        if len(self.pe_support_list) < 1: return
        end3p_pos_list1 = list() 
        end3p_pos_list2 = list()
        for pe_support in self.pe_support_list: 
            end3p_pos_list1.append(pe_support.aligned_read1.end3p_pos)
            end3p_pos_list2.append(pe_support.aligned_read2.end3p_pos)

        self.median_end3p_pos1 = np.median(end3p_pos_list1)
        self.median_end3p_pos2 = np.median(end3p_pos_list2)

        self.bk1_pos = self.median_end3p_pos1 + self.map_orientation1 * median_readpair_gap_distance / 2
        self.bk2_pos = self.median_end3p_pos2 + self.map_orientation2 * median_readpair_gap_distance / 2 
        return

class SRSupportGroup:

    def __init__(self, sr_support_list = None, endtype1 = None, endtype2 = None):

        self.sr_support_list = list()
        self.n_sr_support = 0
        self.endtype1 = None
        self.endtype2 = None
        self.bk1_pos = -1
        self.bk2_pos = -1

        if sr_support_list != None:
            self.sr_support_list = sr_support_list
            self.n_sr_support = len(self.sr_support_list)
            self.endtype1 = endtype1 
            self.endtype2 = endtype2 
            self.bk1_pos = None
            self.bk2_pos = None
            self.predict_bk_pos()

    def predict_bk_pos(self):

        major_clip_pos_list1 = list()
        major_clip_pos_list2 = list()
        for sr_support in self.sr_support_list:
            major_clip_pos_list1.append(sr_support.major_clip_pos1)
            major_clip_pos_list2.append(sr_support.major_clip_pos2)

        self.bk1_pos = np.median(major_clip_pos_list1)
        self.bk2_pos = np.median(major_clip_pos_list2)

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    import global_distribution
    myprint('estimation global parameters')
    global_distribution.estimate_global_distribution(args, dbo_args, endpoint_args)

    refine_breakpoints(args, dbo_args, endpoint_args)

    return

def get_isize_distribution(args):

    out_dir = os.path.join(args.out_dir, 'refinebreakpoints')
    
    os.system('mkdir -p %s' % out_dir)
    tmp_file = os.path.join(out_dir, 'isize.txt') 
    tmp_sh_file = os.path.join(out_dir, 'tmp.sh') 
    tmp_sh_fp = open(tmp_sh_file, 'w')
    tmp_sh_fp.write('#!/bin/bash\n')
    cmd = '%s view %s | cut -f 9 | head -1000000 > %s\n' % (args.samtools, args.bam, tmp_file) 
    tmp_sh_fp.write(cmd)
    tmp_sh_fp.close()

    os.system('sh %s ' % tmp_sh_file)
    tmp_fp = open(tmp_file, 'r') 
    isize_list = list()
    while 1:
        line = tmp_fp.readline()
        if not line: break
        isize = int(line.strip())
        if isize > 0: isize_list.append(isize)
    tmp_fp.close()

    percentiles = [0, 25, 50, 75, 100]
    q = np.percentile(isize_list, percentiles)
    args.median_isize = int(q[2]) 
    args.max_isize_cutoff = int(q[3] + 1.5 * (q[3] - q[1]))
    myprint('median insert size is: %d, isize cut-off is: %d' % (args.median_isize, args.max_isize_cutoff))
    os.system('rm -rf %s' % tmp_file)

    return  

def get_read_length(args):

    out_dir = os.path.join(args.out_dir, 'refinebreakpoints')

    tmp_sh_file = os.path.join(out_dir, 'tmp.sh') 
    tmp_sh_fp = open(tmp_sh_file, 'w')
    tmp_sh_fp.write('#!/bin/bash\n')
    tmp_file = os.path.join(out_dir, 'tmp.sam') 
    cmd = '%s view %s | head -10000 > %s\n' % (args.samtools, args.bam, tmp_file)
    tmp_sh_fp.write(cmd)
    tmp_sh_fp.close()
    os.system('sh %s ' % tmp_sh_file)

    tmp_fp = open(tmp_file, 'r') 

    rd_length_list1 = list()
    rd_length_list2 = list()

    while 1:
        line = tmp_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        flag = int(line[1])
        read_length = len(line[9]) 
        if flag & 0x40: rd_length_list1.append(read_length)
        if flag & 0x80: rd_length_list2.append(read_length)
    tmp_fp.close()

    args.rd1_length = int(np.median(rd_length_list1))
    args.rd2_length = int(np.median(rd_length_list2))

    myprint('read1 length is %d, read2 length is %d' % (args.rd1_length, args.rd2_length))
    os.system('rm -rf %s' % tmp_file)

    return

def refine_breakpoints(args, dbo_args, endpoint_args):
    
    myprint('refine breakpoints')
    in_bam = args.bam
    in_bam_prefix = os.path.splitext(in_bam)[0]

    bai_file1 = in_bam_prefix + '.bai'  
    bai_file2 = in_bam_prefix + '.bam.bai' 

    if os.path.isfile(bai_file1) == False and os.path.isfile(bai_file2) == False:
        myprint('bam index file not found, indexing using samtools')
        cmd = args.samtools + ' index ' + in_bam
        os.system(cmd)

    get_isize_distribution(args)  
    get_read_length(args)

    svcall_bedpe_list = read_object_file(args.quantified_bk_pair_file, QuantifiedBKCand)

    myprint ('number of SV calls: %d' % len(svcall_bedpe_list))

    refined_svcall_list = list()
    out_file = args.refinedbedpe_file
    out_fp = open(out_file, 'w')
    for svcall_bedpe in svcall_bedpe_list:
        if svcall_bedpe.score > 0 and svcall_bedpe.logp_nosv_two_mol > svcall_bedpe.logp_nosv_one_mol:
            refined_svcall = refine1sv(args, dbo_args, endpoint_args, svcall_bedpe) 
            refined_svcall_list.append(refined_svcall)
            out_fp.write(refined_svcall.output() + endl)

    out_fp.close()

    return

def get_refined_sv_call(args, svcall_bedpe, pe_support_group, sr_support_group): 

    n_pe_support = pe_support_group.n_pe_support
    n_sr_support = sr_support_group.n_sr_support
    n_support_reads = n_pe_support + n_sr_support 
    n_support_fragment = svcall_bedpe.num_fragment_support

    refined_bk1 = -1 
    refined_bk2 = -1 

    if sr_support_group.bk1_pos != None and sr_support_group.bk2_pos != None:
        refined_bk1 = sr_support_group.bk1_pos
        refined_bk2 = sr_support_group.bk2_pos
    elif pe_support_group.bk1_pos != None and pe_support_group.bk2_pos != None:    
        refined_bk1 = pe_support_group.bk1_pos
        refined_bk2 = pe_support_group.bk2_pos

    median_isize = args.median_isize
    median_inter_readpair_gap_distance = args.gap_distance750
    p0 = float(median_isize/median_inter_readpair_gap_distance)

    if n_support_reads > n_support_fragment:
        logp = 0 
    else:
        logp = binom.logcdf(n_support_reads, n_support_fragment, p0) / np.log(10)  

    refine_score = logp  

    attr_list = svcall_bedpe.attr_list() + [refine_score, n_pe_support, n_sr_support, n_support_reads, refined_bk1, refined_bk2] 

    refinded_sv_call = RefinedQuantBKCand(attr_list)

    return refinded_sv_call
    
def refine1sv(args, dbo_args, endpoint_args, svcall_bedpe):
   
    max_split_pos_variation = 50
    in_bam = args.bam
    search_range = args.gap_distance_cutoff
    
    endtype1 = svcall_bedpe.endtype1
    endtype2 = svcall_bedpe.endtype2
    support_barcode_set = set(svcall_bedpe.support_barcodes.strip(',').split(','))

    out_dir = os.path.join(args.out_dir, 'refinebreakpoints')

    os.system('mkdir -p %s' % out_dir)

    myprint('current sv call is: %s' % (svcall_bedpe.output_svcall()))
    search_start1 = svcall_bedpe.start1 - search_range
    search_end1   = svcall_bedpe.end1 + search_range 

    search_start2 = svcall_bedpe.start2 - search_range
    search_end2   = svcall_bedpe.end2 + search_range 

    bam_name  = os.path.split(in_bam)[1]

    region1_sam_file = os.path.join(out_dir, bam_name + '.%s_%d_%d.sam' % (svcall_bedpe.chrm1, search_start1, search_end1))
    cmd1 = '%s view %s %s:%d-%d > %s' % (args.samtools, in_bam, svcall_bedpe.chrm1, search_start1, search_end1, region1_sam_file) 
    os.system(cmd1)

    region2_sam_file = os.path.join(out_dir, bam_name + '.%s_%d_%d.sam' % (svcall_bedpe.chrm2, search_start2, search_end2))
    cmd2 = '%s view %s %s:%d-%d > %s' % (args.samtools, in_bam, svcall_bedpe.chrm2, search_start2, search_end2, region2_sam_file) 
    os.system(cmd2)

    region1_read_list = list() 
    region1_sam_fp = open(region1_sam_file, 'r')
    while 1:
        line = region1_sam_fp.readline()
        if not line: break
        aligned_read = AlignedRead(line)
        if aligned_read.bcd not in support_barcode_set: continue
        region1_read_list.append(aligned_read)

    region1_sam_fp.close()

        
    region2_read_list = list() 
    region2_sam_fp = open(region2_sam_file, 'r')
    while 1:
        line = region2_sam_fp.readline()
        if not line: break
        aligned_read = AlignedRead(line)
        if aligned_read.bcd not in support_barcode_set: continue
        region2_read_list.append(aligned_read)

    region2_sam_fp.close()

    myprint ('number of reads in region1 and region2: %d, %d' % (len(region1_read_list), len(region2_read_list)))

    pe_support_list = list()
    sr_support_list = list() 

    for aligned_read1 in region1_read_list:
        for aligned_read2 in region2_read_list:
            if aligned_read1.readname != aligned_read2.readname: continue
            if aligned_read1.flag & 0x40 != aligned_read2.flag & 0x40:
                pe_support = PairedEndSupport(aligned_read1, aligned_read2)
                if pe_support.endtype1 == endtype1 and pe_support.endtype2 == endtype2: pe_support_list.append(pe_support)
            else:
                if not aligned_read1.has_clip() or not aligned_read2.has_clip(): continue
                sr_support = SplitReadSupport(aligned_read1, aligned_read2)
                if sr_support.endtype1 == endtype1 and sr_support.endtype2 == endtype2: sr_support_list.append(sr_support)
    
    pe_support_group = analyze_pe_support_list(args, pe_support_list, svcall_bedpe) 
    sr_support_group = analyze_sr_support_list(sr_support_list, svcall_bedpe, max_split_pos_variation) 

    refined_sv_call = get_refined_sv_call(args, svcall_bedpe, pe_support_group, sr_support_group)

    os.system('rm ' + region1_sam_file)
    os.system('rm ' + region2_sam_file)

    return refined_sv_call

def analyze_sr_support_list(sr_support_list, svcall_bedpe, max_split_pos_variation):

    if len(sr_support_list) == 0: return SRSupportGroup()

    myprint ('analysis of split-read support' )
    edge_list  = list()
    row = list()
    col = list()
    data = list()
    
    for i in range(0, len(sr_support_list)):
        for j in range(i+1, len(sr_support_list)):
            sr_support1 = sr_support_list[i]
            sr_support2 = sr_support_list[j]
            dist1 = abs(sr_support1.major_clip_pos1 - sr_support2.major_clip_pos1)
            dist2 = abs(sr_support1.major_clip_pos2 - sr_support2.major_clip_pos2) 
            if dist1 < max_split_pos_variation and dist2 < max_split_pos_variation:
                graph_weight = dist1 + dist2
                edge = (i, j, graph_weight)
                edge_list.append(edge)

    for edge in edge_list:
        row.append (edge[0])
        col.append (edge[1])
        data.append (edge[2])
            
    n_node = len(sr_support_list)
    sr_support_csr_matrix = csr_matrix((data, (row, col)), shape=[n_node, n_node])
    n_components, label_list = connected_components(sr_support_csr_matrix, directed = False)
    component_element_db = [0] * n_components
    for i in range(0, len(component_element_db)): 
        component_element_db[i] = list()
    # component_element_db[component_id] = index of pe_support
    for i in range(0, len(label_list)):
        component_element_db[label_list[i]].append(i)
    
    sr_support_group_list = list()
    for component_id in range(0, len(component_element_db)):
        sr_support_sublist = list()
        for index in component_element_db[component_id]:
            sr_support = sr_support_list[index] 
            sr_support_sublist.append(sr_support)

        sr_support_group = SRSupportGroup(sr_support_sublist, svcall_bedpe.endtype1, svcall_bedpe.endtype2)
        sr_support_group_list.append(sr_support_group)

    if len(sr_support_group_list) == 0: return SRSupportGroup()

    max_group_id = 0
    for i in range(1, len(sr_support_group_list)): 
        if sr_support_group_list[i].n_sr_support > sr_support_group_list[max_group_id]: max_group_id = i

    return sr_support_group_list[max_group_id]
    

def analyze_pe_support_list (args, pe_support_list, svcall_bedpe):

    if len(pe_support_list) == 0: return PESupportGroup()

    myprint ('analysis of PE support')
    max_isize_cutoff = args.max_isize_cutoff
    median_readpair_gap_distance = args.median_isize - (args.rd1_length + args.rd2_length)

    if median_readpair_gap_distance < 0: median_readpair_gap_distance = 0

    endtype1 = svcall_bedpe.endtype1
    endtype2 = svcall_bedpe.endtype2
    edge_list  = list()
    row = list()
    col = list()
    data = list()

    for i in range(0, len(pe_support_list)):
        for j in range(i+1, len(pe_support_list)):
            pe_support1 = pe_support_list[i]
            pe_support2 = pe_support_list[j]
            dist1 = abs(pe_support1.aligned_read1.end3p_pos - pe_support2.aligned_read1.end3p_pos)
            dist2 = abs(pe_support1.aligned_read2.end3p_pos - pe_support2.aligned_read2.end3p_pos)
            map_ori_match = False
            if pe_support1.aligned_read1.map_orientation() == pe_support2.aligned_read1.map_orientation() and pe_support1.aligned_read2.map_orientation() == pe_support2.aligned_read2.map_orientation(): map_ori_match = True
            if dist1 < max_isize_cutoff and dist2 < max_isize_cutoff and map_ori_match == True:
                graph_weight = dist1 + dist2 
                edge = (i, j, graph_weight) 
                edge_list.append(edge)

    for edge in edge_list:
        row.append (edge[0])
        col.append (edge[1])
        data.append (edge[2])

    n_node = len(pe_support_list)
    pe_support_csr_matrix = csr_matrix((data, (row, col)), shape=[n_node, n_node])
    n_components, label_list = connected_components(pe_support_csr_matrix, directed=False)
    component_element_db = [0] * n_components
    for i in range(0, len(component_element_db)): 
        component_element_db[i] = list()
    # component_element_db[component_id] = index of pe_support
    for i in range(0, len(label_list)):
        component_element_db[label_list[i]].append(i)
    
    pe_support_group_list = list()
    for component_id in range(0, len(component_element_db)):
        pe_support_sublist = list()
        for index in component_element_db[component_id]:
            pe_support = pe_support_list[index] 
            pe_support_sublist.append(pe_support)

        pe_support_group = PESupportGroup(pe_support_sublist, endtype1, endtype2, median_readpair_gap_distance)
        pe_support_group_list.append(pe_support_group)

    if len(pe_support_group_list) == 0: return PESupportGroup()

    max_group_id = 0
    for i in range(1, len(pe_support_group_list)): 
        if pe_support_group_list[i].n_pe_support > pe_support_group_list[max_group_id]: max_group_id = i

    return pe_support_group_list[max_group_id]
    

if __name__ == '__main__':
    main()
