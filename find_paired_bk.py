#!/usr/bin/env python

import math
import numpy as np
from scipy.spatial import * # KDTree
from scipy.sparse import csr_matrix # csr_matrix
from scipy.sparse.csgraph import connected_components # connected_components
from my_utils import *
from fragment import *
from bedpe import *
from bed import *
import bisect
import gc


class TBed:
    def __init__(self, attr_list = None):
        if attr_list == None or len(attr_list) < 3: 
            self.tid = None
            self.start = None
            self.end = None
            return

        self.tid, self.start, self.end = attr_list[0:3]
        self.tid = int(self.tid)
        self.start = int(self.start)
        self.end = int(self.end)

    def key(self):
        return self.tid * FIX_LENGTH + self.start

    def extend_interval(self, interval_length):
        mean = (self.start + self.end)/2
        self.start = mean - interval_length/2
        self.end = mean + interval_length
        if self.start < 0: self.start = 0

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    find_paired_bk(args, dbo_args, endpoint_args)

    return 


def read_and_extend_tbed_file(in_tbed_file, max_gap_distance):

    in_tbed_fp = open(in_tbed_file, 'r')
    lines = list(in_tbed_fp)
    in_tbed_fp.close()

    tbed_list = list()
    for line in lines:
        line = line.strip().split(tab)
        tbed = TBed(line)
        tbed.extend_interval(max_gap_distance)
        tbed_list.append(tbed)

    return tbed_list

def merge_bk_cand_file(args, dbo_args, endpoint_args, max_gap_distance):

    if args.only_method1 == True:
        endpoint_tbed_list = list()
    else:
        endpoint_tbed_list = read_and_extend_tbed_file(endpoint_args.bk_file, max_gap_distance)

    if args.only_method2 == True:
        dbo_tbed_list = list()
    else:
        dbo_tbed_list = read_and_extend_tbed_file(dbo_args.bk_file, max_gap_distance)

    total_tbed_list = endpoint_tbed_list + dbo_tbed_list
    if len(total_tbed_list) == 0:
        return list()
    
    myprint ('number of raw tbed elements: %d' % len(total_tbed_list))
    total_tbed_list.sort(key=lambda x: x.key())
    merge_tbed_list = list() 
    merge_tbed_list.append(total_tbed_list[0])
    for i in range(1, len(total_tbed_list)):
        if total_tbed_list[i].tid == merge_tbed_list[-1].tid and total_tbed_list[i].start - merge_tbed_list[-1].end < 1:
            merge_tbed_list[-1].end = total_tbed_list[i].end
        else:
            merge_tbed_list.append(total_tbed_list[i])

    myprint ('number of merged tbed elements: %d' % len(merge_tbed_list))
    return merge_tbed_list
    
def create_nodes_for_frm_list(same_bcd_frm_list):

    node_list33 = list()
    node_list55 = list()
    node_list53 = list()
    node_list35 = list()
    same_bcd_frm_list.sort(key = lambda frm: frm.key_start())

    for i in range(0, len(same_bcd_frm_list)):
        for j in range(i+1, len(same_bcd_frm_list)):
            frm1 = same_bcd_frm_list[i]
            frm2 = same_bcd_frm_list[j]

            node_list33.append((frm1.key_end(),   frm2.key_end()))

            node_list55.append((frm1.key_start(), frm2.key_start()))

            node_list53.append((frm1.key_start(), frm2.key_end()))

            node_list35.append((frm1.key_end(),   frm2.key_start()))

    return node_list33, node_list55, node_list53, node_list35


def ouput_node_list2file(node_list, out_fp):
    for node in node_list:
        out_fp.write('%d\t%d\n' %  (node[0], node[1]))
    return

def build_graph_from_frm_list(args, all_potential_frm_list, min_support_fragments, max_gap_distance):

    myprint('building graph')

    all_potential_frm_list.sort(key = lambda frm: frm.bcd)
    if check_file_exists(args.node33_file) == False or check_file_exists(args.node55_file) == False or check_file_exists(args.node53_file) == False or check_file_exists(args.node35_file)== False:
        myprint('creating nodes')
        node33_fp = open(args.node33_file, 'w')
        node55_fp = open(args.node55_file, 'w')
        node53_fp = open(args.node53_file, 'w')
        node35_fp = open(args.node35_file, 'w')
        same_bcd_frm_list = list() 
        for frm in all_potential_frm_list:
            if len(same_bcd_frm_list) == 0:
                same_bcd_frm_list.append(frm)
                continue
            if frm.bcd == same_bcd_frm_list[0].bcd:
                same_bcd_frm_list.append(frm)
            else:
                node_list33, node_list55, node_list53, node_list35 = create_nodes_for_frm_list(same_bcd_frm_list) 
                ouput_node_list2file(node_list33, node33_fp)
                ouput_node_list2file(node_list55, node55_fp)
                ouput_node_list2file(node_list53, node53_fp)
                ouput_node_list2file(node_list35, node35_fp)
                same_bcd_frm_list = list()
                same_bcd_frm_list.append(frm)

        node33_fp.close()
        node55_fp.close()
        node53_fp.close()
        node35_fp.close()

        del node_list33, node_list55, node_list53, node_list35
        gc.collect()
    else:
        myprint ('node file existed. skipped creating nodes')

    out_file = args.bk_cand_pair_file

    out_fp = open(out_file, 'w')

    build_graph_from_node_list(args, args.node33_file, args.node_cluster33_file, min_support_fragments, max_gap_distance, out_fp, '3p_end', '3p_end')
    gc.collect()

    build_graph_from_node_list(args, args.node55_file, args.node_cluster55_file, min_support_fragments, max_gap_distance, out_fp, '5p_end', '5p_end')
    gc.collect()

    build_graph_from_node_list(args, args.node53_file, args.node_cluster53_file, min_support_fragments, max_gap_distance, out_fp, '5p_end', '3p_end')
    gc.collect()

    build_graph_from_node_list(args, args.node35_file, args.node_cluster35_file, min_support_fragments, max_gap_distance, out_fp, '3p_end', '5p_end')
    gc.collect()

    return


def read_node_list_file(node_list_file):
    
    node_list = list()
    node_list_fp = open(node_list_file, 'r')
    while 1: 
        line = node_list_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        node = (int(line[0]), int(line[1]))
        node_list.append(node)

    node_list_fp.close()

    return node_list

def read_node_cluster_file(node_cluster_file):
    node_cluster_list = list()
    node_cluster_fp = open(node_cluster_file, 'r')
    max_node_cluster_id = 0
    while 1:
        line = node_cluster_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        node_cluster_id = int(line[0])
        if max_node_cluster_id < node_cluster_id: 
            max_node_cluster_id = node_cluster_id

    n_cluster = max_node_cluster_id + 1
    node_cluster_fp.seek(0, 0)
    node_cluster_list = [0] * n_cluster
    for i in range(0, len(node_cluster_list)):
        node_cluster_list[i] = list()

    while 1:
        line = node_cluster_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        node_cluster_id = int(line[0])
        node = (int(line[1]), int(line[2]))
        if node_cluster_id > len(node_cluster_list)-1: print node_cluster_id, len(node_cluster_list)
        node_cluster_list[node_cluster_id].append(node)
    
    node_cluster_fp.close()
    return node_cluster_list

def build_graph_from_node_list(args, node_list_file, output_node_cluster_file, min_support_fragments, max_gap_distance, out_fp, endtype1, endtype2):
    
    if check_file_exists(output_node_cluster_file) == False:

        myprint('reading node file:%s' % node_list_file)
        node_list = read_node_list_file(node_list_file)

        edge_list = list()
        distance_buffer = max_gap_distance * 1.415
        myprint ('building KD tree, distance buffer = %d' % distance_buffer)
        tree = cKDTree(node_list, leafsize = 10000)

        myprint ('searching nearby nodes')
        for i in range(0, len(node_list)):
            if i % 1000000 == 1: myprint ('finished search for %d nodes' % i)
            node1 = node_list[i]
            index_list = tree.query_ball_point(node1, distance_buffer)
            if len(index_list) < min_support_fragments: continue
            nearby_node_index_list = list()
            for j in index_list:
                node2 = node_list[j]
                if abs(node1[0] - node2[0]) < max_gap_distance and abs(node1[1] - node2[1]) < max_gap_distance:
                    nearby_node_index_list.append(j)
            if len(nearby_node_index_list) < min_support_fragments: continue
            for j in nearby_node_index_list: 
                edge = (i, j)    
                edge_list.append(edge)
            
        row = list()
        col = list()
        data = list()
         
        for edge in edge_list:
            row.append (edge[0])
            col.append (edge[1])
            data.append (1)

        n_node = len(node_list)

        myprint ('get connected components')
        n_components, label_list, component_node_index_db = get_connected_components(n_node, row, col, data, True, 'strong')

        node_cluster_list = [0] * n_components
        for i in range(0, n_components):
            node_cluster_list[i] = list()
            for index in component_node_index_db[i]:
                node_cluster_list[i].append(node_list[index])
                
        ## output node cluster file ## 
        output_node_cluster_fp = open(output_node_cluster_file, 'w')
        for i in range(0, len(node_cluster_list)):
            node_cluster = node_cluster_list[i]
            for j in range(0, len(node_cluster)):
                node = node_cluster[j]
                output_node_cluster_fp.write('%d\t%d\t%d\n' % (i, node[0], node[1])) 
        output_node_cluster_fp.close()

        del node_list, edge_list, row, col, data, component_node_index_db

    else:
        myprint ('reading node cluster file: %s' % output_node_cluster_file)
        node_cluster_list = read_node_cluster_file(output_node_cluster_file) 

    ## output paired_bk_cand_file ##
    paired_bk_cand_list = list()
    for node_cluster in node_cluster_list:
        if len(node_cluster) < min_support_fragments: continue
        paired_bk_cand = convert_node_cluster_to_paired_bk_cand(node_cluster, max_gap_distance, args.tid2chrname, endtype1, endtype2)
        paired_bk_cand_list.append(paired_bk_cand)

    for paired_bk_cand in paired_bk_cand_list:
        out_fp.write(paired_bk_cand.output() + endl)

    return 

def get_connected_components(n_node, row, col, data, is_directed = False, connection_type = 'week'):

    node_csr_matrix = csr_matrix((data, (row, col)), shape=[n_node, n_node])
    n_components, label_list = connected_components(node_csr_matrix, directed = is_directed, connection = connection_type)
    component_node_index_db = [0] * n_components
    for i in range(0, len(component_node_index_db)): 
        component_node_index_db[i] = list()
    # component_node_index_db[component_id] = index of node
    for i in range(0, len(label_list)):
        component_node_index_db[label_list[i]].append(i)

    return n_components, label_list, component_node_index_db 


def get_tid_pos_from_key(key):

    tid = int(key / FIX_LENGTH)
    pos = key % FIX_LENGTH
    return tid, pos

def convert_node_cluster_to_paired_bk_cand(node_cluster, max_gap_distance, tid2chrname, endtype1, endtype2):
    
    xlist = list()
    ylist = list()
    for node in node_cluster: 
        xlist.append(node[0]) 
        ylist.append(node[1]) 
    
    bin_size = 100
    win_size = 50 
    xkey_start, xkey_end = get_max_density_region(xlist, win_size, bin_size) 
    ykey_start, ykey_end = get_max_density_region(ylist, win_size, bin_size)

    xtid, xstart = get_tid_pos_from_key(xkey_start)
    xtid, xend   = get_tid_pos_from_key(xkey_end)

    ytid, ystart = get_tid_pos_from_key(ykey_start) 
    ytid, yend   = get_tid_pos_from_key(ykey_end) 

    xchr = tid2chrname[xtid]
    ychr = tid2chrname[ytid]

    svtype = 'UNK'
     
    if xtid == ytid:
        svlength = str(abs(ystart-xstart))
    else:
        svlength = 'NA' 

    core_node_cluster = list()
    for node in node_cluster:
        if node[0] >= xkey_start and node[0] <= xkey_end and node[1] >= ykey_start and node[1] <= ykey_end:
            core_node_cluster.append(node)

    attr_list = [xchr, xstart, xend, ychr, ystart, yend, svtype, svlength, endtype1, endtype2, len(core_node_cluster)] 
    paired_bk_cand = PairedBkCand(attr_list) 
    return paired_bk_cand

def get_max_density_region(alist, win_size, bin_size):

    amin = min(alist)
    blist = list()
    for a in alist:
        b = int((a-amin)/bin_size)
        blist.append(b)
    bmax = max(blist)
    pmf = [0] * (bmax+1) 
    for b in blist:
        pmf[b] += 1

    csum = [0] * (bmax+1)
    csum[0] = pmf[0]
    for i in range(1, len(pmf)):
        csum[i] = csum[i-1] + pmf[i]
   
    if len(csum) < win_size:  
        return amin, max(alist)

    window_sum = [0] * len(csum)
    window_sum[win_size-1] = csum[win_size-1]
    for right_window_idx in range(win_size, len(csum)):
        window_sum[right_window_idx] = csum[right_window_idx] - csum[right_window_idx-win_size]

    max_win_sum_idx = 0
    max_win_sum = window_sum[max_win_sum_idx]
    for i in range(1, len(window_sum)):
        if max_win_sum < window_sum[i]: 
            max_win_sum_idx = i
            max_win_sum = window_sum[max_win_sum_idx]

    aright = (max_win_sum_idx+1) * bin_size + amin 
    aleft = aright - win_size * bin_size
    return aleft - bin_size, aright + bin_size

def find_paired_bk(args, dbo_args, endpoint_args):

    max_gap_distance = min(10000, args.gap_distance_cutoff)

    bcd22_file = endpoint_args.bcd22_file
    if args.all_to_all == False:
        merge_bk_tbed_list = merge_bk_cand_file(args, dbo_args, endpoint_args, max_gap_distance)
    else:
        merge_bk_tbed_list = list()


    myprint('reading bcd22 file:%s' % bcd22_file)
    all_potential_frm_list = get_supp_fragment(args, merge_bk_tbed_list, bcd22_file, endpoint_args.min_frag_length)

    myprint('searching paired breakpoints')
    build_graph_from_frm_list(args, all_potential_frm_list, args.min_support_fragments, max_gap_distance)
 
    return

def get_supp_fragment(args, merge_bk_tbed_list, bcd22_file, min_frag_length):

    bcd22_frm_list = read_bcd22_file_core(bcd22_file, min_frag_length)

    myprint('total number of fragments: %d' % (len(bcd22_frm_list)))
    if args.all_to_all: return bcd22_frm_list

    start_sorted_frm_list = sorted(bcd22_frm_list, key = lambda frm: frm.key_start())
    start_key_list = list()
    for frm in start_sorted_frm_list:
        start_key_list.append(frm.key_start())

    end_sorted_frm_list = sorted(bcd22_frm_list, key = lambda frm: frm.key_end())
    end_key_list = list()
    for frm in end_sorted_frm_list:
        end_key_list.append(frm.key_end())

    all_potential_frm_list = list()
    for tbed in merge_bk_tbed_list:
        ## get all fragments of which the endpoints in the interval
        searchstart_index = bisect.bisect_left(start_key_list, tbed.tid * FIX_LENGTH + tbed.start) 
        searchend_index = bisect.bisect_right(start_key_list, tbed.tid * FIX_LENGTH + tbed.end)
        for i in range(searchstart_index, searchend_index):
            all_potential_frm_list.append(start_sorted_frm_list[i])

        searchstart_index = bisect.bisect_left(end_key_list, tbed.tid * FIX_LENGTH + tbed.start)
        searchend_index = bisect.bisect_right(end_key_list, tbed.tid * FIX_LENGTH + tbed.end)
        for i in range(searchstart_index, searchend_index):
            all_potential_frm_list.append(end_sorted_frm_list[i])


    all_potential_frm_bcd_set = set()
    all_potential_frm_id_set = set()
    for frm in all_potential_frm_list: 
        all_potential_frm_bcd_set.add(frm.bcd)
        all_potential_frm_id_set.add(frm.frag_id)

    for frm in bcd22_frm_list:
        if frm.bcd in all_potential_frm_bcd_set and frm.frag_id not in all_potential_frm_id_set:
            all_potential_frm_list.append(frm)
            all_potential_frm_id_set.add(frm.frag_id)

    if len(all_potential_frm_list) == 0: return list()

    all_potential_frm_list.sort(key = lambda frm: frm.frag_id)

    rmdup_all_potential_frm_list = list()
    rmdup_all_potential_frm_list.append(all_potential_frm_list[0])

    for i in range(1, len(all_potential_frm_list)):
        if all_potential_frm_list[i].frag_id != rmdup_all_potential_frm_list[-1].frag_id:
            rmdup_all_potential_frm_list.append(all_potential_frm_list[i])

    myprint('potential candidate fragments: %d' % (len(bcd22_frm_list), len(rmdup_all_potential_frm_list)))
    del all_potential_frm_bcd_set, all_potential_frm_id_set, bcd22_frm_list, start_sorted_frm_list, end_sorted_frm_list, start_key_list, end_key_list
    return rmdup_all_potential_frm_list
    
if __name__ == '__main__':
    main()

