#!/usr/bin/env python


import math
import numpy as np
from scipy.spatial import * # KDTree
from scipy.sparse import csr_matrix # csr_matrix
from scipy.sparse.csgraph import connected_components # connected_components
import bisect
import gc

try:
    from scripts.my_utils import *
except ImportError:
    from my_utils import *

try:
    from scripts.fragment import *
except ImportError:
    from fragment import *

try:
    from scripts.bedpe import *
except ImportError:
    from bedpe import *

try:
    from scripts.bed import *
except ImportError:
    from bed import *


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

class EndpointNode:
    def __init__(self, x, y, x_frag_id, y_frag_id):
        self.x = int(x)
        self.y = int(y)
        self.x_frag_id = int(x_frag_id)
        self.y_frag_id = int(y_frag_id)

    def output(self):
        outstring = '%d\t%d\t%d\t%d' % (self.x, self.y, self.x_frag_id, self.y_frag_id)
        return outstring

    def csvoutput(self):
        outstring = '%d,%d,%d,%d' % (self.x, self.y, self.x_frag_id, self.y_frag_id)
        return outstring

class FragmentPair:

    def __init__(self, frm_pair_id, frm1_id, frm2_id):
        self.frm_pair_id = int(frm_pair_id)
        self.frm1_id = int(frm1_id)
        self.frm2_id = int(frm2_id)

    def id_key(self):
        return '%d\t%d' % (self.frm1_id, self.frm2_id)

def main():

    args, dbo_args, endpoint_args = parse_user_arguments()

    find_paired_bk(args, dbo_args, endpoint_args)

    return 


def create_nodes_for_frm_list(same_bcd_frm_list, min_frag_length):

    node_list33 = list()
    node_list55 = list()
    node_list53 = list()
    node_list35 = list()
    same_bcd_frm_list.sort(key = lambda frm: frm.key_start())

    for i in range(0, len(same_bcd_frm_list)):
        for j in range(i+1, len(same_bcd_frm_list)):
            frm1 = same_bcd_frm_list[i]
            frm2 = same_bcd_frm_list[j]
            if frm1.length + frm2.length < min_frag_length: continue

            node33 = EndpointNode(frm1.key_end(), frm2.key_end(), frm1.frag_id, frm2.frag_id) 
            node55 = EndpointNode(frm1.key_start(), frm2.key_start(), frm1.frag_id, frm2.frag_id) 
            node53 = EndpointNode(frm1.key_start(), frm2.key_end(), frm1.frag_id, frm2.frag_id) 
            node35 = EndpointNode(frm1.key_end(), frm2.key_start(), frm1.frag_id, frm2.frag_id) 

            node_list33.append(node33)
            node_list55.append(node55)
            node_list53.append(node53)
            node_list35.append(node35)

    return node_list33, node_list55, node_list53, node_list35


def ouput_node_list2file(node_list, out_fp):
    for node in node_list:
        out_fp.write(node.output() + endl)
    return

def build_graph_from_fragments (args, dbo_args, endpoint_args):

    if args.run_from_begining == True or (check_file_exists(args.node33_file) == False or check_file_exists(args.node55_file) == False or check_file_exists(args.node53_file) == False or check_file_exists(args.node35_file)== False):

        myprint('building nodes from fragments')
        myprint('reading bcd22 file:%s' % endpoint_args.bcd22_file)
        if args.is_wgs == True:
            min_frag_length = endpoint_args.min_frag_length
        else:
            min_frag_length = 0
        all_potential_frm_list = read_bcd22_file_core(endpoint_args.bcd22_file, min_frag_length) # all fragments that are longer than min_frag_length
        myprint('total number of fragments: %d' % (len(all_potential_frm_list)))
        all_potential_frm_list.sort(key = lambda frm: frm.bcd)
        myprint('writing to node file')
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
                node_list33, node_list55, node_list53, node_list35 = create_nodes_for_frm_list(same_bcd_frm_list, endpoint_args.min_frag_length) 
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

        del node_list33, node_list55, node_list53, node_list35, same_bcd_frm_list
        del all_potential_frm_list
        gc.collect()
    else:
        myprint ('node file existed. skipped creating nodes')

    gc.collect()

    max_gap_distance = args.gap_distance_cutoff

    if args.run_from_begining == True or (check_file_exists(args.node33_candidate_file) == False or check_file_exists(args.node55_candidate_file) == False or check_file_exists(args.node53_candidate_file) == False or check_file_exists(args.node35_candidate_file)== False): 

        myprint ('removing sparse nodes, min_support_fragments is %d' % args.min_support_fragments )
        cmd = '%s %s %s %d %s %d' % (args.remove_sparse_nodes, args.node33_file, args.node33_candidate_file, max_gap_distance, args.faidx_file, args.min_support_fragments) 
        run_system_cmd(cmd)
        cmd = '%s %s %s %d %s %d' % (args.remove_sparse_nodes, args.node55_file, args.node55_candidate_file, max_gap_distance, args.faidx_file, args.min_support_fragments) 
        run_system_cmd(cmd)
        cmd = '%s %s %s %d %s %d' % (args.remove_sparse_nodes, args.node35_file, args.node35_candidate_file, max_gap_distance, args.faidx_file, args.min_support_fragments) 
        run_system_cmd(cmd)
        cmd = '%s %s %s %d %s %d' % (args.remove_sparse_nodes, args.node53_file, args.node53_candidate_file, max_gap_distance, args.faidx_file, args.min_support_fragments) 
        run_system_cmd(cmd)

        if os.path.exists(args.node33_candidate_file) == False: 
            myprint ('ERROR! %s was not generated' % args.node33_candidate_file)
            sys.exit()
        if os.path.exists(args.node55_candidate_file) == False: 
            myprint ('ERROR! %s was not generated' % args.node55_candidate_file)
            sys.exit()
        if os.path.exists(args.node35_candidate_file) == False: 
            myprint ('ERROR! %s was not generated' % args.node35_candidate_file)
            sys.exit()
        if os.path.exists(args.node53_candidate_file) == False: 
            myprint ('ERROR! %s was not generated' % args.node53_candidate_file)
            sys.exit()

    else:
        myprint ('node candidate file existed. skipped removing sparse nodes')

    if args.rm_temp_files:
        args.temp_file_list.append(args.node55_file)
        args.temp_file_list.append(args.node53_file)
        args.temp_file_list.append(args.node35_file)
        for temp_file in args.temp_file_list:
            if check_file_exists(temp_file): os.remove(temp_file)

        args.temp_file_list = list()
        

    myprint ('clustering nodes, max distance for connecting two nodes is: %d' % max_gap_distance) 

    clustering_nodes(args, dbo_args, endpoint_args, args.node33_candidate_file, args.node_cluster33_file, max_gap_distance,  'R_end', 'R_end')
    gc.collect()

    clustering_nodes(args, dbo_args, endpoint_args, args.node55_candidate_file, args.node_cluster55_file, max_gap_distance,  'L_end', 'L_end')
    gc.collect()

    clustering_nodes(args, dbo_args, endpoint_args, args.node53_candidate_file, args.node_cluster53_file, max_gap_distance,  'L_end', 'R_end')
    gc.collect()

    clustering_nodes(args, dbo_args, endpoint_args, args.node35_candidate_file, args.node_cluster35_file, max_gap_distance,  'R_end', 'L_end')
    gc.collect()

    if args.rm_temp_files:
        args.temp_file_list = list()
        args.temp_file_list.append(args.node33_candidate_file)
        args.temp_file_list.append(args.node55_candidate_file)
        args.temp_file_list.append(args.node53_candidate_file)
        args.temp_file_list.append(args.node35_candidate_file)
        for temp_file in args.temp_file_list:
            if check_file_exists(temp_file): os.remove(temp_file)

        args.temp_file_list = list()

    return

def run_system_cmd(cmd):

    myprint('Running CMD: %s' % cmd)
    ret = os.system(cmd)
    if ret != 0: 
        myprint('ERROR: Failed to run command: %s' % cmd)
        myprint('Return value is: %d' % ret)
        sys.exit()

    return 

def get_lines_from_file(input_file):

    in_fp = open(input_file, 'r')
    lines = list(in_fp)
    in_fp.close()
    return lines

def read_node_list_file(node_list_file, black_region_key_set, alt_tid_set):
    
    node_list = list()
    node_list_fp = open(node_list_file, 'r')

    num_filtered_nodes = 0
    while 1: 
        line = node_list_fp.readline()
        if not line: break

        line = line.strip().split(tab)
        node = EndpointNode(line[0], line[1], line[2], line[3])

        tid1, pos1 = get_tid_pos_from_key(node.x) 
        tid2, pos2 = get_tid_pos_from_key(node.y) 

        if tid1 in alt_tid_set: continue
        if tid2 in alt_tid_set: continue

        key1 = node.x - node.x % 100
        key2 = node.y - node.y % 100
        if key1 in black_region_key_set or key2 in black_region_key_set:  
            num_filtered_nodes += 1
            continue

        node_list.append(node)

    node_list_fp.close()

    myprint('number of nodes in blacklist region: %s' % num_filtered_nodes)

    return node_list


def read_node_cluster_file(node_cluster_file):
    node_cluster_list = list()
    node_cluster_fp = open(node_cluster_file, 'r')
    while 1:
        line = node_cluster_fp.readline()
        if not line: break
        line = line.strip().split(';')
        one_node_cluster = list()
        for item in line:
            item = item.split(',')
            node = EndpointNode(item[0], item[1], item[2], item[3])
            one_node_cluster.append(node)
        node_cluster_list.append(one_node_cluster)

    node_cluster_fp.close()

    return node_cluster_list

def get_tid_list_from_node_list(node_list):

    tid_set = set()
    for node in node_list:
        tid = int(node.x / FIX_LENGTH)
        tid_set.add(tid)

        tid = int(node.y / FIX_LENGTH)
        tid_set.add(tid)

    tid_list = list(tid_set)
    tid_list.sort()
    return tid_list

def read_black_region_file(black_region_bed_file, chrname2tid):

    black_region_key_set = set()

    black_region_fp = open(black_region_bed_file, 'r')
    while 1:
        line = black_region_fp.readline()
        if not line: break 
        line = line.strip().split(tab)
        chrom = line[0]

        if chrom not in chrname2tid: continue

        start = int(line[1])
        end   = int(line[2])
        tid   = chrname2tid[chrom]

        for i in range(start, end, 100): 
            pos = i - i % 100
            key = tid * FIX_LENGTH + pos

            black_region_key_set.add(key)

    black_region_fp.close()

    return black_region_key_set


def clustering_nodes (args, dbo_args, endpoint_args, node_list_file, output_node_cluster_file, max_gap_distance, endtype1, endtype2):
    
    myprint ('min support fragment pairs is: %d' % args.min_support_fragments)

    myprint ('reading black region bed file')

    black_region_key_set = read_black_region_file(args.black_region_bed_file, args.chrname2tid)

    if args.run_from_begining == True or check_file_exists(output_node_cluster_file) == False:

        myprint('reading node candidate file:%s' % node_list_file)
        node_list = read_node_list_file(node_list_file, black_region_key_set, args.alt_tid_set)
        myprint('number of nodes in node candidate file: %d' % len(node_list))

        output_node_cluster_fp = open(output_node_cluster_file, 'w')
        output_node_cluster_fp.close()

        tid_list = get_tid_list_from_node_list(node_list)
        n_tid = len(tid_list)
        for i in range(0, n_tid):
            for j in range(0, n_tid):
                tid1 = tid_list[i]
                tid2 = tid_list[j]
                chr1 = args.tid2chrname[tid1]
                chr2 = args.tid2chrname[tid2]
                sub_node_list = list()
                for node in node_list:
                    node_tid1 = int(node.x / FIX_LENGTH)
                    node_tid2 = int(node.y / FIX_LENGTH)
                    if tid1 == node_tid1 and tid2 == node_tid2: 
                        sub_node_list.append(node)
                if len(sub_node_list) == 0: continue
                myprint ('clustering nodes for two chrom: %s, %s, number of nodes is: %d' % (chr1, chr2, len(sub_node_list)) )
                clustering_nodes_for_two_chrom(args, dbo_args, endpoint_args, output_node_cluster_file, max_gap_distance, sub_node_list)
    else:
        myprint ('node cluster file existed: %s, skipped clustering' % output_node_cluster_file)
        return

def clustering_nodes_for_two_chrom(args, dbo_args, endpoint_args, output_node_cluster_file, max_gap_distance, node_list):

    if len(node_list) == 0:
        return

    node_list.sort(key = lambda node: node.x)

    batch_node_list = list() 
    output_string = ''
    batch_node_list.append(node_list[0])

    for i in range(1, len(node_list)):
        if len(batch_node_list) < 1000 or node_list[i].x - node_list[i-1].x <= max_gap_distance+1:
            batch_node_list.append(node_list[i])
        else:
            output_string += clustering_nodes_for_one_batch(args, dbo_args, endpoint_args, max_gap_distance, batch_node_list)
            batch_node_list = list()
            batch_node_list.append(node_list[i])
    
    if len(batch_node_list) > 0: 
        output_string += clustering_nodes_for_one_batch(args, dbo_args, endpoint_args, max_gap_distance, batch_node_list)

    output_node_cluster_fp = open(output_node_cluster_file, 'a')
    output_node_cluster_fp.write(output_string) 
    output_node_cluster_fp.close()

    return

def clustering_nodes_for_one_batch(args, dbo_args, endpoint_args, max_gap_distance, node_list):

    if len(node_list) == 0: return  ''

    node_coord_list = list()
    for node in node_list:
        node_coord_list.append((node.x, node.y))

    row = list()
    col = list()
    data = list()

    distance_buffer = max_gap_distance * 1.415
    #myprint ('building KD tree, distance buffer = %d' % distance_buffer)
    tree = cKDTree(node_coord_list, leafsize = 10000)

    #myprint ('searching nearby nodes')
    for i in range(0, len(node_list)):

        node1 = node_list[i]
        index_list = tree.query_ball_point((node1.x, node1.y), distance_buffer)
        if len(index_list) < args.min_support_fragments: continue
        nearby_node_index_list = list()
        for j in index_list:
            node2 = node_list[j]
            if abs(node1.x - node2.x) < max_gap_distance and abs(node1.y - node2.y) < max_gap_distance:
                nearby_node_index_list.append(j)
        if len(nearby_node_index_list) < args.min_support_fragments: continue

        for j in nearby_node_index_list: 
            if i <= j:
                edge = (i, j)
            else:
                edge = (j, i)

            row.append (edge[0])
            col.append (edge[1])
            data.append (1)
        
    #myprint ('finished connecting edges, number of edges is: %d' % len(row))

    n_node = len(node_list)

    #myprint ('getting connected components')
    n_components, label_list, component_node_index_db = get_connected_components(n_node, row, col, data, False, 'weak')
    node_cluster_list = [0] * n_components
    for i in range(0, n_components):
        node_cluster_list[i] = list()
        for index in component_node_index_db[i]:
            node_cluster_list[i].append(node_list[index])
            
    ## output ## 
    total_out_string = ''
    for i in range(0, len(node_cluster_list)):
        node_cluster = node_cluster_list[i]
        if len(node_cluster) < args.min_support_fragments: continue

        output_string = ''
        for j in range(0, len(node_cluster)):
            node = node_cluster[j]
            output_string += node.csvoutput() + ';'

        output_string = output_string.rstrip(';') + endl
        total_out_string += output_string

    del row, col, data, component_node_index_db

    gc.collect()

    return total_out_string 

def get_connected_components(n_node, row, col, data, is_directed = False, connection_type = 'weak'):

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

def convert_node_cluster_to_paired_bk_cand(args, dbo_args, endpoint_args, bcd22_frm_list, node_cluster, max_gap_distance, endtype1, endtype2):

    supp_frm_list1 = list()
    supp_frm_list2 = list()

    bcd22_frag_id_list = list()
    for frm in bcd22_frm_list:
        bcd22_frag_id_list.append(frm.frag_id)

    supp_bcd_list = list()
    for node in node_cluster:
        frm_idx1 = bisect.bisect_left(bcd22_frag_id_list, node.x_frag_id)
        frm_idx2 = bisect.bisect_left(bcd22_frag_id_list, node.y_frag_id)

        supp_frm_list1.append(bcd22_frm_list[frm_idx1])
        supp_frm_list2.append(bcd22_frm_list[frm_idx2])
        supp_bcd_list.append(bcd22_frm_list[frm_idx1].bcd)

    supp_bcd_set = set(supp_bcd_list)
    supp_bcd_string = ''
    for bcd in supp_bcd_set:
        supp_bcd_string += bcd + '|'
    supp_bcd_string.strip('|')

    supp_frm_with_pe_list1 = list()
    supp_frm_with_pe_list2 = list()
    supp_frm_without_pe_list1 = list()
    supp_frm_without_pe_list2 = list()

    for i in range(0, len(supp_frm_list1)):
        if exist_read_pair_support(supp_frm_list1[i], supp_frm_list2[i], endtype1, endtype2):
            supp_frm_with_pe_list1.append(supp_frm_list1[i])
            supp_frm_with_pe_list2.append(supp_frm_list2[i])
        else:
            supp_frm_without_pe_list1.append(supp_frm_list1[i])
            supp_frm_without_pe_list2.append(supp_frm_list2[i])

    bin_size = 50

    xbk_pos, x_total_score, x_total_n_supp, x_withpe_score, x_withpe_n_supp, x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp, ybk_pos, y_total_score, y_total_n_supp, y_withpe_score, y_withpe_n_supp, y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp = predict_breakpoint_position (args, dbo_args, endpoint_args, supp_frm_with_pe_list1, supp_frm_with_pe_list2, supp_frm_without_pe_list1, supp_frm_without_pe_list2, endtype1, endtype2, bin_size)

    xtid, xstart = get_tid_pos_from_key(xbk_pos)
    ytid, ystart = get_tid_pos_from_key(ybk_pos)

    n_supp = min(x_total_n_supp, y_total_n_supp)
    score = x_total_score + y_total_score
    xchr = args.tid2chrname[xtid]
    ychr = args.tid2chrname[ytid]

    svtype = 'UNK'
     
    if xtid == ytid:
        svlength = str(abs(ystart-xstart))
    else:
        svlength = 'NA' 

    info = 'x_total_score=%.2f;x_total_n_supp=%d;x_withpe_score=%.2f;x_withpe_n_supp=%d;x_withoutpe_R_score=%.2f;x_withoutpe_L_score=%.2f;x_withoutpe_n_R_supp=%d;x_withoutpe_n_L_supp=%d;' % (x_total_score, x_total_n_supp, x_withpe_score, x_withpe_n_supp, x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp)
    info += 'y_total_score=%.2f;y_total_n_supp=%d;y_withpe_score=%.2f;y_withpe_n_supp=%d;y_withoutpe_R_score=%.2f;y_withoutpe_L_score=%.2f;y_withoutpe_n_R_supp=%d;y_withoutpe_n_L_supp=%d;' % (y_total_score, y_total_n_supp, y_withpe_score, y_withpe_n_supp, y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp)
    info += 'x_num_supp_frm_withpe=%d;x_num_supp_frm_withoutpe=%d;y_num_supp_frm_withpe=%d;y_num_supp_frm_withoutpe=%d;supp_bcd=%s' % (len(supp_frm_with_pe_list1), len(supp_frm_without_pe_list1), len(supp_frm_with_pe_list2), len(supp_frm_without_pe_list2), supp_bcd_string)

    attr_list = [xchr, xstart, xstart+1, ychr, ystart, ystart+1, svtype, svlength, endtype1, endtype2, n_supp, score, info] 
    paired_bk_cand = PairedBkCand(attr_list) 

    return paired_bk_cand

def predict_breakpoint_position (args, dbo_args, endpoint_args, supp_frm_with_pe_list1, supp_frm_with_pe_list2, supp_frm_without_pe_list1, supp_frm_without_pe_list2, endtype1, endtype2, bin_size):


    x_with_pe_list = list()
    y_with_pe_list = list()

    for i in range(0, len(supp_frm_with_pe_list1)):

        frm1 = supp_frm_with_pe_list1[i]
        frm2 = supp_frm_with_pe_list2[i]

        if endtype1 == 'R_end': 
            x_with_pe_list.append(frm1.key_end())
        else: 
            x_with_pe_list.append(frm1.key_start())

        if endtype2 == 'R_end':
            y_with_pe_list.append(frm2.key_end())
        else:
            y_with_pe_list.append(frm2.key_start())

    xR_without_pe_list = list()
    xL_without_pe_list = list()
    yR_without_pe_list = list()
    yL_without_pe_list = list()
    for i in range(0, len(supp_frm_without_pe_list1)):
        frm1 = supp_frm_without_pe_list1[i]
        frm2 = supp_frm_without_pe_list2[i]
        xR_without_pe_list.append(frm1.key_end())
        xL_without_pe_list.append(frm1.key_start())
        yR_without_pe_list.append(frm2.key_end())
        yL_without_pe_list.append(frm2.key_start())

    if endtype1 == 'R_end': 
        x_total_list = xR_without_pe_list + x_with_pe_list
    else:
        x_total_list = xL_without_pe_list + x_with_pe_list

    if endtype2 == 'R_end': 
        y_total_list = yR_without_pe_list + y_with_pe_list
    else:
        y_total_list = yL_without_pe_list + y_with_pe_list


    xmin = min(x_total_list) - 1000
    xmax = max(x_total_list) + 1000
    ymin = min(y_total_list) - 1000
    ymax = max(y_total_list) + 1000

    mean_gap_withoutpe = 1.0 / args.read_per_bp_genome 
    mean_gap_withpe = 200.0 

    xscore_list = list()
    for x_key in range(xmin, xmax, bin_size): 
        x_withpe_score, x_withpe_n_supp = calculate_withpe_score (x_key, x_with_pe_list, mean_gap_withpe, endtype1)
        x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp = calculate_withoutpe_score(x_key, xR_without_pe_list, xL_without_pe_list, mean_gap_withoutpe)  
        x_total_score = x_withpe_score + x_withoutpe_R_score + x_withoutpe_L_score
        x_total_n_supp = x_withpe_n_supp + x_withoutpe_n_R_supp + x_withoutpe_n_L_supp

        xscore_list.append( (x_key, x_total_score, x_total_n_supp, x_withpe_score, x_withpe_n_supp, x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp) )

    yscore_list = list()
    for y_key in range(ymin, ymax, bin_size):
        y_withpe_score, y_withpe_n_supp = calculate_withpe_score (y_key, y_with_pe_list, mean_gap_withpe, endtype2)
        y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp = calculate_withoutpe_score(y_key, yR_without_pe_list, yL_without_pe_list, mean_gap_withoutpe)  
        y_total_score  = y_withpe_score + y_withoutpe_R_score + y_withoutpe_L_score
        y_total_n_supp = y_withpe_n_supp + y_withoutpe_n_R_supp + y_withoutpe_n_L_supp
        yscore_list.append( (y_key, y_total_score, y_total_n_supp, y_withpe_score, y_withpe_n_supp, y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp) )

    max_xscore_index = get_max_index_from_score_list(xscore_list)
    max_yscore_index = get_max_index_from_score_list(yscore_list)

    xbk_pos, x_total_score, x_total_n_supp, x_withpe_score, x_withpe_n_supp, x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp = xscore_list[max_xscore_index] 
    ybk_pos, y_total_score, y_total_n_supp, y_withpe_score, y_withpe_n_supp, y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp = yscore_list[max_yscore_index]

    return xbk_pos, x_total_score, x_total_n_supp, x_withpe_score, x_withpe_n_supp, x_withoutpe_R_score, x_withoutpe_L_score, x_withoutpe_n_R_supp, x_withoutpe_n_L_supp, ybk_pos, y_total_score, y_total_n_supp, y_withpe_score, y_withpe_n_supp, y_withoutpe_R_score, y_withoutpe_L_score, y_withoutpe_n_R_supp, y_withoutpe_n_L_supp


def get_max_index_from_score_list(score_list):

    max_index = 0
    max_score = score_list[max_index][1]
    for index in range(0, len(score_list)):
        score = score_list[index][1]
        if score > max_score:
            max_score = score
            max_index = index

    return max_index

def calculate_withoutpe_score(bk_pos, eR_without_pe_list, eL_without_pe_list, mean_gap):

    max_score = 1.0
    p = 1.0 / float(mean_gap)
    max_gap = math.log(0.01) / math.log(1-p)

    gap_R_list = list()
    gap_L_list = list()

    for pos in eR_without_pe_list:
        gap = bk_pos - pos
        gap_R_list.append(gap)

    for pos in eL_without_pe_list:
        gap = pos - bk_pos 
        gap_L_list.append(gap)

    total_R_score = 0
    total_L_score = 0

    n_R_supp = 0
    n_L_supp = 0

    for i in range(0, len(gap_R_list)):
        gapR = gap_R_list[i]
        gapL = gap_L_list[i]

        scoreR = convert_gap_to_score(gapR, max_gap, max_score)
        scoreL = convert_gap_to_score(gapL, max_gap, max_score)

        total_R_score += scoreR
        total_L_score += scoreL

        if scoreR > 0: n_R_supp += 1
        if scoreL > 0: n_L_supp += 1

    return total_R_score, total_L_score, n_R_supp, n_L_supp 

def convert_gap_to_score(gap, max_gap, max_score):
    if gap > max_gap: 
        score = 0
    elif gap >= 0 and gap <= max_gap:
        score = float(max_gap - gap) / max_gap * max_score 
    elif gap >= -200.0 and gap < 0: 
        score = max_score/200.0 * gap + max_score 
    else: 
        score = 0   
    return score

def calculate_withpe_score(bk_pos, pos_list, mean_gap, endtype):
    
    if len(pos_list) == 0: return 0.0, 0

    max_score = 5.0
    p = 1.0 / float(mean_gap)
    max_gap = math.log(0.01) / math.log(1-p) 
    n_supp = 0
    if max_gap > 600: max_gap = 600.0 
    if max_gap < 200: max_gap = 200.0 

    gap_list = list()
    for pos in pos_list:
        if endtype == 'R_end':
            gap = bk_pos - pos
        else:
            gap = pos - bk_pos 
        gap_list.append(gap)

    total_score = 0
    for i in range(0, len(gap_list)):
        gap = gap_list[i]
        if gap > max_gap: 
            score = 0
        elif gap >= 0 and gap <= max_gap:
            score = float(max_gap - gap) / max_gap * max_score 
        elif gap >= -100.0 and gap < 0: 
            score = max_score/100.0 * gap + max_score 
        else: 
            score = 0   
        total_score += score
        if score > 0: n_supp += 1

    return total_score, n_supp
    

def find_paired_bk(args, dbo_args, endpoint_args):

    myprint('searching paired breakpoints')

    build_graph_from_fragments(args, dbo_args, endpoint_args)

    get_paired_bk_from_node_clusters(args, dbo_args, endpoint_args)

    return

def get_paired_bk_from_node_clusters(args, dbo_args, endpoint_args):

    max_gap_distance = args.gap_distance_cutoff
        
    out_file = args.bk_cand_pair_file

    out_fp = open(out_file, 'w')

    paired_bk_cand_list33 = get_paired_bk_from1type_node_clusters(args, dbo_args, endpoint_args, 'R_end', 'R_end', args.node_cluster33_file, max_gap_distance, out_fp)  
    paired_bk_cand_list55 = get_paired_bk_from1type_node_clusters(args, dbo_args, endpoint_args, 'L_end', 'L_end', args.node_cluster55_file, max_gap_distance, out_fp)  
    paired_bk_cand_list53 = get_paired_bk_from1type_node_clusters(args, dbo_args, endpoint_args, 'L_end', 'R_end', args.node_cluster53_file, max_gap_distance, out_fp)  
    paired_bk_cand_list35 = get_paired_bk_from1type_node_clusters(args, dbo_args, endpoint_args, 'R_end', 'L_end', args.node_cluster35_file, max_gap_distance, out_fp)  

    out_fp.close()

    return

def get_paired_bk_from1type_node_clusters(args, dbo_args, endpoint_args, endtype1, endtype2, node_cluster_file, max_gap_distance, out_fp):

    node_cluster_list = read_node_cluster_file(node_cluster_file)

    frm_id_set = set()
    for node_cluster in node_cluster_list:
        for node in node_cluster:
            frm_id_set.add(node.x_frag_id)
            frm_id_set.add(node.y_frag_id)

    bcd22_frm_list = extract_frm_from_bcd22_file(endpoint_args.bcd22_file, frm_id_set)

    bcd22_frm_list.sort(key = lambda frm: frm.frag_id)

    myprint ('number of candidate fragments: %d' % len(bcd22_frm_list))

    paired_bk_cand_list = list()

    for node_cluster in node_cluster_list:

        paired_bk_cand = convert_node_cluster_to_paired_bk_cand(args, dbo_args, endpoint_args, bcd22_frm_list, node_cluster, max_gap_distance, endtype1, endtype2)
        paired_bk_cand_list.append(paired_bk_cand)
        out_fp.write(paired_bk_cand.output() + endl)

    return paired_bk_cand_list


if __name__ == '__main__':
    main()
