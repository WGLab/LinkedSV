#!/usr/bin/env python

import os
import sys
from scipy.spatial import *
from scipy.sparse import csr_matrix # csr_matrix
from scipy.sparse.csgraph import connected_components # connected_components
from my_utils import *
import gc


tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<weird_reads.txt> <out_file> <faidx_file>'
argc  = 2 


class ShortReadSupport: 
    def __init__(self, attr_list) :
        self.chrm1, self.pos1, self.chrm2, self.pos2, self.endtype1, self.endtype2, self.mapq, self.supp_type, self.aligned_read1, self.aligned_read2, self.tid1, self.tid2 = attr_list[0:12]

        self.pos1 = int(self.pos1)
        self.pos2 = int(self.pos2)
        self.mapq = int(self.mapq)
        self.tid1 = int(self.tid1)
        self.tid2 = int(self.tid2)

        self.key1  = FIX_LENGTH * self.tid1 + self.pos1
        self.key2  = FIX_LENGTH * self.tid2 + self.pos2

    def output(self):
        outstring = '%s\t%d\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s' % (self.chrm1, self.pos1, self.chrm2, self.pos2, self.endtype1, self.endtype2, self.mapq, self.supp_type, self.aligned_read1, self.aligned_read2)
        return outstring

def main():

    if len(arg) < argc:
        print usage
        sys.exit()

    in_weird_reads_file = arg.pop(0)
    out_file            = arg.pop(0)
    faidx_file          = arg.pop(0)


    tid2chrname_list, chrname2tid_dict = get_chrnames(faidx_file)
    max_distance = 200
    min_n_short_read_supp = 2
    max_n_short_read_supp = 100 

    cluster_weired_reads(in_weird_reads_file, out_file, tid2chrname_list, chrname2tid_dict, max_distance, min_n_short_read_supp, max_n_short_read_supp)


def cluster_weired_reads(in_weird_reads_file, out_file, tid2chrname_list, chrname2tid_dict, max_distance, min_n_short_read_supp, max_n_short_read_supp):


    myprint ('reading weird reads file: %s' % in_weird_reads_file)
    short_read_support_list33, short_read_support_list55, short_read_support_list53, short_read_support_list35 = read_weird_reads_file(in_weird_reads_file, chrname2tid_dict)

    out_fp = open(out_file, 'w')

    cluster_weired_reads1type(short_read_support_list33, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)
    cluster_weired_reads1type(short_read_support_list55, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)
    cluster_weired_reads1type(short_read_support_list53, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)
    cluster_weired_reads1type(short_read_support_list35, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)

    out_fp.close()

    return

def cluster_weired_reads1type(short_read_support_list, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp):
    

    max_n_chr = len(tid2chrname_list)
    all_coord_lists = [0] * max_n_chr 
    chr_short_read_support_lists = [0] * max_n_chr
    for tid1 in range(0, max_n_chr):
        all_coord_lists[tid1] = [0] * max_n_chr
        chr_short_read_support_lists[tid1] = [0] * max_n_chr
        for tid2 in range(0, max_n_chr):
            all_coord_lists[tid1][tid2] = list() 
            chr_short_read_support_lists[tid1][tid2] = list() 

    for short_read_support in short_read_support_list:
        tid1 = short_read_support.tid1 
        tid2 = short_read_support.tid2
        all_coord_lists[tid1][tid2].append( (short_read_support.key1, short_read_support.key2) )
        chr_short_read_support_lists[tid1][tid2].append(short_read_support) 

    del short_read_support_list

    for tid1 in range(0, max_n_chr):
        for tid2 in range(tid1, max_n_chr):
            if len(all_coord_lists[tid1][tid2]) < 1: continue 
            myprint('tid1=%d, tid2=%d, %s, %s, %d reads' % (tid1, tid2, chr_short_read_support_lists[tid1][tid2][0].endtype1, chr_short_read_support_lists[tid1][tid2][0].endtype2, len(all_coord_lists[tid1][tid2]) ) )
            cluster_one_region(chr_short_read_support_lists[tid1][tid2], all_coord_lists[tid1][tid2], out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)
            chr_short_read_support_lists[tid1][tid2] = list()
            all_coord_lists[tid1][tid2] = list()
            gc.collect()

def cluster_one_region(short_read_support_list, coord_list, out_fp, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp):


    if len(coord_list) < 1: return
    edge_list = list() 
    distance_buffer = max_distance * 1.415
    tree = cKDTree(coord_list, leafsize = 10000)

    for i in range(0, len(short_read_support_list)):

        if i > 0 and i % 10000 == 0: myprint ('finished searching for %d weird reads' % i)

        node1 = (short_read_support_list[i].key1, short_read_support_list[i].key2)
        index_list = tree.query_ball_point( node1, distance_buffer )

        if len(index_list) < min_n_short_read_supp: continue
        if len(index_list) > max_n_short_read_supp: continue

        nearby_node_index_list = list()
        for j in index_list:
            node2 = (short_read_support_list[j].key1, short_read_support_list[j].key2)
            if abs(node1[0] - node2[0]) < max_distance and abs(node1[1] - node2[1]) < max_distance:
                nearby_node_index_list.append(j)
            if len(nearby_node_index_list) < min_n_short_read_supp: continue
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

    n_node = len(short_read_support_list)

    myprint ('get connected components')
    n_components, label_list, component_node_index_db = get_connected_components(n_node, row, col, data, False, 'weak')
    node_cluster_list = [0] * n_components
    for i in range(0, n_components):
        node_cluster_list[i] = list()
        for index in component_node_index_db[i]:
            node_cluster_list[i].append(short_read_support_list[index])

    myprint ('output clusters of weird reads')
    for i in range(0, len(node_cluster_list)):
        node_cluster = node_cluster_list[i]
        if len(node_cluster) < min_n_short_read_supp: continue
        if len(node_cluster) > max_n_short_read_supp: continue
        output_string = '#cluster\t%d\n' % (len(node_cluster))
        for j in range(0, len(node_cluster)):
            node = node_cluster[j]
            output_string += node.output() + endl 
        out_fp.write(output_string)

    del edge_list, row, col, data, component_node_index_db, label_list, node_cluster_list
    gc.collect()

    return

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

def get_weired_readname_dict (chrname2tid_dict, weired_reads_file):

    weired_reads_fp = open(weired_reads_file, 'r')
    weired_readname_dict = dict()
    while 1:
        line = weired_reads_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        if len(line) < 10: continue
        tid1 = chrname2tid_dict[line[0]]
        tid2 = chrname2tid_dict[line[2]]
        short_read_support = ShortReadSupport(line + [tid1, tid2])

        readname = short_read_support.aligned_read1.split('@')[0]
        weired_readname_dict[readname] = short_read_support

    weired_reads_fp.close()
    return weired_readname_dict


def read_weird_reads_file(in_weird_reads_file, chrname2tid_dict):

    short_read_support_list33 = list()
    short_read_support_list55 = list()
    short_read_support_list53 = list()
    short_read_support_list35 = list()

    in_weird_reads_fp = open(in_weird_reads_file, 'r')

    while 1:
        line = in_weird_reads_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        if len(line) < 10:
            continue
        tid1 = chrname2tid_dict[line[0]]
        tid2 = chrname2tid_dict[line[2]]

        short_read_support = ShortReadSupport(line + [tid1, tid2])
        if abs(short_read_support.key1 - short_read_support.key2) > 50000: continue 

        if short_read_support.endtype1 == '3p_end' and short_read_support.endtype2 == '3p_end': 
            short_read_support_list33.append(short_read_support)

        elif short_read_support.endtype1 == '5p_end' and short_read_support.endtype2 == '5p_end': 
            short_read_support_list55.append(short_read_support)

        elif short_read_support.endtype1 == '5p_end' and short_read_support.endtype2 == '3p_end': 
            short_read_support_list53.append(short_read_support)

        elif short_read_support.endtype1 == '3p_end' and short_read_support.endtype2 == '5p_end': 
            short_read_support_list35.append(short_read_support)

    in_weird_reads_fp.close()


    return short_read_support_list33, short_read_support_list55, short_read_support_list53, short_read_support_list35

    


if __name__ == '__main__':
    main()
