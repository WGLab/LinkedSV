#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


import os
import sys
from scipy.spatial import *
from scipy.sparse import csr_matrix # csr_matrix
from scipy.sparse.csgraph import connected_components # connected_components
import gc

try:
    from scripts.my_utils import *
except ImportError:
    from my_utils import *



tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<weird_reads.txt> <out_file> <faidx_file>'
argc  = 3 



class ShortReadSupport: 
    def __init__(self):
        self.tid1 = self.tid2 = self.start1 = self.end1 = self.start2 = self.end2 = -1
        self.mapq1 = self.mapq2 = 0
        
        self.read_id = self.bcd = '' 
        self.hap_type = self.flag1 = self.flag2 = -1

    def init_from_two_lines(self, line1, line2):
        line1 = line1.strip().split(tab)
        line2 = line2.strip().split(tab)
        if len(line1) < 13:
            my_utils.myprint('ERROR! This line is less than 13 coloumns: %s' % tab.join(line1)) 
            return
        
        if len(line2) < 13:
            my_utils.myprint('ERROR! This line is less than 13 coloumns: %s' % tab.join(line2)) 
            return

        if line1[6] != line2[6]:
            my_utils.myprint('ERROR! line1 and line2 have different read id!')
            my_utils.myprint('line1: %s' % tab.join(line1))
            my_utils.myprint('line2: %s' % tab.join(line2))
            sys.exit()
            return
        
        if line1[0] != line2[0]:
            my_utils.myprint('ERROR! line1 and line2 have different tid!')
            my_utils.myprint('line1: %s' % tab.join(line1))
            my_utils.myprint('line2: %s' % tab.join(line2))
            sys.exit()
            return
        
        if int(line1[1]) > int(line2[1]):
            tmp   = line1
            line1 = line2
            line2 = tmp

        self.tid1, self.start1, self.end1, self.mapq1 = line1[0:4]
        self.tid2, self.start2, self.end2, self.mapq2 = line2[0:4]

        self.tid1    = int(self.tid1)
        self.start1  = int(self.start1)
        self.end1    = int(self.end1)
        self.mapq1   = int(self.mapq1)

        self.tid2    = int(self.tid2)
        self.start2  = int(self.start2)
        self.end2    = int(self.end2)
        self.mapq2   = int(self.mapq2)

        self.flag1 = int(line1[7])
        self.flag2 = int(line2[7])
        self.read_id = line1[6]
        self.bcd = line1[4]
        self.hap_type = int(line1[5])

        if self.start1 > self.start2:
            my_utils.myprint('ERROR! start1 > start2')
            sys.exit()

        return

    def key1(self):
        return my_utils.FIX_LENGTH * self.tid1 + self.end1
        
    def key2(self):
        return my_utils.FIX_LENGTH * self.tid2 + self.start2

    def inner_size(self):
        return self.start2 - self.end1
        
    def output(self):
        outstring = '%s\t%d\t%d\t%d\t%s\t%s\t%d' % (self.tid1, self.end1, self.start2, self.start2-self.end1, self.bcd, self.read_id, self.hap_type)
        return outstring
    def output_info(self):
        outstring = '%s|%d|%d|%d|%s|%s|%d' % (self.tid1, self.end1, self.start2, self.start2-self.end1, self.bcd, self.read_id, self.hap_type)
        return outstring

    def pos1(self):
        return self.end1
    def pos2(self):
        return self.start2

def main():

    if len(arg) < argc:
        print (usage)
        sys.exit()

    in_weird_reads_file = arg.pop(0)
    out_file            = arg.pop(0)
    faidx_file          = arg.pop(0)

    cluster_weird_reads(in_weird_reads_file, out_file, faidx_file)
    return

def cluster_weird_reads(in_weird_reads_file, out_file, faidx_file):

    tid2chrname_list, chrname2tid_dict = my_utils.get_chrnames(faidx_file)
    max_distance = 300
    min_n_short_read_supp = 2
    max_n_short_read_supp = 1000
    min_sv_length = 1000

    my_utils.myprint('reading file: %s' % in_weird_reads_file)
    short_read_support_list35 = read_weird_reads_file(in_weird_reads_file, chrname2tid_dict, min_sv_length)
    my_utils.myprint('finished reading file: %s' % in_weird_reads_file)

    out_fp = open(out_file, 'w')
    out_fp.write('')
    out_fp.close()

    my_utils.myprint('clustering discordant reads')
    cluster_weird_reads1type(short_read_support_list35, out_file, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)

    return

def cluster_weird_reads1type(short_read_support_list, out_file, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp):
    
    coord_list = list()
  
    for short_read_support in short_read_support_list:
        coord_list.append( (short_read_support.key1(), short_read_support.key2()) )

    cluster_one_region(short_read_support_list, coord_list, out_file, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp)
    gc.collect()

    return

def cluster_one_region(short_read_support_list, coord_list, out_file, min_n_short_read_supp, max_distance, tid2chrname_list, chrname2tid_dict, max_n_short_read_supp):

    if len(coord_list) < 1: return
    edge_list = list()

    distance_buffer = max_distance * 1.415
    tree = cKDTree(coord_list, leafsize = 10000)

    for i in range(0, len(short_read_support_list)):

        if i > 0 and i % 100000 == 0: my_utils.myprint ('finished searching for %d weird reads' % i)

        node1 = (short_read_support_list[i].key1(), short_read_support_list[i].key2())
        index_list = tree.query_ball_point( node1, distance_buffer )

        if len(index_list) > max_n_short_read_supp: continue

        nearby_node_index_list = list()
        for j in index_list:
            if i == j: continue
            node2 = (short_read_support_list[j].key1(), short_read_support_list[j].key2())
            if abs(node1[0] - node2[0]) < max_distance and abs(node1[1] - node2[1]) < max_distance:
                nearby_node_index_list.append(j)

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

    my_utils.myprint ('get connected components')
    n_components, label_list, component_node_index_db = get_connected_components(n_node, row, col, data, False, 'weak')
    node_cluster_list = [0] * n_components
    for i in range(0, n_components):
        node_cluster_list[i] = list()
        for index in component_node_index_db[i]:
            node_cluster_list[i].append(short_read_support_list[index])

    my_utils.myprint ('output clusters of weird reads')
    out_fp = open(out_file, 'w')
    for i in range(0, len(node_cluster_list)): # for i-th cluster
        node_cluster = node_cluster_list[i]
        if len(node_cluster) < min_n_short_read_supp: continue
        if len(node_cluster) > max_n_short_read_supp: continue
        mean_start_pos = mean_end_pos = 0
        hap_type_cnt = [0] * 3
        output_info_string = 'SVTYPE=DEL'
        for j in range(0, len(node_cluster)):
            short_read_support = node_cluster[j]
            output_info_string += ';' + short_read_support.output_info() 
            mean_start_pos += short_read_support.pos1()
            mean_end_pos   += short_read_support.pos2()
            hap_type_cnt[short_read_support.hap_type] += 1

        num_pe_supp = len(node_cluster)
        mean_start_pos = int( 0.5 + (float(mean_start_pos)) / num_pe_supp)
        mean_end_pos   = int( 0.5 + (float(mean_end_pos))   / num_pe_supp)
        
        tid = node_cluster[0].tid1
        chrom = tid2chrname_list[tid]
        if len(node_cluster) >= 5:
            flt = 'PASS'
        else:
            flt = 'LowQual'
        sv_size = mean_end_pos - mean_start_pos
        sv_type = 'DEL'
        out_fp.write('%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n' % (chrom, mean_start_pos, mean_start_pos+1, chrom, mean_end_pos, mean_end_pos+1, sv_type, flt, sv_size, num_pe_supp, hap_type_cnt[0], hap_type_cnt[1], hap_type_cnt[2], output_info_string))

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

def get_weird_readname_dict (chrname2tid_dict, weird_reads_file):

    weird_reads_fp = open(weird_reads_file, 'r')
    weird_readname_dict = dict()
    while 1:
        line = weird_reads_fp.readline()
        if not line: break
        line = line.strip().split(tab)
        if len(line) < 10: continue
        tid1 = chrname2tid_dict[line[0]]
        tid2 = chrname2tid_dict[line[2]]
        short_read_support = ShortReadSupport(line + [tid1, tid2])

        readname = short_read_support.aligned_read1.split('@')[0]
        weird_readname_dict[readname] = short_read_support

    weird_reads_fp.close()
    return weird_readname_dict


def read_weird_reads_file(in_weird_reads_file, chrname2tid_dict, min_sv_length):

    short_read_support_list35 = list()

    in_weird_reads_fp = open(in_weird_reads_file, 'r')

    while 1:
        line1 = in_weird_reads_fp.readline()
        line2 = in_weird_reads_fp.readline()
        if not line1: break
        if not line2: break

        short_read_support = ShortReadSupport()
        short_read_support.init_from_two_lines(line1, line2)

        if short_read_support.inner_size() > 45000: continue 
        if short_read_support.inner_size() < min_sv_length - 200: continue
        short_read_support_list35.append(short_read_support)

    in_weird_reads_fp.close()


    return short_read_support_list35

    


if __name__ == '__main__':
    main()
