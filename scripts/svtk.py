#!/usr/bin/env python

import os
import sys

from scipy.spatial import *
from scipy.sparse import csr_matrix # csr_matrix
from scipy.sparse.csgraph import connected_components # connected_components

try:
    from scripts import my_utils
except ImportError:
    import my_utils

tab  = '\t'
endl = '\n'

class SVCall:
    def __init__(self):
        
        self.chrm1 = ''
        self.pos1 = -1

        self.chrm2 = ''
        self.pos2 = -1
 
        self.sv_type = 'UNK'
        self.filter = ''
        self.sv_size = 0
        self.score = 0.0
        self.method = ''
        self.num_pe_support = 0
        self.num_frm_support = 0
        self.assembled = False
        self.is_precise = False
        self.aux_info = ''
        self.tid1 = -1
        self.tid2 = -1
        
    def key1(self):
        if self.tid1 < 0: 
            my_utils.myprint('ERROR! tid1 < 0')
            sys.exit()
        return self.tid1 * my_utils.FIX_LENGTH + self.pos1
   
    def key2(self):
        if self.tid2 < 0: 
            my_utils.myprint('ERROR! tid2 < 0')
            sys.exit()
        return self.tid2 * my_utils.FIX_LENGTH + self.pos2
  
    def size(self):
        if self.chrm1 == self.chrm2:
            return self.pos2 - self.pos1
        else:
            my_utils.myprint('ERROR! chrm1 != chrm2!')
            sys.exit()
            

    def read_bedpe_line(self, line, chrname2tid_dict = None):
        line = line.strip().split(tab)
        if len(line) < 10: 
            my_utils.myprint ('ERROR! number of columns is less than 10. The line is:')
            my_utils.myprint (tab.join(line))
            sys.exit()
            
        self.chrm1, self.pos1 = line[0:2]
        self.chrm2, self.pos2 = line[3:5]
        self.pos1  = int(self.pos1)
        self.pos2  = int(self.pos2)
    
        self.sv_type, self.sv_id, self.sv_size, self.score, self.filter, self.aux_info = line[6:12]
        self.sv_size = int(self.sv_size)
        self.score   = float(self.score)

        if chrname2tid_dict != None: 
            if self.chrm1 in chrname2tid_dict:
                self.tid1 = chrname2tid_dict[self.chrm1]
            if self.chrm2 in chrname2tid_dict:
                self.tid2 = chrname2tid_dict[self.chrm2]
        
        aux_list = self.aux_info.split(tab)
        for aux in aux_list:
            if aux == 'SVMETHOD=local_assembly':
                self.assembled = True
            if aux == 'PRECISE':
                self.is_precise = True
            if aux == 'IMPRECISE':
                self.is_precise = False

    def output_bedpe_line(self):

        outstring = '%s\t%d\t%d\t%s\t%d\t%d\t' % (self.chrm1, self.pos1, self.pos1+1, self.chrm2, self.pos2, self.pos2+1)      
        outstring += '%s\t%s\t%d\t%.f\t%s\t%s' % (self.sv_type, self.sv_id, self.sv_size, int(self.score), self.filter, self.aux_info)
        
        return outstring
        
        
def read_sv_bedpe_file(in_sv_bedpe_file, chrname2tid_dict = None):
    
    sv_list = list()
    in_sv_bedpe_fp = open(in_sv_bedpe_file, 'r')
    while 1:
        line = in_sv_bedpe_fp.readline()
        if not line: break
        if line[0] == '#': continue
        
        svcall = SVCall()
        svcall.read_bedpe_line(line, chrname2tid_dict)
        sv_list.append(svcall)
        
    in_sv_bedpe_fp.close()
    
    return sv_list
    
def output_sv_list(svcall_list, out_file):

    out_fp = open(out_file, 'w')
    for svcall in svcall_list:
        out_fp.write(svcall.output_bedpe_line() + endl)
        
    out_fp.close()

    return

def remove_redundantsv(svcall_list, overlap_fraction = 0.5, max_distance = 10000, coordinates = 'list1'):

    merged_svcall_list = list()
    coord_list = list()
    distance_buffer = max_distance * 1.415
    for svcall in svcall_list:
        node = (svcall.key1(), svcall.key2())
        coord_list.append(node)

    if len(coord_list) == 0: 
        return svcall_list
    tree = cKDTree(coord_list, leafsize = 10000)
    edge_list = list()

    for i in range(0, len(svcall_list)):
        svcall= svcall_list[i]
        node = (svcall.key1(), svcall.key2())

        index_list = tree.query_ball_point(node, distance_buffer)
        for j in index_list:
            if j == i: continue
            svcall2 = svcall_list[j]
            if is_same_sv(svcall, svcall2, overlap_fraction, max_distance):
                edge_list.append((i, j))

    row = list()
    col = list()
    data = list()
    for edge in edge_list:
        row.append (edge[0])
        col.append (edge[1])
        data.append (1) 
    n_node = len(coord_list)

    n_components, label_list, component_node_index_db = get_connected_components(n_node, row, col, data, False, 'weak')
    sv_group_list = [0] * n_components
    for i in range(0, n_components):
        sv_group_list[i] = list()
        for index in component_node_index_db[i]:
            sv_group_list[i].append(svcall_list[index])

    for i in range(0, len(sv_group_list)):
        merged_svcall_list.append(sv_group_list[i][0])

    return merged_svcall_list
        


def merge2svcallset(svcall_list1, svcall_list2, overlap_fraction = 0.5, max_distance = 10000, coordinates = 'list1'):

    '''
    overlap_fraction: minimal fraction of reciprocal overlap
    coordinates = 'list1'   means use coordiants of svcall_list1
    coordinates = 'list2'   means use coordiants of svcall_list2
    coordinates = 'mean'    means use mean coordiant of svcall_list1 and svcall_list2
    coordinates = 'outer'   means use outer coordiant
    coordinates = 'innter'  means use inner coordiant
    coordinates = 'precise' means use coordiant of the sv call which has precise breakpoint position
    '''

    merged_svcall_list = list()

    if coordinates != 'list1': return merged_svcall_list

    if len(svcall_list1) == 0: return svcall_list2

    if len(svcall_list2) == 0: return svcall_list1


    coord_list1 = list()
    distance_buffer = max_distance * 1.415

    for svcall in svcall_list1:
        node1 = (svcall.key1(), svcall.key2())
        coord_list1.append(node1)

    tree = cKDTree(coord_list1, leafsize = 10000)

    for svcall in svcall_list1:
        merged_svcall_list.append(svcall)

    for i in range(0, len(svcall_list2)):
        sv2 = svcall_list2[i]
        node2 = (sv2.key1(), sv2.key2())

        index_list = tree.query_ball_point(node2, distance_buffer )
        same_sv_index_list = list()
        for j in index_list:
            sv1 = svcall_list1[j]
            if is_same_sv(sv1, sv2, overlap_fraction, max_distance):
                same_sv_index_list.append(j)

        if len(same_sv_index_list) == 0:
            merged_svcall_list.append(sv2)

    return merged_svcall_list


def is_same_sv(sv1, sv2, overlap_fraction, max_distance):

    if sv1.sv_type != sv2.sv_type: return False
    if sv1.chrm1 != sv2.chrm1: return False
    if sv1.chrm2 != sv2.chrm2: return False

    if abs(sv1.pos1-sv2.pos1) > max_distance: return False
    if abs(sv1.pos2-sv2.pos2) > max_distance: return False

    if sv1.chrm1 == sv1.chrm2: 
        ovl_len = min(sv1.pos2, sv2.pos2) - max(sv1.pos1, sv2.pos1)
        ovl_frac = float(ovl_len) / max(sv1.size(), sv2.size())
        if ovl_frac >= overlap_fraction: 
            return True
        else:
            return False
    else:
        return True



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

