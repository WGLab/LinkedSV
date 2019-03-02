#!/usr/bin/env python

from my_utils import *
from bed import *
from bedpe import *
from scipy import spatial
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def main():
    
    args, dbo_args, endpoint_args = parse_user_arguments()

    merge_quantified_calls(args, dbo_args, endpoint_args)

def merge_quantified_calls(args, dbo_args, endpoint_args):
    
    myprint('merging candidate calls')

    quantified_svcall_list = read_object_file(args.quantified_bk_pair_file, QuantifiedBKCand)

    edge_list = list()

    myprint('building edges for candidate calls')

    frm_id_set_list = list()
    for i in range(0, len(quantified_svcall_list)):
        frm_id_set = quantified_svcall_list[i].all_frm_id_set()
        frm_id_set_list.append(frm_id_set)

    for i in range(0, len(quantified_svcall_list)):
        for j in range(i+1, len(quantified_svcall_list)):
            frm_id_set1 = frm_id_set_list[i] 
            frm_id_set2 = frm_id_set_list[j] 
            shared_frm_id_set = frm_id_set1.intersection(frm_id_set2) 
            n_frm_id1 = len(frm_id_set1)
            n_frm_id2 = len(frm_id_set2)
            n_shared_frm = len(shared_frm_id_set)
            if n_shared_frm >= min(n_frm_id1, n_frm_id2) / 2:
                edge_list.append((i, j))
                edge_list.append((j, i))

    row = list()
    col = list()
    data = list()

    for edge in edge_list:
        row.append (edge[0])
        col.append (edge[1])
        data.append (1) 

    n_node = len(quantified_svcall_list)
    myprint('connected components')
    bedpe_csr_matrix = csr_matrix((data, (row, col)), shape=[n_node, n_node])
    n_components, label_list = connected_components(bedpe_csr_matrix, directed = False)

    component_element_db = [0] * n_components

    for i in range(0, len(component_element_db)):
        component_element_db[i] = list()
    # component_element_db[component_id] = list of bedpe index 
    for i in range(0, len(label_list)):
        component_element_db[label_list[i]].append(i)
    
    merged_call_list = list()
    for component_id in range(0, len(component_element_db)):
        bedpe_merge_group = list()
        for index in component_element_db[component_id]: 
            bedpe_merge_group.append(quantified_svcall_list[index])
     
        merged_call = merge1call_group(bedpe_merge_group)
        merged_call_list.append(merged_call)

    merged_call_bedpe_file = args.merged_bedpe_file

    merged_call_bedpe_fp = open(merged_call_bedpe_file, 'w')
    for merged_call in merged_call_list:
        if merged_call.score < 20: continue
        merged_call_bedpe_fp.write(merged_call.output_core() + endl)

    merged_call_bedpe_fp.close()

    return

def merge1call_group(bedpe_merge_group):

    if len(bedpe_merge_group) == 1: return bedpe_merge_group[0]

    integrated_score_list = list()
    for bedpe in bedpe_merge_group:
        #score = bedpe.n_support_reads * int(1e12) + bedpe.score * int(1e4)  + bedpe.type_score
        score = bedpe.score * int(1e4)  + bedpe.type_score
        integrated_score_list.append(score)

    max_score_index = 0
    for i in range(1, len(integrated_score_list)):
        if integrated_score_list[i] > integrated_score_list[max_score_index]: max_score_index = i 

    return bedpe_merge_group[max_score_index]


if __name__ == '__main__':
    main()
