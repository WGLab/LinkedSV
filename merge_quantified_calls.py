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

    #refined_svcall_list = read_object_file(args.refinedbedpe_file, RefinedQuantBKCand)
    quantified_svcall_list = read_object_file(args.quantified_bk_pair_file, QuantifiedBKCand)
    refined_svcall_list = quantified_svcall_list

    edge_list = list()
    myprint('building edges for candidate calls')

    for i in range(0, len(quantified_svcall_list)):
        svcall1 = quantified_svcall_list[i]
        if svcall1.score < 5: continue
        frm_id_set11 = svcall1.frm_id_set1() 
        frm_id_set12 = svcall1.frm_id_set2() 
        if i > 0 and i % 10 == 0: myprint ('processed %d calls' % i)
        for j in range(i+1, len(quantified_svcall_list)):
            svcall2 = quantified_svcall_list[j]
            if svcall2.score < 5: continue
            frm_id_set21 = svcall2.frm_id_set1() 
            frm_id_set22 = svcall2.frm_id_set2() 

            shared_frm_id_set1 = frm_id_set11.intersection(frm_id_set21)
            shared_frm_id_set2 = frm_id_set12.intersection(frm_id_set22)

            if len(shared_frm_id_set1) > 0.66 * float(min(len(frm_id_set11), len(frm_id_set21))) and len(shared_frm_id_set2) > 0.66 * float(min(len(frm_id_set12), len(frm_id_set22))):
                edge_list.append((i, j))
                edge_list.append((j, i))
            else:
                bk1_nearby = 0
                bk2_nearby = 0
                if svcall1.endtype1 == svcall2.endtype1 and abs(svcall1.start1-svcall2.start1) < 5000: bk1_nearby = 1
                if svcall1.endtype1 != svcall2.endtype1 and abs(svcall1.start1-svcall2.start1) < 50000: bk1_nearby = 1

                if svcall1.endtype2 == svcall2.endtype2 and abs(svcall1.start2-svcall2.start2) < 5000: bk2_nearby = 1
                if svcall1.endtype2 != svcall2.endtype2 and abs(svcall1.start2-svcall2.start2) < 50000: bk2_nearby = 1

                if bk1_nearby and bk2_nearby and shared_frm_id_set1 > 3 and shared_frm_id_set2 > 3: 
                    edge_list.append((i, j))
                    edge_list.append((j, i))

    row = list()
    col = list()
    data = list()

    for edge in edge_list:
        row.append (edge[0])
        col.append (edge[1])
        data.append (1) 

    n_node = len(refined_svcall_list)
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
    merged_call_bedpe_file = args.merged_bedpe_file
    merged_call_bedpe_fp = open(merged_call_bedpe_file, 'w')

    for component_id in range(0, len(component_element_db)):
        bedpe_merge_group = list()
        for index in component_element_db[component_id]:
            bedpe_merge_group.append(refined_svcall_list[index])

        retained_index_list = merge1call_group(bedpe_merge_group)
        for index in retained_index_list:
            merged_call = bedpe_merge_group[index]
            if merged_call.score > 5: merged_call_bedpe_fp.write(merged_call.output_core() + endl)

    merged_call_bedpe_fp.close()
 
    return

def merge1call_group(bedpe_merge_group):

    max_score_index = 0
    for i in range(1, len(bedpe_merge_group)):
        if bedpe_merge_group[i].score > bedpe_merge_group[max_score_index].score: 
            max_score_index = i

    max_score = bedpe_merge_group[max_score_index].score
    retained_index_list = list()

    for i in range(0, len(bedpe_merge_group)):
        if bedpe_merge_group[i].score > float(max_score) * 0.75:
            retained_index_list.append(i)

    return retained_index_list 


if __name__ == '__main__':
    main()
