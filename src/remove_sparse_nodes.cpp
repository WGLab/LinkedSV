#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "tk.h"

#define FIX_LENGTH (int64_t)1e10


int usage(FILE * fp)
{
    fprintf (fp, "Usage: remove_sparse_nodes <input_node_file> <output_node_file> <distance_limit> <faidx_file> <min_num_nodes>\n");
    return 0;
}



NODE_LIST *** read_node_file(const char *in_node_file, int n_chr)
{
    FILE * in_node_fp;
    NODE_LIST *** all_nodes_list;
    NODE new_node;
    char * line;
    int node_tid1, node_tid2;

    char *newline_pos;
    all_nodes_list = (NODE_LIST***)calloc(n_chr, sizeof (NODE_LIST**));
    if (all_nodes_list == NULL){
        fprintf(stderr, "ERROR! Failed to allocate memory!\n");
        exit(1);
    }
    for (int i = 0; i < n_chr; i++)
    {
        all_nodes_list[i] = (NODE_LIST**)calloc(n_chr, sizeof(NODE_LIST*));
        if (all_nodes_list[i] == NULL)
        {
            fprintf(stderr, "ERROR! Failed to allocate memory!\n");
            exit(1);
        }
        for (int j = 0; j < n_chr; j++)
        {
            all_nodes_list[i][j] = init_node_list (1000);
        }
    }
    line = (char *) calloc (LINE_MAX, sizeof(char)); 
    in_node_fp = fopen(in_node_file, "r");
    if (NULL == in_node_fp){
        fprintf(stderr, "ERROR! Failed to open file: %s\n", in_node_file);
        exit(1);
    }
    while (fgets(line, LINE_MAX, in_node_fp)) {
        newline_pos = strchr(line, '\n');
        if (newline_pos != NULL){
            *newline_pos = '\0';
        }
        sscanf (line, "%lld\t%lld\t%lld\t%lld\n", &new_node.key1, &new_node.key2, &new_node.frag_id1, &new_node.frag_id2);    
        node_tid1 = (int)(int64_t)(new_node.key1 / FIX_LENGTH);
        node_tid2 = (int)(int64_t)(new_node.key2 / FIX_LENGTH);
        append_node_list (all_nodes_list[node_tid1][node_tid2], new_node);
    }

    fclose(in_node_fp);
    free(line);
    return all_nodes_list;
}

int remove_sparse_nodes(const char *in_node_file, const char *out_node_file, int distance_limit, const char *faidx_file, int min_num_nodes)
{

    INT_LIST * chr_length_list;
    NODE_LIST *** all_nodes_list;
    NODE_LIST * work_nodes_list;
    int n_chr;
    int64_t tid1;
    int64_t tid2;
    int64_t node_tid1, node_tid2;
    int64_t node_pos1, node_pos2;
    int idx1, idx2;
    int n_bin1, n_bin2;
    FILE * out_node_fp;
    NODE node;
    int ** squre_node_count;
    int sum_node_count;

    distance_limit += 1;

    chr_length_list = get_chr_length_from_faidx_file(faidx_file);
    n_chr = chr_length_list->size;

    // fprintf(stderr, "reading input node file: %s\n", in_node_file);
    all_nodes_list = read_node_file(in_node_file, n_chr);

    // fprintf(stderr, "counting nodes in squre regions\n");
    out_node_fp = fopen(out_node_file, "w");
    if (NULL == out_node_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", out_node_file);
        exit(1);
    }

    for (tid1 = 0; tid1 < n_chr; tid1++)
    {
        // fprintf(stderr, "processed %lld chromosome(s)\n", tid1+1);
        for (tid2 = 0; tid2 < n_chr; tid2++)
        {
            
            work_nodes_list = all_nodes_list[tid1][tid2];
            n_bin1 = (int) (chr_length_list->data_list[tid1] / (distance_limit) ) + 1;
            n_bin2 = (int) (chr_length_list->data_list[tid2] / (distance_limit) ) + 1;
            squre_node_count = (int **) calloc (n_bin1, sizeof(int*));
            if (squre_node_count == NULL){
                fprintf(stderr, "ERROR! Failed to allocate memory for squre_node_count!"); 
                exit(1);
            }
            for (int i = 0; i < n_bin1; i++) {
                squre_node_count[i] = (int *) calloc(n_bin2, sizeof(int));
                if (squre_node_count[i] == NULL){
                    fprintf(stderr, "ERROR! Failed to allocate memory for squre_node_count[%d]!", i); 
                    exit(1);
                }
            }

            for (int i = 0; i < work_nodes_list->size; i++)
            {
                node = work_nodes_list->data_list[i];
                node_tid1 = (int64_t)(node.key1 / FIX_LENGTH);
                node_tid2 = (int64_t)(node.key2 / FIX_LENGTH);

                if (node_tid1 != tid1 || node_tid2 != tid2){ continue; }

                node_pos1 = node.key1 % FIX_LENGTH;
                node_pos2 = node.key2 % FIX_LENGTH;
                idx1 = (int) node_pos1 / distance_limit;
                idx2 = (int) node_pos2 / distance_limit;
                if (idx1 >= n_bin1 || idx2 >= n_bin2){
                    fprintf(stderr, "WARNING! Node coordinate larger than chromosome length! Skipped this node.\n");
                }else{
                    squre_node_count[idx1][idx2] += 1;
                }
            }

            for (int i = 0; i < work_nodes_list->size; i++)
            {
                node = work_nodes_list->data_list[i];
                node_tid1 = (int64_t)(node.key1 / FIX_LENGTH);
                node_tid2 = (int64_t)(node.key2 / FIX_LENGTH);
                if (node_tid1 != tid1 || node_tid2 != tid2){ continue; }
                node_pos1 = node.key1 % FIX_LENGTH;
                node_pos2 = node.key2 % FIX_LENGTH;
                idx1 = (int) node_pos1 / distance_limit;
                idx2 = (int) node_pos2 / distance_limit;

                if (idx1 == 0 || idx1 + 1 >= n_bin1) { continue;}
                if (idx2 == 0 || idx2 + 1 >= n_bin2) { continue;}
                sum_node_count = 0;
                for (int a = idx1-1; a <= idx1+1; a++)
                {
                    for (int b = idx2-1; b <= idx2+1; b++) 
                    {
                        sum_node_count += squre_node_count[a][b];
                    }                      
                }
                if (sum_node_count > min_num_nodes){
                    fprintf(out_node_fp, "%lld\t%lld\t%lld\t%lld\n", node.key1, node.key2, node.frag_id1, node.frag_id2);
                }

            }
            for (int i = 0; i < n_bin1; i++) {
                free(squre_node_count[i]);
            }

            free(squre_node_count);
        }
    }
    fclose(out_node_fp);

    return 0;
}

int main (int argc, char * argv[])
{
    char * in_node_file;
    char * out_node_file;
    int distance_limit;
    char * faidx_file;
    int min_num_nodes;

    if (argc < 6){
        usage(stderr);
        return 1;
    }

    in_node_file = argv[1];
    out_node_file = argv[2];
    distance_limit = atoi(argv[3]);
    faidx_file = argv[4];
    min_num_nodes = atoi(argv[5]);

    // fprintf(stderr, "input file is: %s\n", in_node_file);
    // fprintf(stderr, "output file is: %s\n", out_node_file);
    // fprintf(stderr, "distance limit is: %d\n", distance_limit);
    // fprintf(stderr, "faidx file is: %s\n", faidx_file);

    remove_sparse_nodes(in_node_file, out_node_file, distance_limit, faidx_file, min_num_nodes);

    return 0;
}
