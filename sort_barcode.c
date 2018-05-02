#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "utils.h"


static BC_INFO_SET * read_bc_info(char * bcd_file, CHR_INFO * chr_info)
{
    FILE * fp = NULL;
    char * line = NULL;
    BC_INFO_SET * bc_info_set = NULL;
    char * bcd_seq = NULL;
    int i;
    int tid;
    int pos;
    int endpos;
    int mapq;
    int hap_type;
    int n_item;

    line  = (char *) calloc(LINE_MAX, sizeof(char));
    bcd_seq = (char *) calloc(LINE_MAX, sizeof(char));

    fp = fopen(bcd_file, "r");
    if (fp == NULL) {
        fprintf (stderr, "[%s]: ERROR! Can not open file: %s \n", __func__, bcd_file);
        return NULL;
    }  

    fprintf (stderr, "[%s] Start reading bcd file: %s \n", __func__, bcd_file);
    n_item = 0;
    while (fgets(line, LINE_MAX, fp)){
        if (line[0] != '\n' && line[0] != '\0'){
            n_item++;
        }
    }
    
    fseek(fp, 0, SEEK_SET);
    bc_info_set = init_bc_info_set(n_item);
    i = 0;
    while (fgets(line, LINE_MAX, fp)){
        sscanf(line, "%d\t%d\t%d\t%d\t%s\t%d\n", &tid, &pos, &endpos, &mapq, bcd_seq, &hap_type);
        bc_info_set->bc_info[i].tid = tid; 
        bc_info_set->bc_info[i].pos = pos; 
        bc_info_set->bc_info[i].map_rd_len = (int16_t)(endpos - pos); 
        bc_info_set->bc_info[i].bcd_int = bcdseq2bcdint(bcd_seq);
        bc_info_set->bc_info[i].mapq = mapq;
        bc_info_set->bc_info[i].hap_type = hap_type;
        bc_info_set->bc_info[i].tid_prev = 0;
        bc_info_set->bc_info[i].pos_prev = -1;
        bc_info_set->bc_info[i].tid_next = INT32_MAX;
        bc_info_set->bc_info[i].pos_next = INT32_MAX;
        i++;
    }

    fclose(fp);
    free(line);
    free(bcd_seq);
    fprintf (stderr, "[%s] Finished reading from bcd file. %d items read.\n", __func__, bc_info_set->n_item);
    return bc_info_set;
}

static inline int cmp_bc_seq (const void * a, const void * b)
{
    BC_INFO1 * bc_info1;
    BC_INFO1 * bc_info2;

    bc_info1 = (BC_INFO1 *) a;
    bc_info2 = (BC_INFO1 *) b;

    if (bc_info1->bcd_int < bc_info2->bcd_int){
        return -1;
    }else if (bc_info1->bcd_int > bc_info2->bcd_int){
        return 1;
    }else{
        if (bc_info1->tid < bc_info2->tid){
            return -1;
        }else if (bc_info1->tid > bc_info2->tid) {
            return 1; 
        }else{
            if (bc_info1->pos < bc_info2->pos){
                return -1;
            }else if (bc_info1->pos > bc_info2->pos){
                return 1;
            }else{
                return 0;
            }
        }
    }
}

static inline int cmp_map_pos (const void * a, const void * b)
{
    BC_INFO1 * bc_info1;
    BC_INFO1 * bc_info2;

    bc_info1 = (BC_INFO1 *) a;
    bc_info2 = (BC_INFO1 *) b;

    if (bc_info1->tid < bc_info2->tid){
        return -1;
    }else if (bc_info1->tid > bc_info2->tid) {
        return 1; 
    }else{
        if (bc_info1->pos < bc_info2->pos){
            return -1;
        }else if (bc_info1->pos > bc_info2->pos){
            return 1;
        }else{
            return 0;
        }
    }
}
static int output_files(BC_INFO_SET * bc_info_set, char * bcd21_file, CHR_INFO * chr_info)
{
    FILE * fp;
    int i;
    int mapq, hap_type;
    int32_t dist_prev, dist_next;
    int32_t endpos;
    char * bcd_seq;
    BC_INFO1 * bc_info;

    if (NULL == bc_info_set || NULL == bc_info_set->bc_info)
    {
        fprintf (stderr, "[%s] ERROR! bc_info_set is NULL!\n", __func__);
        return 1;
    }

    bc_info = bc_info_set->bc_info;
    fp = fopen(bcd21_file, "w");
    if (NULL == fp){
        fprintf (stderr, "[%s] ERROR! Can not open file for write: \n", __func__, bcd21_file);
        return 1;
    }
    
    bcd_seq = (char *) calloc (LINE_MAX, sizeof(char));
    mapq = hap_type = 0;
    dist_prev = dist_next = 0;
    fprintf (stderr, "[%s] Start output file: %s .\n", __func__, bcd21_file);
    for (i = 0; i < bc_info_set->n_item; i++){
        mapq = (int) bc_info[i].mapq;
        hap_type = (int) bc_info[i].hap_type;
        if (bc_info[i].tid == bc_info[i].tid_prev){
            dist_prev = bc_info[i].pos - bc_info[i].pos_prev;
        }else{
            dist_prev = DIST_INTRA_CHR;
        }
        if (bc_info[i].tid == bc_info[i].tid_next){
            dist_next = bc_info[i].pos_next - bc_info[i].pos;
        }else{
            dist_next = DIST_INTRA_CHR;
        }
        bcdint2bcdseq (bc_info[i].bcd_int, bcd_seq);
        endpos = bc_info[i].pos + (int32_t) bc_info[i].map_rd_len;
        fprintf (fp, "%d\t%d\t%d\t%d\t%s\t%d\t%d:%d\t%d:%d\t%d\t%d\n", bc_info[i].tid, bc_info[i].pos, endpos, mapq, bcd_seq, hap_type, bc_info[i].tid_prev, bc_info[i].pos_prev, bc_info[i].tid_next, bc_info[i].pos_next, dist_prev, dist_next);
    }

    return 0;

}

static int get_adjacent_position(BC_INFO_SET * bc_info_set, CHR_INFO *  chr_info)
{
    int i;
    uint32_t curr_bcd_int;
    BC_INFO1 * bc_info;

    if (NULL == bc_info_set || NULL == bc_info_set->bc_info)
    {
        fprintf (stderr, "[%s] ERROR! bc_info_set is NULL!\n", __func__);
        return 1;
    }

    fprintf (stderr, "[%s] get adjacent position. \n", __func__);
    bc_info = bc_info_set->bc_info;
    i = 0;
    curr_bcd_int = bc_info[i].bcd_int;
    for (i = 1; i < bc_info_set->n_item; i++)
    {
        if (bc_info[i].bcd_int == curr_bcd_int){
            bc_info[i].tid_prev = bc_info[i-1].tid;
            bc_info[i].pos_prev = bc_info[i-1].pos;
            bc_info[i-1].tid_next = bc_info[i].tid;
            bc_info[i-1].pos_next = bc_info[i].pos;
        }else{
            curr_bcd_int = bc_info[i].bcd_int;
        }
    } 
    return 0;
}

static int usage(FILE * fp) 
{
    fprintf (fp, "Usage: sort_barcode  <in.bcd> <out_prefix> <faidx> \n");
    return 0;
}


int main (int argc, char * argv[])
{
    char * bcd_file;
    char * bcd21_file;
    char * faidx;
    char * out_prefix;

    if (argc < 4){
        usage(stderr);
        return 1;
    }

    bcd_file = argv[1];
    out_prefix = argv[2];
    faidx = argv[3];

    CHR_INFO * chr_info;
    BC_INFO_SET * bc_info_set;
    bcd21_file = (char *) calloc (LINE_MAX, sizeof(char));
    sprintf (bcd21_file, "%s.bcd21", out_prefix); // sorted by barcode

    chr_info = get_chr_info (faidx);

    bc_info_set =  read_bc_info(bcd_file, chr_info);

    fprintf (stderr, "[%s] Start sorting ... \n", __func__);
    qsort(bc_info_set->bc_info, bc_info_set->n_item, sizeof(BC_INFO1), cmp_bc_seq); // sort by barcode seq and read pos;
    fprintf (stderr, "[%s] Finished sorting ... \n", __func__);

    get_adjacent_position(bc_info_set, chr_info);
    
    output_files (bc_info_set, bcd21_file, chr_info);
    


    return 0;

}
