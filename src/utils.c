#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "utils.h"

static int debug = 1;
BC_INFO_SET * init_bc_info_set(int n_item)
{
    BC_INFO_SET * bc_info_set;
    bc_info_set = (BC_INFO_SET *) calloc (1, sizeof(BC_INFO_SET));
    if (bc_info_set == NULL){
        fprintf (stderr, "[%s] ERROR! No enough memory!\n", __func__);
        return NULL;
    }
    bc_info_set->n_item = n_item;
    if (debug){
        fprintf (stderr, "number of items in bc_info_set: %d\n", bc_info_set->n_item);
        //fprintf (stderr, "sizeof(BC_INFO1): %d\n", sizeof(BC_INFO1));
    }
    bc_info_set->bc_info = (BC_INFO1 *) calloc(bc_info_set->n_item , sizeof(BC_INFO1));
    if (bc_info_set->bc_info == NULL){
        fprintf (stderr, "[%s] ERROR! No enough memory!\n", __func__);
        free (bc_info_set);
        return NULL;
    }
    return bc_info_set;
}

CHR_INFO * get_chr_info (const char * faidx) {
    FILE * fp ;
    char * line;
    CHR_INFO * chr_info;
    int i;
    int total_len, a, b;

    line     = (char *) calloc(LINE_MAX, sizeof(char));
    chr_info = (CHR_INFO *) calloc(1, sizeof(CHR_INFO));

    fp = fopen(faidx, "r");
    if (fp == NULL) {
        fprintf (stderr, "[%s]: ERROR! Can not open file: %s \n", __func__, faidx);
        return NULL;
    }

    chr_info->n_chr = 0;
    // count how many lines in the faidx file
    while (fgets(line, LINE_MAX, fp)){
        if (line[0] != '\n'){
            chr_info->n_chr++;
        }
    }

    // alloc memory for chr_length and chr_name
    chr_info->chr_length = (int *) calloc(chr_info->n_chr, sizeof(int));
    chr_info->chr_name = (char **) calloc(chr_info->n_chr, sizeof(char*));
    for (i = 0; i < chr_info->n_chr; i++)
    {
        chr_info->chr_name[i] = (char *)calloc(LINE_MAX, sizeof(char));
    }

    fseek(fp, 0, SEEK_SET);
    i = 0;

    // read from the faidx file
    while (fgets(line, LINE_MAX, fp)){
        sscanf(line, "%s %d %d %d %d", chr_info->chr_name[i], &chr_info->chr_length[i], &total_len, &a, &b);
        i++;
    }

    fclose(fp);
    free (line);
    fprintf (stderr, "[%s] Finish reading chr info from file: %s.\n", __func__, faidx);
    return chr_info;
}

inline int64_t tidpos2linearpos(int32_t tid, int32_t pos, CHR_INFO * chr_info) {
    if (chr_info == NULL){
        fprintf (stderr, "[%s]: ERROR! chr_info is NULL.\n", __func__);
        return -1;
    }
    int i;
    int64_t linearpos;

    i = 0;
    linearpos = pos;
    for (i = 0; i < tid; i++) {
        linearpos += chr_info->chr_length[i];
    }
    return linearpos;
}

inline int linearpos2tidpos(int64_t linearpos, int32_t * p_tid, int32_t * p_pos, CHR_INFO * chr_info)
{
    if (chr_info == NULL){
        fprintf (stderr, "[%s]: ERROR! chr_info is NULL.\n", __func__);
        *p_tid = -1; 
        *p_pos = -1;
        return 1;
    }

    int32_t i = 0;
    while (linearpos > chr_info->chr_length[i])
    {
        linearpos -= chr_info->chr_length[i];
        i++;
    }
    *p_tid = i;
    *p_pos = (int32_t) linearpos;

    return 0;
}

inline uint32_t bcdseq2bcdint(const char * bcd_seq)
{
    uint32_t bcd_int;
    int i;
    uint8_t base2int[256] = {10};
    base2int['A'] = 0;
    base2int['C'] = 1;
    base2int['G'] = 2;
    base2int['T'] = 3;
    if (bcd_seq == NULL){
        fprintf (stderr, "[%s]: ERROR! barcode sequence is NULL. \n", __func__);
        return 0;
    }
    bcd_int = base2int[bcd_seq[0]];
    for (i = 1; i < 16; i++)
    {
        bcd_int = bcd_int <<2 | base2int[bcd_seq[i]];  
        if (base2int[bcd_seq[i]] > 3) {
            return UINT32_MAX;
        }
    }

    return bcd_int;
}

inline int bcdint2bcdseq (uint32_t bcd_int, char * bcd_seq) // bcd_seq must have at least 19 bytes spaces 
{
    int i;
    if (NULL == bcd_seq){
        fprintf (stderr, "[%s] ERROR! bcd_seq is NULL!\n", __func__);
        return 1;
    }
    char int2base[4] = {'A', 'C', 'G', 'T'};
    bcd_seq[16] = '-';
    bcd_seq[17] = '1';
    bcd_seq[18] = '\0';
    for (i = 15; i >= 0; i--)
    {
        bcd_seq[i] = int2base[(bcd_int&3)];
        bcd_int = bcd_int >> 2;
    }
    return 0;
}






