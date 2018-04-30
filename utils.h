#ifndef MY_UTILS_H
#define MY_UTILS_H

#include <stdio.h>
#include <stdint.h>

#define LINE_MAX 4096
#define MAX_TID 1<<32
#define DIST_INTRA_CHR INT32_MAX

// CHR_INFO
typedef struct chr_info_t {
    char ** chr_name;
    int * chr_length;
    int n_chr;
} CHR_INFO;

CHR_INFO * get_chr_info (const char * faidx); // read chr info from faidx file

// BC_INFO_SET
typedef struct bc_info1_t{
	int32_t tid;
	int32_t pos;
	int16_t map_rd_len;  // mapped length of the read
    uint8_t mapq;
    uint8_t hap_type;
    uint32_t bcd_int;
	int32_t tid_prev;
	int32_t pos_prev;
	int32_t tid_next;
	int32_t pos_next;
	int32_t flag;
}BC_INFO1;

typedef struct bc_info_set_t{
    int n_item;
    BC_INFO1 * bc_info;
}BC_INFO_SET;

BC_INFO_SET * init_bc_info_set(int n_item);



// position conversion
int64_t tidpos2linearpos(int32_t tid, int32_t pos, CHR_INFO * chr_info); //convert tid,pos to linear pos

int linearpos2tidpos (int64_t linearpos, int32_t * p_tid, int32_t * p_pos, CHR_INFO * chr_info); // convert linear pos to tid,pos


// bcd_seq, bcd_int conversion
inline uint32_t bcdseq2bcdint(const char * bcd_seq);  // convert barcode sequence to uint32_t
inline int bcdint2bcdseq (uint32_t bcd_int, char * bcd_seq); //convert bcd_int to barcode sequence

#endif
