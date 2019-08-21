#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "tk.h"

#define MAX_READ_PER_CHROM_PER_BARCODE 300000

static int process_one_barcode_in_one_chr(int tid, INT_LIST * read_pos_list, INT_LIST ** wg_bcd_depth_list, int bin_size, INT_LIST * idx_list)
{
	int idx;

	idx_list->size = 0;

	for (int i = 0; i < read_pos_list->size; i++ ) {
		idx = read_pos_list->data_list[i] / bin_size;
		append_int_list(idx_list, idx);
	}

	sort_int_list(idx_list);

	for (int i = 0; i < idx_list->size; i++)
	{
		if (i > 0 && idx_list->data_list[i] == idx_list->data_list[i-1]){ continue; }
		idx = idx_list->data_list[i];
		if (idx >= 0 && idx < wg_bcd_depth_list[tid]->size){
			wg_bcd_depth_list[tid]->data_list[idx] += 1;
		}
	}

	return 0;
}

static int cal_barcode_depth_from_bcd21 (const char * bcd21_file, const char * output_file, const char * faidx_file, int bin_size, int mapq_cutoff)
{
	gzFile bcd21_fp;
	FILE * out_fp;
	int n_chr;
	int n_bin;
	int chr_len;
	char * line;
	char * curr_bcd;
	char * new_bcd;

	int tid;
	int start_pos;
	int end_pos;
	int mapq;
	int avg_pos;
	int curr_tid;
	int new_tid;
	double high_depth, total_depth;
	int flag;
	int hap_type;
	char * read_id;

	INT_LIST ** wg_high_mapq_bcd_depth_list;
	INT_LIST ** wg_total_bcd_depth_list;
	INT_LIST * high_mapq_read_pos_list;
	INT_LIST * total_read_pos_list;

	INT_LIST * idx_list_high;
	INT_LIST * idx_list_total;

	INT_LIST * chr_length_list;

	line = (char *) calloc (LINE_MAX, sizeof(char));
	curr_bcd = (char *) calloc (128, sizeof(char));
	new_bcd = (char *) calloc (128, sizeof(char));
	read_id = (char *) calloc (LINE_MAX, sizeof(char));

	chr_length_list = get_chr_length_from_faidx_file(faidx_file);
	n_chr = chr_length_list->size;    

	wg_high_mapq_bcd_depth_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
	wg_total_bcd_depth_list     = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

	for (int tid = 0; tid < n_chr; tid++)
	{
		chr_len = chr_length_list->data_list[tid]; 
		n_bin = chr_len/bin_size + 2;
		wg_high_mapq_bcd_depth_list[tid] = init_int_list(n_bin); 
		wg_high_mapq_bcd_depth_list[tid]->size = n_bin;

		wg_total_bcd_depth_list[tid] = init_int_list(n_bin); 
		wg_total_bcd_depth_list[tid]->size = n_bin;

	}

	
	high_mapq_read_pos_list = init_int_list( MAX_READ_PER_CHROM_PER_BARCODE);
	total_read_pos_list     = init_int_list( MAX_READ_PER_CHROM_PER_BARCODE);

	idx_list_high   = init_int_list( MAX_READ_PER_CHROM_PER_BARCODE);
	idx_list_total  = init_int_list( MAX_READ_PER_CHROM_PER_BARCODE);


	bcd21_fp = gzopen(bcd21_file, "r");
	if (Z_NULL == bcd21_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bcd21_file);
		exit(1);
	}

	curr_tid = -1;
	int bcd_cnt = 0;
	while (gzgets(bcd21_fp, line, LINE_MAX))
	{
		if (line[0] == '#') { continue; }
		sscanf(line, "%d\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%*s\n", &new_tid, &start_pos, &end_pos, &mapq, new_bcd, &hap_type, read_id, &flag);
		if (flag & (256 + 1024 + 2048) ){ continue; }

		avg_pos = (start_pos + end_pos) / 2;
		if (total_read_pos_list->size == 0)
		{
			curr_tid = new_tid;
			strcpy(curr_bcd, new_bcd);
			append_int_list(total_read_pos_list, avg_pos);
			if (mapq >= mapq_cutoff) {
				append_int_list(high_mapq_read_pos_list, avg_pos);
			}
			continue;
		}
		if (curr_tid != new_tid || strcmp(curr_bcd, new_bcd) != 0) {
			process_one_barcode_in_one_chr(curr_tid, total_read_pos_list, wg_total_bcd_depth_list, bin_size, idx_list_total);
			process_one_barcode_in_one_chr(curr_tid, high_mapq_read_pos_list, wg_high_mapq_bcd_depth_list, bin_size, idx_list_high);
			if (strcmp(curr_bcd, new_bcd) != 0){
				bcd_cnt ++;
				if (bcd_cnt % 1000000 == 0){
					fprintf(stderr, "processed %d barcodes\n", bcd_cnt);
				}
			}
			curr_tid = new_tid;
			strcpy(curr_bcd, new_bcd);
			reset_int_list(total_read_pos_list);
			reset_int_list(high_mapq_read_pos_list);
			append_int_list(total_read_pos_list, avg_pos);
			if (mapq >= mapq_cutoff) {
				append_int_list(high_mapq_read_pos_list, avg_pos);
			}
		}else{
			append_int_list(total_read_pos_list, avg_pos);
			if (mapq >= mapq_cutoff) {
				append_int_list(high_mapq_read_pos_list, avg_pos);
			}
		}
	}
	if (total_read_pos_list->size > 0){
		process_one_barcode_in_one_chr(curr_tid, total_read_pos_list, wg_total_bcd_depth_list, bin_size, idx_list_total);
		process_one_barcode_in_one_chr(curr_tid, high_mapq_read_pos_list, wg_high_mapq_bcd_depth_list, bin_size, idx_list_high);
	}

	out_fp = fopen(output_file, "w");
	if (NULL == out_fp){
		fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
		exit(1);
	}

	fprintf(out_fp, "#tid\tstart_pos\tend_pos\thigh_qual_barcode_depth\ttotal_barcode_depth\n");
	for (int tid = 0; tid < n_chr; tid++)
	{
		n_bin = wg_high_mapq_bcd_depth_list[tid]->size;
		for (int idx = 0; idx < n_bin; idx ++)
		{
			start_pos = idx * bin_size;
			end_pos = start_pos + bin_size;
			high_depth  = (double) wg_high_mapq_bcd_depth_list[tid]->data_list[idx] / (double) bin_size * 100.0;
			total_depth = (double) wg_total_bcd_depth_list[tid]->data_list[idx]     / (double) bin_size * 100.0;

			fprintf(out_fp, "%d\t%d\t%d\t%.2f\t%.2f\n", tid, start_pos, end_pos, high_depth, total_depth); 
		}
	}

	fclose(out_fp);
	gzclose(bcd21_fp);

	for (int tid = 0; tid < n_chr; tid++)
	{
		free_int_list(wg_high_mapq_bcd_depth_list[tid]);
		free_int_list(wg_total_bcd_depth_list[tid]);
	}

	free(wg_high_mapq_bcd_depth_list);
	free(wg_total_bcd_depth_list);

	free_int_list(high_mapq_read_pos_list);
	free_int_list(total_read_pos_list);
	free_int_list(idx_list_high);
	free_int_list(idx_list_total);

	free(line);
	free(curr_bcd);
	free(new_bcd);
	free(read_id);

	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 6){ 
		fprintf (stderr, "Usage: cal_barcode_depth_from_bcd21  <input.bcd21.gz> <output_file> <faidx_file> <bin_size> <mapq_cutoff>\n");
		return 1;
	}

	char * bcd21_file;
	char * output_file;
	char * faidx_file;
	int bin_size;
	int mapq_cutoff;
	int ret;


	bcd21_file = argv[1]; 
	output_file = argv[2];
	faidx_file = argv[3];
	bin_size = atoi(argv[4]);
	mapq_cutoff = atoi(argv[5]);


	if (bin_size < 10){
		fprintf(stderr, "ERROR! bin_size should be at least 10. Your bin_size value is: %d\n", bin_size);
		exit(1);
	}   

	if (mapq_cutoff < 0 || mapq_cutoff > 60){
		fprintf(stderr, "ERROR! mapq_cutoff should be in [0, 60]. Your mapq_cutoff value is: %d\n", mapq_cutoff);
		exit(1);
	}   

	ret = cal_barcode_depth_from_bcd21 (bcd21_file, output_file, faidx_file, bin_size, mapq_cutoff);

	return ret;
}
