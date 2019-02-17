#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "tk.h"

static int convert_read_pos_to_segm(INT_LIST * read_pos_list, INT_LIST * seg_start_list, INT_LIST * seg_end_list, int shift_size, int win_size)
{
	int start; 
	int end;
	int curr_pos;
	int n_bin_per_win;

	seg_start_list->size = 0;
	seg_end_list->size = 0;

	if (read_pos_list->size == 0){
		return 0;
	}

	n_bin_per_win = win_size/shift_size;
	start = read_pos_list->data_list[0];
	end = start;
	for (int i = 1; i < read_pos_list->size; i++)
	{
		curr_pos = read_pos_list->data_list[i];
		if (curr_pos/shift_size - end/shift_size < n_bin_per_win){
			end = curr_pos;
		}else{
			append_int_list(seg_start_list, start);
			append_int_list(seg_end_list, end);
			start = curr_pos;
			end = start;
		}
	}
	append_int_list(seg_start_list, start);
	append_int_list(seg_end_list, end);
	return 0;
}

static int process_one_barcode_in_one_chr(int tid, INT_LIST * read_pos_list, INT_LIST ** wg_win1_cnt_list, INT_LIST ** wg_win2_cnt_list, INT_LIST ** wg_ovl_cnt_list, int shift_size, int win_size, INT_LIST * chr_length_list, INT_LIST * seg_start_list, INT_LIST * seg_end_list, INT_LIST * idx_list_win1, INT_LIST * idx_list_win2, INT_LIST * idx_list_ovl)
{
	int n_bin_per_win;
	int n_chr;
	int n_segm;
	int start_bin_idx;
	int end_bin_idx;
	int win1_min_idx;
	int win2_min_idx;
	int win1_max_idx;
	int win2_max_idx;
	int start_pos;
	int end_pos;
	int idx;

	// preprocess INT_LIST that are reused multiple times

	seg_start_list->size = 0;
	seg_end_list->size = 0;
	idx_list_win1->size = 0;
	idx_list_win2->size = 0;
	idx_list_ovl->size = 0;


	n_bin_per_win = win_size/shift_size;
	n_chr = chr_length_list->size;

	convert_read_pos_to_segm(read_pos_list, seg_start_list, seg_end_list, shift_size, win_size);
	n_segm = seg_start_list->size;

	for (int i = 0; i < n_segm; i++)
	{
		start_pos = seg_start_list->data_list[i];
		end_pos   = seg_end_list->data_list[i];
		start_bin_idx = start_pos / shift_size;
		end_bin_idx = end_pos / shift_size;
		win2_min_idx = start_bin_idx - (n_bin_per_win-1);
		win2_max_idx = end_bin_idx;

		win1_min_idx = start_bin_idx + 1;
		win1_max_idx = end_bin_idx + 1 + (n_bin_per_win-1);

		for (int idx = win1_min_idx; idx < win1_max_idx+1; idx++) {
			append_int_list(idx_list_win1, idx);
		}
		for (int idx = win2_min_idx; idx < win2_max_idx+1; idx++) {
			append_int_list(idx_list_win2, idx);
		}
	}

	sort_int_list(idx_list_win1);
	sort_int_list(idx_list_win2);

	int i = 0;
	int j = 0;

	while (i < idx_list_win1->size && j < idx_list_win2->size)
	{
		if (idx_list_win1->data_list[i] > idx_list_win2->data_list[j]) {
			j++;
		}else if (idx_list_win1->data_list[i] < idx_list_win2->data_list[j]) {
			i++;
		}else{
			append_int_list(idx_list_ovl, idx_list_win1->data_list[i]);
			i++;
			j++;
		}
	}

	for (int i = 0; i < idx_list_win1->size; i++)
	{
		if (i > 0 && idx_list_win1->data_list[i] == idx_list_win1->data_list[i-1]){ continue; }
		idx = idx_list_win1->data_list[i];
		if (idx >= 0 && idx < wg_win1_cnt_list[tid]->size){
			wg_win1_cnt_list[tid]->data_list[idx] += 1;
		}
	}
	for (int i = 0; i < idx_list_win2->size; i++)
	{
		if (i > 0 && idx_list_win2->data_list[i] == idx_list_win2->data_list[i-1]){ continue; }
		idx = idx_list_win2->data_list[i];
		if (idx >= 0 && idx < wg_win2_cnt_list[tid]->size){
			wg_win2_cnt_list[tid]->data_list[idx] += 1;
		}
	}
	for (int i = 0; i < idx_list_ovl->size; i++)
	{
		if (i > 0 && idx_list_ovl->data_list[i] == idx_list_ovl->data_list[i-1]){ continue; }
		idx = idx_list_ovl->data_list[i];
		if (idx >= 0 && idx < wg_ovl_cnt_list[tid]->size){
			wg_ovl_cnt_list[tid]->data_list[idx] += 1;
		}
	}

	return 0;
}

static int cal_twin_win_bcd_cnt(const char * bcd21_file, const char * output_file, const char * faidx_file, int shift_size, int win_size)
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
	int mid_pos;
	int w1_cnt, w2_cnt, ovl_cnt;

	INT_LIST ** wg_win1_cnt_list;
	INT_LIST ** wg_win2_cnt_list;
	INT_LIST ** wg_ovl_cnt_list;
	INT_LIST * read_pos_list;

	INT_LIST * seg_start_list;
	INT_LIST * seg_end_list;
	INT_LIST * idx_list_win1;
	INT_LIST * idx_list_win2;
	INT_LIST * idx_list_ovl;

	INT_LIST * chr_length_list;

	line = (char *) calloc (LINE_MAX, sizeof(char));
	curr_bcd = (char *) calloc (128, sizeof(char));
	new_bcd = (char *) calloc (128, sizeof(char));

	chr_length_list = get_chr_length_from_faidx_file(faidx_file);
	n_chr = chr_length_list->size;    

	wg_win1_cnt_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
	wg_win2_cnt_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
	wg_ovl_cnt_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

	for (int tid = 0; tid < n_chr; tid++)
	{
		chr_len = chr_length_list->data_list[tid]; 
		n_bin = chr_len/shift_size + 2;
		wg_win1_cnt_list[tid] = init_int_list(n_bin); 
		wg_win2_cnt_list[tid] = init_int_list(n_bin); 
		wg_ovl_cnt_list[tid] = init_int_list(n_bin); 
		wg_win1_cnt_list[tid]->size = n_bin;
		wg_win2_cnt_list[tid]->size = n_bin;
		wg_ovl_cnt_list[tid]->size = n_bin;
	}

	seg_start_list = init_int_list(read_pos_list->capacity);
	seg_end_list = init_int_list(read_pos_list->capacity);
	idx_list_win1 = init_int_list(read_pos_list->capacity);
	idx_list_win2 = init_int_list(read_pos_list->capacity);
	idx_list_ovl = init_int_list(read_pos_list->capacity);

	bcd21_fp = gzopen(bcd21_file, "r");
	if (Z_NULL == bcd21_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bcd21_file);
		exit(1);
	}


	read_pos_list = init_int_list( (int) 3e5 );
	curr_tid = -1;
	while (gzgets(bcd21_fp, line, LINE_MAX))
	{
		if (line[0] == '#') { continue; }
		sscanf(line, "%d\t%d\t%d\t%d\t%s\t%*s\n", &new_tid, &start_pos, &end_pos, &mapq, new_bcd);
		avg_pos = (start_pos + end_pos) / 2;
		if ( read_pos_list->size == 0)
		{
			curr_tid = new_tid;
			strcpy(curr_bcd, new_bcd);
			append_int_list(read_pos_list, avg_pos);
			continue;
		}
		if (curr_tid != new_tid || strcmp(curr_bcd, new_bcd) != 0) {
			process_one_barcode_in_one_chr(curr_tid, read_pos_list, wg_win1_cnt_list, wg_win2_cnt_list, wg_ovl_cnt_list, shift_size, win_size, chr_length_list, seg_start_list, seg_end_list, idx_list_win1, idx_list_win2, idx_list_ovl);
			curr_tid = new_tid;
			strcpy(curr_bcd, new_bcd);
			reset_int_list(read_pos_list);
			append_int_list(read_pos_list, avg_pos);
		}else{
			append_int_list(read_pos_list, avg_pos);
		}
	}
	if (read_pos_list->size > 0){
		process_one_barcode_in_one_chr(curr_tid, read_pos_list, wg_win1_cnt_list, wg_win2_cnt_list, wg_ovl_cnt_list, shift_size, win_size, chr_length_list, seg_start_list, seg_end_list, idx_list_win1, idx_list_win2, idx_list_ovl);
	}

	out_fp = fopen(output_file, "w");
	if (NULL == out_fp){
		fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
		exit(1);
	}

	// fprintf(out_fp, "#shift_size=%d\n", shift_size);
	for (int tid = 0; tid < n_chr; tid++)
	{
		n_bin = wg_win1_cnt_list[tid]->size;
		// fprintf(out_fp, ">tid=%d;n_bin=%d\n", tid, n_bin);
		for (int idx = 0; idx < n_bin; idx ++)
		{
			mid_pos = idx * shift_size;
			w1_cnt = wg_win1_cnt_list[tid]->data_list[idx];
			w2_cnt = wg_win2_cnt_list[tid]->data_list[idx];
			ovl_cnt = wg_ovl_cnt_list[tid]->data_list[idx];

			fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\n", tid, mid_pos, mid_pos + win_size, w1_cnt, w2_cnt, ovl_cnt); 
		}
	}

	fclose(out_fp);

	for (int tid = 0; tid < n_chr; tid++)
	{
		free_int_list(wg_win1_cnt_list[tid]);
		free_int_list(wg_win2_cnt_list[tid]);
		free_int_list(wg_ovl_cnt_list[tid]);
	}

	free_int_list(read_pos_list);
	free(wg_win1_cnt_list);
	free(wg_win2_cnt_list);
	free(wg_ovl_cnt_list);

	free_int_list(idx_list_win1);
	free_int_list(idx_list_win2);
	free_int_list(seg_start_list);
	free_int_list(seg_end_list);
	free_int_list(idx_list_ovl);

	free(line);
	free(curr_bcd);
	free(new_bcd);

	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 5){ 
		fprintf (stderr, "Usage: cal_twin_win_bcd_cnt <input.bcd21.gz> <output_file> <faidx_file>\n");
		return 1;
	}

	char * bcd21_file;
	char * output_file;
	char * faidx_file;
	int shift_size;
	int win_size; 
	int ret;

	bcd21_file  = argv[1]; 
	output_file = argv[2];
	faidx_file  = argv[3];

	shift_size = 100;
	win_size = 100 * 100;

	if (shift_size < 10){
		fprintf(stderr, "ERROR! shift_size should be at least 10. Your shift_size value is: %d\n", shift_size);
		exit(1);
	}

	if (win_size % shift_size != 0){
		win_size -= win_size % shift_size;
	}
	ret = cal_twin_win_bcd_cnt(bcd21_file, output_file, faidx_file, shift_size, win_size);

	return ret;
}
