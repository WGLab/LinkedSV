#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "tk.h"

#define MAX_READ_PER_CHROM_PER_BARCODE 300000
#define MIN_MAPQ 20
#define FIX_LENGTH (uint64_t)1e10

typedef struct {
	int32_t tid;
	int32_t idx; // idx = (int) (pos / bin_size)
	uint32_t uint_bcd;
} READ_CORE;

typedef struct {
	READ_CORE * data_list;
	int32_t size;
	int32_t capacity;
} READ_CORE_LIST;

READ_CORE_LIST * init_read_core_list(int capacity)
{
	READ_CORE_LIST * new_read_core_list;
	new_read_core_list = (READ_CORE_LIST *) calloc(1, sizeof(READ_CORE_LIST));
	new_read_core_list->size = 0;
	new_read_core_list->capacity = capacity;
	new_read_core_list->data_list = (READ_CORE *) calloc(new_read_core_list->capacity, sizeof(READ_CORE));

	if (new_read_core_list->data_list == NULL){
		fprintf(stderr, "ERROR! Failed to alloc memory for READ_CORE_LIST (capacity=%d)\n", new_read_core_list->capacity);
		exit(1);
	}   

	return new_read_core_list; 

}

int append_read_core_list (READ_CORE_LIST * read_core_list,  const READ_CORE * read_core)
{
	if (read_core_list->capacity == read_core_list->size)
	{   
		read_core_list->capacity  = (int)(read_core_list->capacity * 2);
		read_core_list->data_list = (READ_CORE *) realloc (read_core_list->data_list, read_core_list->capacity*sizeof(READ_CORE));
		if (read_core_list->data_list == NULL){
			fprintf(stderr, "ERROR! Failed to realloc memory for read_core_list (capacity=%d)\n", read_core_list->capacity);
			exit(1);
		}
	}

	//read_core_list->data_list[read_core_list->size].tid = read_core->tid;
	//read_core_list->data_list[read_core_list->size].pos = read_core->pos;
	//read_core_list->data_list[read_core_list->size].uint_bcd = read_core->uint_bcd;

	memcpy( &(read_core_list->data_list[read_core_list->size]), read_core, sizeof(READ_CORE));

	read_core_list->size += 1;

	return 0;
}

int free_read_core_list(READ_CORE_LIST * read_core_list)
{
	free(read_core_list->data_list);
	free(read_core_list);
	return 0;
}

INT_LIST ** get_target_region_from_bedpe_list(BEDPE_CORE_LIST * in_bedpe_list, INT_LIST * chr_length_list, int bin_size)
{
	INT_LIST ** wg_target_region_list;
	int n_chr = chr_length_list->size;    
	int tid, start_idx, end_idx;
	int chr_len;
	int n_bin;
	int target_length;

	wg_target_region_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

	for (int tid = 0; tid < n_chr; tid++)
	{
		chr_len = chr_length_list->data_list[tid]; 
		n_bin = chr_len/bin_size + 2;
		wg_target_region_list[tid] = init_int_list(n_bin); 
		wg_target_region_list[tid]->size = n_bin;
	}

	for (int i = 0; i < in_bedpe_list->size; i++)
	{
		tid = in_bedpe_list->data_list[i].tid1;	
		start_idx = in_bedpe_list->data_list[i].start1 / bin_size;
		end_idx = in_bedpe_list->data_list[i].end1 / bin_size + 1;	
		if (tid < 0 || tid >= n_chr || start_idx < 0 || end_idx > wg_target_region_list[tid]->size) {
			fprintf(stderr, "WARNING! Input bedpe_file and faidx_file don't match! current tid | start | end: %d | %d | %d\n", tid, in_bedpe_list->data_list[i].start1, in_bedpe_list->data_list[i].end1);
			if (tid < 0 || tid > n_chr ) { continue; } 
			if (start_idx < 0) { start_idx = 0; }
			if (end_idx > wg_target_region_list[tid]->size ) { end_idx = wg_target_region_list[tid]->size; }
		}
		for (int idx = start_idx; idx < end_idx; idx ++) {
			wg_target_region_list[tid]->data_list[idx] += 1;
		}

		tid = in_bedpe_list->data_list[i].tid2;	
		start_idx = in_bedpe_list->data_list[i].start2 / bin_size;
		end_idx = in_bedpe_list->data_list[i].end2 / bin_size + 1;	
		if (tid < 0 || tid > n_chr || start_idx < 0 || end_idx > wg_target_region_list[tid]->size) {
			fprintf(stderr, "WARNING! Input bedpe_file and faidx_file don't match! current tid | start | end: %d | %d | %d\n", tid, in_bedpe_list->data_list[i].start1, in_bedpe_list->data_list[i].end1);
			if (tid < 0 || tid > n_chr ) { continue; } 
			if (start_idx < 0) { start_idx = 0; }
			if (end_idx > wg_target_region_list[tid]->size ) { end_idx = wg_target_region_list[tid]->size; }
		}
		for (int idx = start_idx; idx < end_idx; idx ++) {
			wg_target_region_list[tid]->data_list[idx] += 1;
		}
	}

	target_length = 0;
	for (int tid = 0; tid < n_chr; tid ++)
	{
		for (int idx = 0; idx < wg_target_region_list[tid]->size; idx++)
		{
			if (wg_target_region_list[tid]->data_list[idx] > 0)
			{
				target_length += 1;
			}
		}
	}

	target_length = target_length * bin_size;
	fprintf(stderr, "target region length is: %d bp\n", target_length);
	return wg_target_region_list;
}

int print_wg_readcore_list(READ_CORE_LIST ** wg_readcore_list, int n_chr)
{
	READ_CORE * readcore;
	FILE * temp_fp;
	temp_fp = open_file("wg_readcore_list.txt", "w");

	for (int tid = 0; tid < n_chr; tid++)
	{
		for (int i = 0; i < wg_readcore_list[tid]->size; i++)
		{
			readcore = wg_readcore_list[tid]->data_list + i;
			fprintf(temp_fp, "%d\t%d\t%u\n", readcore->tid, readcore->idx, readcore->uint_bcd);
		}
	}

	fclose(temp_fp);

	return 0;
}

READ_CORE_LIST ** get_reads_from_bcd21_file(const char * bcd21_file, INT_LIST ** wg_target_region_list, INT_LIST * chr_length_list, int bin_size, int show_low_qual_mapping)
{
	READ_CORE_LIST ** wg_readcore_list;
	gzFile bcd21_fp;
	char * line;
	char * bcd;
	int mapq, flag, hap_type;
	int32_t tid, start_pos, end_pos, avg_pos;
	uint32_t uint_bcd; 
	READ_CORE * readcore;
	READ_CORE * last_readcore;

	int n_chr;
	int idx;
	int n_reads;
	int i;
	int init_list_size = (int)1e5;

	line = (char *) calloc (LINE_MAX, sizeof(char));
	bcd = (char *) calloc (64, sizeof(char));
	readcore = (READ_CORE *)calloc(1, sizeof(READ_CORE));


	n_chr = chr_length_list->size;
	wg_readcore_list = (READ_CORE_LIST **) calloc(n_chr, sizeof(READ_CORE_LIST *));

	for (int tid = 0; tid < n_chr; tid++) {
		wg_readcore_list[tid] = init_read_core_list(init_list_size); 
	}
	
	bcd21_fp = gzopen(bcd21_file, "r");
	if (Z_NULL == bcd21_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bcd21_file);
		exit(1);
	}

	n_reads = 0;
	i = 0;
	while (gzgets(bcd21_fp, line, LINE_MAX))
	{
		if (line[0] == '#') { continue; }
		i++;
		if (i % 100000000 == 0){
			fprintf(stderr, "processed %d million lines in bcd21_file.\n", i/1000000, n_reads);
		}

		sscanf(line, "%d\t%d\t%d\t%d\t%s\t%d\t%*s\t%d\t%*s\n", &tid, &start_pos, &end_pos, &mapq, bcd, &hap_type, &flag);

		if (!show_low_qual_mapping && flag & (256 + 1024 + 2048) ){ continue; }
		if (!show_low_qual_mapping && mapq < MIN_MAPQ) { continue; } 

		avg_pos = (start_pos + end_pos) / 2;
		idx = avg_pos / bin_size;

		if (wg_target_region_list[tid]->data_list[idx] == 0) { continue; } // not in target regions

		if (tid >= n_chr){
			fprintf(stderr, "ERROR! bcd21_file and faidx_file don't match! current tid | pos: %d | %d\n", tid, end_pos);
			exit(1);
		}

		if (idx >= wg_target_region_list[tid]->size){
			fprintf(stderr, "ERROR! bcd21_file and faidx_file don't match! current tid | pos: %d | %d\n", tid, end_pos);
			exit(1);
		}

		uint_bcd = bcdseq2bcdint(bcd);
		readcore->tid = tid;
		readcore->idx = idx;
		readcore->uint_bcd = uint_bcd;

		last_readcore = wg_readcore_list[tid]->data_list + wg_readcore_list[tid]->size - 1; 
		if (idx == last_readcore->idx && uint_bcd == last_readcore->uint_bcd){ continue; }

		append_read_core_list (wg_readcore_list[tid], readcore);
		n_reads++;
	}

	free(line);
	free(bcd);
	free(readcore);
	gzclose(bcd21_fp);

	fprintf(stderr, "Finished reading input file: %s\n", bcd21_file);
	// print_wg_readcore_list(wg_readcore_list, n_chr);

	return wg_readcore_list;
}

READ_CORE_LIST * get_readcore_list_in_region(READ_CORE_LIST ** wg_readcore_list, int n_chr, int tid, int start_idx, int end_idx)
{
	READ_CORE_LIST * region_readcore_list;
	READ_CORE * readcore;

	if (tid >= n_chr) {
		fprintf(stderr, "ERROR! bcd21_file and faidx_file don't match! tid=%d, n_chr=%d\n", tid, n_chr);
		exit(1);
	}

	region_readcore_list = init_read_core_list(1000);

	for (int i = 0; i < wg_readcore_list[tid]->size; i++)
	{
		readcore = wg_readcore_list[tid]->data_list + i;
		if (readcore->idx >= start_idx && readcore->idx < end_idx){
			append_read_core_list (region_readcore_list, readcore);
		}
	}

	return region_readcore_list;

}

READ_CORE_LIST * concatenate_two_readcore_list(READ_CORE_LIST * readcore_list1, READ_CORE_LIST * readcore_list2)
{
	for (int i = 0; i < readcore_list2->size; i++) {
		append_read_core_list(readcore_list1, readcore_list2->data_list + i);
	}

	free_read_core_list(readcore_list2);

	return readcore_list1;
}

static inline int cmpfunc_readcore_bcd_then_pos(const void * r1, const void * r2)
{
	READ_CORE * a = (READ_CORE *) r1;
	READ_CORE * b = (READ_CORE *) r2;

	if (a->uint_bcd == b->uint_bcd){
		if (a->tid == b->tid){
			if (a->idx > b->idx) {
				return 1;
			}else if (a->idx < b->idx){
				return -1;
			}else{
				return 0;
			}
		}else if (a->tid > b->tid){
			return 1;
		}else{
			return -1;
		}
	}else if (a->uint_bcd > b->uint_bcd){
		return 1;
	}else {
		return -1;
	}
}

READ_CORE_LIST * sort_readcore_list_by_barcode_and_then_idx(READ_CORE_LIST * readcore_list)
{
	qsort(readcore_list->data_list, readcore_list->size, sizeof(READ_CORE), cmpfunc_readcore_bcd_then_pos);
	return readcore_list;
}

int process_one_barcode(INT_2D_LIST * num_ovl_bcd, READ_CORE_LIST * one_barcode_readcore_list, const int tid1, const int start_idx1, const int end_idx1, const int tid2, const int start_idx2, const int end_idx2)
{
	READ_CORE * readcore;
	READ_CORE * readcore1;
	READ_CORE * readcore2;

	int start_i1, end_i1, start_i2, end_i2;
	uint64_t key;
	uint64_t start_key1, start_key2, end_key1, end_key2;
	int x, y;

	start_i1 = -10;
	end_i1   = -10;
	start_i2 = -10; // self included
	end_i2   = -10; // self included

	if (tid1 < 0 || tid2 < 0 || start_idx1 < 0 || start_idx2 < 0 || end_idx1 < 0 || end_idx2 < 0) { return 0; }

	//printf("processing one barcode:%u, number of reads:%d\n", one_barcode_readcore_list->data_list[0].uint_bcd, one_barcode_readcore_list->size);
	//printf("debug: tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2=%d,%d,%d,%d,%d,%d\n", tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2);

	sort_readcore_list_by_barcode_and_then_idx(one_barcode_readcore_list);

	/*printf("debug: reads are:");
	for (int i = 0; i < one_barcode_readcore_list->size; i++) {
		printf("%llu, ", (uint64_t) one_barcode_readcore_list->data_list[i].tid * FIX_LENGTH + one_barcode_readcore_list->data_list[i].idx);
	}
	printf("\n");

	*/
	start_key1 = (uint64_t) tid1 * FIX_LENGTH + start_idx1;
	start_key2 = (uint64_t) tid2 * FIX_LENGTH + start_idx2;

	end_key1   = (uint64_t) tid1 * FIX_LENGTH + end_idx1;
	end_key2   = (uint64_t) tid2 * FIX_LENGTH + end_idx2;

	// printf("debug: start_key1, end_key1, start_key2, end_key2=%llu,%llu,%llu,%llu\n", start_key1, end_key1, start_key2, end_key2); 

	for (int i = 0; i < one_barcode_readcore_list->size; i++)
	{
		readcore = one_barcode_readcore_list->data_list + i;
		key = ((uint64_t) readcore->tid) * FIX_LENGTH + readcore->idx;
		// printf("key=%llu\n", key);

		if (start_i1 == -10 && key >= start_key1 && key < end_key1 ) { start_i1 = i; }
		if (start_i2 == -10 && key >= start_key2 && key < end_key2 ) { start_i2 = i; }

		if (key >= start_key1 && key < end_key1 ) { end_i1 = i; }
		if (key >= start_key2 && key < end_key2 ) { end_i2 = i; }
	}

	// printf("debug: start_i1=%d,end_i1=%d,start_i2=%d,end_i2=%d\n", start_i1, end_i1, start_i2, end_i2);

	if (start_i1 == -10 || start_i2 == -10 || end_i1 == -10 || end_i2 == -10 ) { return 0; }


	for (int i = start_i1; i <= end_i1; i++)
	{
		for (int j = start_i2; j <= end_i2; j++)
		{
			readcore1 = one_barcode_readcore_list->data_list + i;
			readcore2 = one_barcode_readcore_list->data_list + j;
			x = readcore1->idx - start_idx1;
			y = readcore2->idx - start_idx2;
			num_ovl_bcd->data[x][y] += 1;
			// printf("x=%d,y=%d\n", x, y);
		}
	}

	return 0;
}

int max_int(int a, int b)
{
	if (a > b){
		return a;
	}else{
		return b;
	}
}

int min_int(int a, int b)
{
	if (a < b){
		return a;
	}else{
		return b;
	}
}

int region_has_overlap (int tid1, int start_idx1, int end_idx1, int tid2, int start_idx2, int end_idx2)
{
	if (tid1 != tid2){ return 0; }
	if (min_int(end_idx1, end_idx2) >= max_int(start_idx1, start_idx2)){
		return 1;
	}else{
		return 0;
	}
}

INT_2D_LIST * cal_overlapping_barcodes_for1bedpe(BEDPE_CORE * bedpe, READ_CORE_LIST ** wg_readcore_list, int bin_size, CHR_INFO * chr_info)
{

	int n_bin_x;
	int n_bin_y;
	int start_idx1, end_idx1;
	int start_idx2, end_idx2;
	int tid, tid1, tid2;
	FILE * bcd21_fp;
	int n_chr;
	uint32_t curr_bcd;
	int min_start_idx, max_end_idx;

	INT_2D_LIST * num_ovl_bcd;

	READ_CORE_LIST * region1_readcore_list;
	READ_CORE_LIST * region2_readcore_list;
	READ_CORE_LIST * merged_readcore_list;
	READ_CORE_LIST * one_barcode_readcore_list; 

	n_chr = chr_info->chr_length_list->size;    
	tid1 = bedpe->tid1;
	start_idx1 = bedpe->start1 / bin_size;
	end_idx1 = bedpe->end1 / bin_size + 1;

	tid2 = bedpe->tid2;
	start_idx2 = bedpe->start2 / bin_size;
	end_idx2 = bedpe->end2 / bin_size + 1;
	

	if (region_has_overlap (tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2))
	{
		min_start_idx = min_int(start_idx1, end_idx1);
		max_end_idx = max_int(end_idx1, end_idx2);
		// printf("debug: %d:%d-%d; %d:%d-%d, merged region: %d-%d\n", tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2, min_start_idx, max_end_idx);
		merged_readcore_list = get_readcore_list_in_region(wg_readcore_list, n_chr, tid1, min_start_idx, max_end_idx); 
	}else{
		region1_readcore_list = get_readcore_list_in_region(wg_readcore_list, n_chr, tid1, start_idx1, end_idx1); 
		region2_readcore_list = get_readcore_list_in_region(wg_readcore_list, n_chr, tid2, start_idx2, end_idx2); 
		// printf("debug: len(region1_readcore_list)=%d, len(region2_readcore_list)=%d\n", region1_readcore_list->size, region2_readcore_list->size);
		merged_readcore_list = concatenate_two_readcore_list(region1_readcore_list, region2_readcore_list);  // region1_readcore_list and region2_readcore_list is destroyed 
		// printf("debug: len(merged_readcore_list)=%d", merged_readcore_list->size); 
	}

	sort_readcore_list_by_barcode_and_then_idx(merged_readcore_list);

	n_bin_x = end_idx1 - start_idx1; 
	n_bin_y = end_idx2 - start_idx2;

	num_ovl_bcd = init_int_2d_list(n_bin_x, n_bin_y);
	num_ovl_bcd->size1 = n_bin_x;
	num_ovl_bcd->size2 = n_bin_y;

	if (merged_readcore_list->size == 0) {
		return num_ovl_bcd;
	}

	one_barcode_readcore_list = init_read_core_list(1000);
	for (int i = 1; i < merged_readcore_list->size; i++)
	{
		if (merged_readcore_list->data_list[i].uint_bcd == merged_readcore_list->data_list[i-1].uint_bcd) {
			append_read_core_list(one_barcode_readcore_list, merged_readcore_list->data_list + i);
		}else{
			process_one_barcode(num_ovl_bcd, one_barcode_readcore_list, tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2);
			one_barcode_readcore_list->size = 0;
		}
	}

	process_one_barcode(num_ovl_bcd, one_barcode_readcore_list, tid1, start_idx1, end_idx1, tid2, start_idx2, end_idx2);

	free_read_core_list(one_barcode_readcore_list);
	free_read_core_list(merged_readcore_list);

	return num_ovl_bcd; 

}


static int cal_2d_overlapping_barcodes (const char * bcd21_file, const char * in_bedpe_file, const char * output_file, const char * faidx_file, int bin_size, int show_low_qual_mapping)
{
	FILE * out_fp;
	int n_chr;
	int n_bin;

	int start_pos;
	int end_pos;
	int mapq;
	int avg_pos;
	int tid;

	int flag;
	CHR_INFO * chr_info;

	BEDPE_CORE_LIST * in_bedpe_list;
	BEDPE_CORE * bedpe;

	READ_CORE_LIST ** wg_readcore_list;
	READ_CORE * readcore;
	INT_2D_LIST * num_ovl_bcd;
	int start_idx1, end_idx1, start_idx2, end_idx2;

	in_bedpe_list = read_bedpe_core_file(in_bedpe_file, faidx_file);

	INT_LIST ** wg_target_region_list;

	chr_info = get_chr_info(faidx_file);

	n_chr = chr_info->chr_length_list->size;    

	wg_target_region_list = get_target_region_from_bedpe_list(in_bedpe_list, chr_info->chr_length_list, bin_size);

	wg_readcore_list = get_reads_from_bcd21_file(bcd21_file, wg_target_region_list, chr_info->chr_length_list, bin_size, show_low_qual_mapping);

	out_fp = open_file(output_file, "w");
	for (int i = 0; i < in_bedpe_list->size; i++)
	{
		bedpe = in_bedpe_list->data_list + i;
		fprintf(stderr, "calculating overlapping barcodes for one bedpe region: %s\t%d\t%d\t", chr_info->chrname_list->data_list[bedpe->tid1], bedpe->start1, bedpe->end1);
		num_ovl_bcd = cal_overlapping_barcodes_for1bedpe(bedpe, wg_readcore_list, bin_size, chr_info);
		start_idx1 = bedpe->start1 / bin_size;
		end_idx1 = bedpe->end1 / bin_size + 1;

		start_idx2 = bedpe->start2 / bin_size;
		end_idx2 = bedpe->end2 / bin_size + 1;
	
		fprintf(out_fp, "##%s:%d-%d; ", chr_info->chrname_list->data_list[bedpe->tid1], bedpe->start1, bedpe->end1);
		fprintf(out_fp, "%s:%d-%d\n",  chr_info->chrname_list->data_list[bedpe->tid2], bedpe->start2, bedpe->end2);

		fprintf(out_fp, "#xmin,xmax,ymin,ymax=%d,%d,%d,%d\n", start_idx1 * bin_size , end_idx1 * bin_size, start_idx2 * bin_size, end_idx2 * bin_size);

		fprintf(stderr, "%s\t%d\t%d\n",  chr_info->chrname_list->data_list[bedpe->tid2], bedpe->start2, bedpe->end2);
		for (int j = 0; j < num_ovl_bcd->size1; j++) {
			for (int k = 0; k < num_ovl_bcd->size2-1; k++) {
				fprintf(out_fp, "%d\t", num_ovl_bcd->data[j][k]);
			}
			int k = num_ovl_bcd->size2 - 1;
			fprintf(out_fp, "%d\n", num_ovl_bcd->data[j][k]);
		}
	}

	fclose(out_fp);

	for (int tid = 0; tid < n_chr; tid++)
	{
		free_int_list(wg_target_region_list[tid]);
	}

	free(wg_target_region_list);
	free_int_2d_list(num_ovl_bcd);

	return 0;
}

int usage(FILE * fp)
{
	fprintf (fp, "Usage: cal_2d_overlapping_barcodes  <input.bcd21.gz> <input.bedpe> <output_file> <faidx_file> <bin_size> <show_low_qual_mapping (1=True, 0=False)>\n");
	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 7){ 
		usage(stderr);
		return 1;
	}

	char * bcd21_file;
	char * in_bedpe_file;
	char * output_file;
	char * faidx_file;
	int bin_size;
	int show_low_qual_mapping;
	int ret;

	bcd21_file = argv[1]; 
	in_bedpe_file = argv[2];
	output_file = argv[3];
	faidx_file = argv[4];
	bin_size = atoi(argv[5]);
	show_low_qual_mapping = atoi(argv[6]);

	if (bin_size < 100){
		fprintf(stderr, "ERROR! bin_size should be at least 100. Your bin_size value is: %d\n", bin_size);
		usage(stderr);
		exit(1);
	}   

	ret = cal_2d_overlapping_barcodes (bcd21_file, in_bedpe_file, output_file, faidx_file, bin_size, show_low_qual_mapping);

	return ret;
}
