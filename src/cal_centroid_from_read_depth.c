#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "tk.h"

int get_read_depth_from_read_depth_file(const char * input_read_depth_file, INT_LIST ** wg_high_qual_read_depth_list, INT_LIST ** wg_total_read_depth_list, int bin_size, int n_chr)
{
    gzFile input_read_depth_fp;
    char * line;
    int tid;
    int start_pos;
    int end_pos;
	int mapq;
	int start_idx, end_idx;
	int map_len;
    int idx;
    int line_cnt = 0;

	double high_qual_depth, total_depth;
	

    line = (char *) calloc (LINE_MAX, sizeof(char));
    input_read_depth_fp = gzopen(input_read_depth_file, "r");

    if (Z_NULL == input_read_depth_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", input_read_depth_file);
        exit(1);
    }

    while (gzgets(input_read_depth_fp, line, LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        sscanf(line, "%d\t%d\t%d\t%lf\t%lf\t%*s\n", &tid, &start_pos, &end_pos, &high_qual_depth, &total_depth);
		if (tid >= n_chr){
			fprintf(stderr, "ERROR! The faidx file and input depth file don't match! tid is more than number of chr\n");
			exit(1);
		}

		idx = start_pos / bin_size;
		if (idx >= wg_high_qual_read_depth_list[tid]->size)
		{
			fprintf(stderr, "ERROR! The faidx file and input depth file don't match! start_pos is more than chr length\n");
			exit(1);
		}
		wg_high_qual_read_depth_list[tid]->data_list[idx] = (int) high_qual_depth;
		wg_total_read_depth_list[tid]->data_list[idx] = (int) total_depth;

    }

    free(line);
    gzclose(input_read_depth_fp);

    return 0;
}


int get_bin_size_from_input_read_depth_file(const char * input_read_depth_file)
{
    gzFile input_read_depth_fp;
	char * line;
	int tid;
	int start_pos, end_pos;
	int bin_size;

    line = (char *) calloc (LINE_MAX, sizeof(char));
    input_read_depth_fp = gzopen(input_read_depth_file, "r");
    if (Z_NULL == input_read_depth_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", input_read_depth_file);
        exit(1);
    }
	bin_size = -1;
    while (gzgets(input_read_depth_fp, line, LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        sscanf(line, "%*s\t%d\t%d\t%*s\n", &start_pos, &end_pos);
		bin_size = end_pos - start_pos;
		break;
	}
	free(line);
	gzclose(input_read_depth_fp);

	return bin_size;
}
		
double cal_centroid_from_depth_list (INT_LIST ** read_depth_list, int tid, int win_start_pos, int win_end_pos, int bin_size)
{
	double centroid;
	double total_depth;
	int start_idx, end_idx;
	if (win_start_pos < 0 || win_end_pos < 0){
		return -1.0;
	}
	start_idx = win_start_pos / bin_size;
	end_idx = win_end_pos / bin_size;

	if (end_idx >= read_depth_list[tid]->size || start_idx >= read_depth_list[tid]->size){
		return -1.0;
	}
	centroid = 0.0;
	total_depth = 0.0;
	for (int idx = start_idx; idx < end_idx; idx++)
	{
		total_depth += read_depth_list[tid]->data_list[idx]; 
		centroid += (double) idx * (double) bin_size * (double) read_depth_list[tid]->data_list[idx];
	}
	if (total_depth == 0.0){
		return -1;
	}
	centroid = centroid / total_depth;

	return centroid;
}

int cal_centroid_of_read_depth (const char * input_read_depth_file, const char * input_twin_window_file, const char * output_file, const char * faidx_file)
{
	FILE * out_fp;
	gzFile input_twin_window_fp;
	INT_LIST ** wg_high_qual_read_depth_list; // 
    INT_LIST ** wg_total_read_depth_list;
    INT_LIST * chr_length_list;
    int n_chr;
    int n_bin;
	int bin_size;
    int chr_len;
	int bin_start_pos;
	double high_qual_read_depth, total_read_depth;
	char * line;

	int tid, win1_start_pos, win1_end_pos, win2_start_pos, win2_end_pos; 
	int n_bcd_win1, n_bcd_win2, n_ovl_bcd;
	double high_qual_centroid1, high_qual_centroid2, all_reads_centroid1, all_reads_centroid2;

    line = (char *) calloc (LINE_MAX, sizeof(char));
	bin_size = get_bin_size_from_input_read_depth_file(input_read_depth_file);
	fprintf (stderr, "bin_size = %d\n", bin_size);

    chr_length_list = get_chr_length_from_faidx_file(faidx_file);
    n_chr = chr_length_list->size;    
    wg_high_qual_read_depth_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
    wg_total_read_depth_list  = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

    for (int tid = 0; tid < n_chr; tid++)
    {
        chr_len = chr_length_list->data_list[tid]; 
        n_bin = chr_len/bin_size + 2;

        wg_high_qual_read_depth_list[tid] = init_int_list(n_bin); 
        wg_total_read_depth_list[tid] = init_int_list(n_bin); 
		wg_high_qual_read_depth_list[tid]->size = n_bin;
		wg_total_read_depth_list[tid]->size = n_bin;
    }

	fprintf(stderr, "reading read depth file: %s\n", input_read_depth_file);
    get_read_depth_from_read_depth_file(input_read_depth_file, wg_high_qual_read_depth_list, wg_total_read_depth_list, bin_size, n_chr);
	fprintf(stderr, "finished reading read depth file: %s\n", input_read_depth_file);

    out_fp = fopen(output_file, "w");
    if (NULL == out_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
        exit(1);
    }

	fprintf(out_fp, "#tid\twin2_start_pos\twin2_end_pos\tn_bcd_win1\tn_bcd_win2\tn_ovl_bcd\thigh_qual_centroid1\thigh_qual_centroid2\tall_reads_centroid1\tall_reads_centroid2\n");

    input_twin_window_fp = gzopen(input_twin_window_file, "r");
    if (Z_NULL == input_twin_window_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", input_read_depth_file);
        exit(1);
    }
	int win_cnt = 0;
    while (gzgets(input_twin_window_fp, line, LINE_MAX))
	{
        if (line[0] == '#') { continue; }
        sscanf(line, "%d\t%d\t%d\t%d\t%d\t%d\n", &tid, &win2_start_pos, &win2_end_pos, &n_bcd_win1, &n_bcd_win2, &n_ovl_bcd);
		
		win1_start_pos = win2_start_pos - (win2_end_pos - win2_start_pos);	
		win1_end_pos = win2_start_pos;

		high_qual_centroid1 = cal_centroid_from_depth_list (wg_high_qual_read_depth_list, tid, win1_start_pos, win1_end_pos, bin_size);
		high_qual_centroid2 = cal_centroid_from_depth_list (wg_high_qual_read_depth_list, tid, win2_start_pos, win2_end_pos, bin_size);
		all_reads_centroid1 = cal_centroid_from_depth_list (wg_total_read_depth_list, tid, win1_start_pos, win1_end_pos, bin_size);
		all_reads_centroid2 = cal_centroid_from_depth_list (wg_total_read_depth_list, tid, win2_start_pos, win2_end_pos, bin_size);
		fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", tid, win2_start_pos, win2_end_pos, n_bcd_win1, n_bcd_win2, n_ovl_bcd, high_qual_centroid1, high_qual_centroid2, all_reads_centroid1, all_reads_centroid2);
		win_cnt += 1;
		if (win_cnt % 10000000 == 0){
			fprintf(stderr, "processed %d windows\n", win_cnt);
		}

	}

    fclose(out_fp);
    gzclose(input_twin_window_fp);

    for (int tid = 0; tid < n_chr; tid++)
    {
        free_int_list(wg_high_qual_read_depth_list[tid]);
        free_int_list(wg_total_read_depth_list[tid]);
    }
    free(wg_high_qual_read_depth_list);
    free(wg_total_read_depth_list);
	free(line);

    return 0;
}

int usage(FILE * fp)
{
	fprintf (fp, "Usage: cal_centroid_from_read_depth <input_read_depth_file> <input_twin_window_file (bcd11_file)> <output_file (bcd12_file)> <faidx_file>\n");
	return 0;
}

int main(int argc, char * argv[])
{
    if (argc < 5){ 
		usage(stderr);
        return 1;
    }

    char * input_read_depth_file;
    char * input_twin_window_file;
    char * output_file;
    char * faidx_file;
    int ret;

    input_read_depth_file  = argv[1]; 
    input_twin_window_file = argv[2]; 
    output_file = argv[3];
    faidx_file  = argv[4];

    ret = cal_centroid_of_read_depth(input_read_depth_file, input_twin_window_file, output_file, faidx_file);

    return ret;
}
