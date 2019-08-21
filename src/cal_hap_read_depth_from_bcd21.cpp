#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "tk.h"

int get_read_depth_from_read_bcd21_file(const char * bcd21_file, INT_LIST ** hap0_wg_high_qual_read_depth_list, INT_LIST ** hap0_wg_low_qual_read_depth_list, INT_LIST ** hap1_wg_high_qual_read_depth_list, INT_LIST ** hap1_wg_low_qual_read_depth_list, INT_LIST ** hap2_wg_high_qual_read_depth_list, INT_LIST ** hap2_wg_low_qual_read_depth_list, int bin_size, int n_chr, int mapq_cutoff)

{
    gzFile bcd21_fp;
    char * line;
    int tid;
    int start_pos;
    int end_pos;
	int mapq;
	int start_idx, end_idx;
	// int map_len;
    int idx;
    int line_cnt = 0;

	int flag;
	int hap_type;
	
    INT_LIST ** read_depth_list; 


    line = (char *) calloc (LINE_MAX, sizeof(char));
    bcd21_fp = gzopen(bcd21_file, "r");

    if (Z_NULL == bcd21_fp) {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bcd21_file);
        exit(1);
    }

    while (gzgets(bcd21_fp, line, LINE_MAX))
    {
        if (line[0] == '#') { continue; }
        sscanf(line, "%d\t%d\t%d\t%d\t%*s\t%d\t%*s\t%d\t%*s\n", &tid, &start_pos, &end_pos, &mapq, &hap_type, &flag);
		if (tid >= n_chr){
			fprintf(stderr, "ERROR! The faidx file and input bcd21 file don't match!\n");
			exit(1);
		}

		if (flag & (256 + 1024 + 2048) ){ continue; }


        

        if (mapq >= mapq_cutoff){
            if (hap_type == 1){
                read_depth_list = hap1_wg_high_qual_read_depth_list;
            }else if (hap_type == 2){
                read_depth_list = hap2_wg_high_qual_read_depth_list;
            }else{
                read_depth_list = hap0_wg_high_qual_read_depth_list;
            }
        }else{
            if (hap_type == 1){
                read_depth_list = hap1_wg_low_qual_read_depth_list;
            }else if (hap_type == 2){
                read_depth_list = hap2_wg_low_qual_read_depth_list;
            }else{
                read_depth_list = hap0_wg_low_qual_read_depth_list;
            }
        }

        /*
        idx = (start_pos + end_pos) / ( 2 * bin_size);
        read_depth_list[tid]->data_list[idx] += end_pos - start_pos;
        */

        start_idx = start_pos / bin_size;
        end_idx = end_pos / bin_size; 
        if (end_idx > start_idx){
            read_depth_list[tid]->data_list[start_idx] += bin_size - start_pos % bin_size;
            read_depth_list[tid]->data_list[end_idx] += end_pos % bin_size;
            for (int i = start_idx + 1; i < end_idx; i++)
            {
                read_depth_list[tid]->data_list[i] += bin_size;
            }
        }else{
            read_depth_list[tid]->data_list[start_idx] += end_pos - start_pos; 
        }

        
		line_cnt += 1;
		if (line_cnt % 100000000 == 0){
			fprintf(stderr, "processed %d million reads from bcd21 file: %s\n", line_cnt/1000000, bcd21_file);
        }
    }

    free(line);
    gzclose(bcd21_fp);

    return 0;
}

int cal_bin_read_depth_from_bcd21(const char * bcd21_file, const char * output_file, const char * faidx_file, int bin_size, int mapq_cutoff)
{
    FILE * out_fp;
    INT_LIST ** hap1_wg_high_qual_read_depth_list;
    INT_LIST ** hap2_wg_high_qual_read_depth_list; 
    INT_LIST ** hap0_wg_high_qual_read_depth_list; 


    INT_LIST ** hap1_wg_low_qual_read_depth_list;
    INT_LIST ** hap2_wg_low_qual_read_depth_list;
    INT_LIST ** hap0_wg_low_qual_read_depth_list;

    INT_LIST * chr_length_list;
    int n_chr;
    int n_bin;
    int chr_len;
	int bin_start_pos;



    chr_length_list = get_chr_length_from_faidx_file(faidx_file);
    n_chr = chr_length_list->size;

    hap1_wg_high_qual_read_depth_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
    hap1_wg_low_qual_read_depth_list  = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

    hap2_wg_high_qual_read_depth_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
    hap2_wg_low_qual_read_depth_list  = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));

    hap0_wg_high_qual_read_depth_list = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));
    hap0_wg_low_qual_read_depth_list  = (INT_LIST **) calloc(n_chr, sizeof(INT_LIST *));


    for (int tid = 0; tid < n_chr; tid++)
    {
        chr_len = chr_length_list->data_list[tid]; 
        n_bin = chr_len/bin_size + 2;

        hap1_wg_high_qual_read_depth_list[tid] = init_int_list(n_bin); 
        hap1_wg_low_qual_read_depth_list[tid] = init_int_list(n_bin); 
		hap1_wg_high_qual_read_depth_list[tid]->size = n_bin;
		hap1_wg_low_qual_read_depth_list[tid]->size = n_bin;

        hap2_wg_high_qual_read_depth_list[tid] = init_int_list(n_bin); 
        hap2_wg_low_qual_read_depth_list[tid] = init_int_list(n_bin); 
        hap2_wg_high_qual_read_depth_list[tid]->size = n_bin;
        hap2_wg_low_qual_read_depth_list[tid]->size = n_bin;

        hap0_wg_high_qual_read_depth_list[tid] = init_int_list(n_bin); 
        hap0_wg_low_qual_read_depth_list[tid] = init_int_list(n_bin); 
        hap0_wg_high_qual_read_depth_list[tid]->size = n_bin;
        hap0_wg_low_qual_read_depth_list[tid]->size = n_bin;

    }

    get_read_depth_from_read_bcd21_file (bcd21_file, hap0_wg_high_qual_read_depth_list, hap0_wg_low_qual_read_depth_list, hap1_wg_high_qual_read_depth_list, hap1_wg_low_qual_read_depth_list, hap2_wg_high_qual_read_depth_list, hap2_wg_low_qual_read_depth_list, bin_size, n_chr, mapq_cutoff);

    out_fp = fopen(output_file, "w");
    if (NULL == out_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
        exit(1);
    }
    
	fprintf(out_fp, "#tid\tstart_pos\tend_pos\tall_hap_high_qual_read_depth\tall_hap_total_read_depth\thap1_high_qual_read_depth\thap1_total_read_depth\thap2_high_qual_read_depth\thap2_total_read_depth\thap0_high_qual_read_depth\thap0_total_read_depth\n");

    double all_hap_high_qual_read_depth, all_hap_total_read_depth;
    double hap1_high_qual_read_depth, hap1_total_read_depth;
    double hap2_high_qual_read_depth, hap2_total_read_depth;
    double hap0_high_qual_read_depth, hap0_total_read_depth;

    for (int tid = 0; tid < n_chr; tid++)
    {
		n_bin = hap1_wg_high_qual_read_depth_list[tid]->size; 
        for (int idx = 0; idx < n_bin; idx ++)
		{
			bin_start_pos = idx * bin_size;

            hap1_high_qual_read_depth = (double) hap1_wg_high_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;
            hap2_high_qual_read_depth = (double) hap2_wg_high_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;
            hap0_high_qual_read_depth = (double) hap0_wg_high_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;

            hap1_total_read_depth = hap1_high_qual_read_depth + (double) hap1_wg_low_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;
            hap2_total_read_depth = hap2_high_qual_read_depth + (double) hap2_wg_low_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;
            hap0_total_read_depth = hap0_high_qual_read_depth + (double) hap0_wg_low_qual_read_depth_list[tid]->data_list[idx] / (double) bin_size;

            all_hap_high_qual_read_depth = hap1_high_qual_read_depth + hap2_high_qual_read_depth + hap0_high_qual_read_depth;
            all_hap_total_read_depth = hap1_total_read_depth + hap2_total_read_depth + hap0_total_read_depth;

            fprintf(out_fp, "%d\t%d\t%d\t", tid, bin_start_pos, bin_start_pos + bin_size);
            fprintf(out_fp, "%.2f\t%.2f\t", all_hap_high_qual_read_depth, all_hap_total_read_depth);
            fprintf(out_fp, "%.2f\t%.2f\t", hap1_high_qual_read_depth, hap1_total_read_depth);
            fprintf(out_fp, "%.2f\t%.2f\t", hap2_high_qual_read_depth, hap2_total_read_depth);
            fprintf(out_fp, "%.2f\t%.2f\n", hap0_high_qual_read_depth, hap0_total_read_depth);
		}
    }

    fclose(out_fp);

    for (int tid = 0; tid < n_chr; tid++)
    {
        free_int_list(hap1_wg_high_qual_read_depth_list[tid]);
        free_int_list(hap1_wg_low_qual_read_depth_list[tid]);

        free_int_list(hap2_wg_high_qual_read_depth_list[tid]);
        free_int_list(hap2_wg_low_qual_read_depth_list[tid]);

        free_int_list(hap0_wg_high_qual_read_depth_list[tid]);
        free_int_list(hap0_wg_low_qual_read_depth_list[tid]);
    }

    free(hap1_wg_high_qual_read_depth_list);
    free(hap2_wg_high_qual_read_depth_list); 
    free(hap0_wg_high_qual_read_depth_list); 

    free(hap1_wg_low_qual_read_depth_list);
    free(hap2_wg_low_qual_read_depth_list);
    free(hap0_wg_low_qual_read_depth_list);

    return 0;
}

int main(int argc, char * argv[])
{
    if (argc < 6){ 
        fprintf (stderr, "Usage: cal_bin_read_depth_from_bcd21 <input.bcd21.gz> <output_file> <faidx_file> <bin_size> <mapq_cutoff>\n");
        return 1;
    }

    char * bcd21_file;
    char * output_file;
    char * faidx_file;
	int mapq_cutoff;
    int bin_size;
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

    ret = cal_bin_read_depth_from_bcd21(bcd21_file, output_file, faidx_file, bin_size, mapq_cutoff);

    return ret;
}
