#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "tk.h"

static int usage(FILE * fp) 
{
    fprintf (fp, "Usage: grid_overlap <input.bcd22> <output_file> <min_num_shared_barcodes> <faidx_file>\n");
    return 0;
}

static int two_chr_grid_overlap(int tid1, int tid2, INT_LIST * chr_length_list, const char * input_bcd22_file, int bin_size, int min_num_shared_barcodes, const char * output_file, int16_t ** n_overlap_bcd)
{
    int tid;
    int start;
    int end;
    int frm_length;
    char bcd[20] = {0};
    char exist_bcd[20] = {0};
    FILE * input_bcd22_fp;
    FILE * out_fp;
    char * line;
    char * temp;
    int tab_cnt;
    int start_idx1, end_idx1, start_idx2, end_idx2;
    INT_LIST * start_list1, *end_list1, * start_list2, *end_list2;
    int max_frm_per_bcd = 100;
    int x, y;
    int n_bin1, n_bin2;
    int k;

    line = (char *) calloc (524288, sizeof(char));
    temp = (char *) calloc (524288, sizeof(char));


    start_list1 = init_int_list (max_frm_per_bcd);
    end_list1 = init_int_list (max_frm_per_bcd);
    start_list2 = init_int_list (max_frm_per_bcd);
    end_list2 = init_int_list (max_frm_per_bcd);

    n_bin1 = chr_length_list->data_list[tid1] / bin_size + 1;
    n_bin2 = chr_length_list->data_list[tid2] / bin_size + 1;
    fprintf(stderr, "n_bin1=%d, n_bin2=%d\n", n_bin1, n_bin2);

    input_bcd22_fp = fopen (input_bcd22_file, "r");
    if (NULL == input_bcd22_fp){
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", input_bcd22_file);
        exit(1);
    }

    out_fp = fopen(output_file, "w");
    if (NULL == out_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
        exit(1);
    }

    while (fgets(line, 524288, input_bcd22_fp)) {
        if (line[0] == '#'){ continue; }
        tab_cnt = 0;
        k = 0;
        while (line[k] != 0 && line[k] != '\n')
        {
            if (line[k] == '\t'){
                tab_cnt += 1;
            }
            if (tab_cnt > 5){
                line[k] = 0;
            }
            k += 1;
        }
        tid = -1;
        start = -1;
        end = -1;
        frm_length = 0;
        bcd[0] = 0;
        sscanf(line, "%d\t%d\t%d\t%d\t%s\t%*s\n", &tid, &start, &end, &frm_length, bcd, temp);    
        if (tid != tid1 && tid != tid2){ continue; }

        if (exist_bcd[0] == 0 ) { strcpy(exist_bcd, bcd); }
        if (strcmp(exist_bcd, bcd) == 0)
        {
            if (tid == tid1){
                append_int_list(start_list1, start);
                append_int_list(end_list1, end);
            }
            if (tid == tid2){
                append_int_list(start_list2, start);
                append_int_list(end_list2, end);
            }
        }else{
            for (int i = 0; i < start_list1->size; i++)
            {
                start_idx1 = (int) ((start_list1->data_list[i]) / bin_size); 
                end_idx1 = (int) ((end_list1->data_list[i]) / bin_size); 
                for (int j = 0; j < start_list2->size; j++) 
                {
                    if (i == j && tid1 == tid2) { continue; }

                    start_idx2 = (int) ((start_list2->data_list[j]) / bin_size); 
                    end_idx2 = (int) ((end_list2->data_list[j]) / bin_size); 

                    if (tid1 == tid2 && ( (end_idx2 - start_idx1) * bin_size < 150000 && (end_idx1 - start_idx1)*bin_size < 150000)) 
                    {
                        // if the two fragments are within 150kb, skip
                        continue; 
                    }

                    for (x = start_idx1; x <= end_idx1; x++)
                    {
                        for (y = start_idx2; y <= end_idx2; y++)
                        {
                            if (n_overlap_bcd[x][y] < INT16_MAX - 2){ 
                                n_overlap_bcd[x][y] += 1; 
                            }
                        }
                    }
                }
            }
            start_list1->size = 0; 
            end_list1->size = 0; 
            start_list2->size = 0;
            end_list2->size = 0;
            for (int i = 0; i < start_list1->size; i++){
                start_list1->data_list[i] = 0;
            }
            for (int i = 0; i < end_list1->size; i++){
                end_list1->data_list[i] = 0;
            }
            for (int i = 0; i < start_list2->size; i++){
                start_list2->data_list[i] = 0;
            }
            for (int i = 0; i < end_list2->size; i++){
                end_list2->data_list[i] = 0;
            }

            strcpy(exist_bcd, bcd);

            if (tid == tid1){
                append_int_list(start_list1, start);
                append_int_list(end_list1, end);
            }
            if (tid == tid2){
                append_int_list(start_list2, start);
                append_int_list(end_list2, end);
            }
        }
    }

    for (int i = 0; i < n_bin1; i++) {
        for (int j = 0; j < n_bin2; j++) {
            if (n_overlap_bcd[i][j] >= min_num_shared_barcodes && ( (i - j)*bin_size > 150000 || (j - i)*bin_size > 150000)) {
                fprintf(out_fp, "%d\t%d\t%d\t%d\t%d\n",tid1, tid2, i*bin_size, j*bin_size, n_overlap_bcd[i][j]);
            }
        }
    }

    // free memory 
    free(line);
    free(start_list1->data_list);
    free(end_list1->data_list);
    free(start_list2->data_list);
    free(end_list2->data_list);

    // close files
    fclose(out_fp);
    fclose(input_bcd22_fp);

    return 0;
}


int whole_genome_grid_overlap(const char * input_bcd22_file, const char * output_file, int min_num_shared_barcodes, const char * faidx_file)
{
    FILE * faidx_fp;
    int bin_size;
    int n_chr;
    INT_LIST * chr_length_list;
    int tid1, tid2;
    int max_chr_length;
    int16_t ** n_overlap_bcd = NULL;
    int max_n_bin;
    
    bin_size = 2000;

    chr_length_list = get_chr_length_from_faidx_file(faidx_file);
    n_chr = chr_length_list->size;
    fprintf(stderr, "number of chr:%d\n", n_chr);

    max_chr_length = 0;
    for (int i = 0; i < n_chr; i++) {
        if (chr_length_list->data_list[i] > max_chr_length) {
            max_chr_length = chr_length_list->data_list[i];
        }
    }
    fprintf(stderr, "max chr length: %d\n", max_chr_length);
    max_n_bin = (int) (max_chr_length / bin_size) + 10;
    fprintf(stderr, "max number of bin: %d\n", max_n_bin);
    n_overlap_bcd = (int16_t **) calloc(max_n_bin, sizeof(int16_t*));
    if (NULL == n_overlap_bcd){
        fprintf(stderr, "ERROR! Failed to alloc memory for n_overlap_bcd!\n");
        exit(1);
    }
    for (int i = 0; i < max_n_bin; i++) {
        n_overlap_bcd[i] = (int16_t *) calloc(max_n_bin, sizeof(int16_t)); 
        if (NULL == n_overlap_bcd[i]){
            fprintf(stderr, "ERROR! Failed to alloc memory for n_overlap_bcd[%d]!\n", i);
            exit(1);
        }
    }
    
    for (tid1 = 0; tid1 < n_chr; tid1++){
        for (tid2 = 0; tid2 < n_chr; tid2++)
        {
            fprintf(stderr, "tid1 = %d, tid2 = %d\n", tid1, tid2);
            fprintf(stderr, "Initializing memory\n", tid1, tid2);
            // if (tid1 != 0 || tid2 != 0) { continue; }
            for (int i = 0; i < max_n_bin; i++)
            {
                for (int j = 0; j < max_n_bin; j++)
                {
                    n_overlap_bcd[i][j] = 0;
                }
            }
            fprintf(stderr, "grid overlap search started\n", tid1, tid2);
            two_chr_grid_overlap(tid1, tid2, chr_length_list, input_bcd22_file, bin_size, min_num_shared_barcodes, output_file, n_overlap_bcd);
            fprintf(stderr, "grid overlap search finished\n", tid1, tid2);
        }
    }

   
    for (int i = 0; i < max_n_bin; i++) {
        free(n_overlap_bcd[i]);
    }
    free(n_overlap_bcd);
    return 0;
}

int main(int argc, char * argv[])
{

    char * input_bcd22_file;
    char * output_file;
    char * faidx_file;
    int min_num_shared_barcodes;

    if (argc < 5){
        usage(stderr);
        return 1;
    }

    input_bcd22_file = argv[1];
    output_file = argv[2];
    min_num_shared_barcodes = atoi(argv[3]);
    faidx_file = argv[4];

    whole_genome_grid_overlap(input_bcd22_file, output_file, min_num_shared_barcodes, faidx_file);

    return 0;
}
