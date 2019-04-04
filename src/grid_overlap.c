#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "tk.h"

typedef struct {
    int16_t tid;
    int32_t start;
    int32_t end;
    char bcd[20];
} FRM_CORE;

typedef struct {
    FRM_CORE * data_list;
    int32_t size;
    int32_t capacity;
} FRM_CORE_LIST;

static int usage(FILE * fp) 
{
    fprintf (fp, "Usage: grid_overlap <input.bcd22> <output_file> <min_num_shared_barcodes> <faidx_file>\n");
    return 0;
}

static int two_chr_grid_overlap(int tid1, int tid2, INT_LIST * chr_length_list, FRM_CORE_LIST * all_frm_list, int bin_size, uint8_t min_num_shared_barcodes, const char * output_file, uint8_t ** n_overlap_bcd)
{
    int tid;
    int start;
    int end;
    int frm_length;
    char * bcd;
    char exist_bcd[20] = {0};
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


    out_fp = fopen(output_file, "a");
    if (NULL == out_fp){
        fprintf(stderr, "ERROR! Failed to open file for writing: %s\n", output_file);
        exit(1);
    }

    for (int frm_id = 0; frm_id < all_frm_list->size; frm_id ++)
    {
        tid = all_frm_list->data_list[frm_id].tid;
        start = all_frm_list->data_list[frm_id].start;
        end = all_frm_list->data_list[frm_id].end;
        bcd = all_frm_list->data_list[frm_id].bcd;

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
                end_idx1   = (int) ((end_list1->data_list[i]) / bin_size); 
                for (int j = 0; j < start_list2->size; j++) 
                {
                    if (i == j && tid1 == tid2) { continue; }

                    start_idx2 = (int) ((start_list2->data_list[j]) / bin_size); 
                    end_idx2   = (int) ((end_list2->data_list[j]) / bin_size); 

                    if (tid1 == tid2 && ( (end_idx2 - start_idx1) * bin_size < 150000 && (end_idx1 - start_idx1)*bin_size < 150000)) 
                    {
                        // if the two fragments are within 150kb, skip
                        continue; 
                    }

                    for (x = start_idx1; x <= end_idx1; x++)
                    {
                        for (y = start_idx2; y <= end_idx2; y++)
                        {
                            if (n_overlap_bcd[x][y] < UINT8_MAX - 1){ 
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
            if (n_overlap_bcd[i][j] >= min_num_shared_barcodes && ( tid1 != tid2 || (i - j)*bin_size > 150000 || (j - i)*bin_size > 150000)) {
                fprintf(out_fp, "%d\t%d\t%d\t%d\t%u\n",tid1, i*bin_size, tid2, j*bin_size, n_overlap_bcd[i][j]);
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

    return 0;
}

FRM_CORE_LIST * read_frm_list_from_bcd22_file(const char * input_bcd22_file)
{
    FRM_CORE_LIST * all_frm_list;
    char * line, *temp;
    int tab_cnt;
    int k;
    char bcd[20] = {0};
    int tid, start, end, frm_length;
    FILE * input_bcd22_fp;

    fprintf(stderr, "start reading input bcd22 file: %s\n", input_bcd22_file);
    line = (char *) calloc (524288, sizeof(char));
    temp = (char *) calloc (524288, sizeof(char));


    input_bcd22_fp = fopen (input_bcd22_file, "r");
    if (NULL == input_bcd22_fp){
        fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", input_bcd22_file);
        exit(1);
    }

    int line_cnt = 0;
    while (fgets(line, 524288, input_bcd22_fp)) {
        if (line[0] == '#'){ continue; }
        line_cnt += 1;
    }

    fprintf(stderr, "total number of fragments in the input bcd22 file is: %d\n", line_cnt);

    all_frm_list = (FRM_CORE_LIST * ) calloc(1, sizeof(FRM_CORE_LIST));
    all_frm_list->capacity = line_cnt + 10;
    all_frm_list->data_list = (FRM_CORE *) calloc (all_frm_list->capacity, sizeof(FRM_CORE));
    all_frm_list->size = 0;

    fseek (input_bcd22_fp, 0, SEEK_SET);  
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

        if (tid >= 0 && start >= 0 && end >= 0){
            all_frm_list->data_list[all_frm_list->size].tid = tid;
            all_frm_list->data_list[all_frm_list->size].start = start;
            all_frm_list->data_list[all_frm_list->size].end = end;
            strcpy(all_frm_list->data_list[all_frm_list->size].bcd, bcd);
            all_frm_list->size += 1;
        }
    }

    free(line);
    free(temp);
    fprintf(stderr, "finished reading input bcd22 file\n");
    return all_frm_list;
}

int whole_genome_grid_overlap(const char * input_bcd22_file, const char * output_file, uint8_t min_num_shared_barcodes, const char * faidx_file)
{
    FILE * faidx_fp;
    int bin_size;
    int n_chr;
    INT_LIST * chr_length_list;
    int tid1, tid2;
    int max_chr_length;
    uint8_t ** n_overlap_bcd = NULL;
    int max_n_bin;

    FRM_CORE_LIST * all_frm_list;
    
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

    // alloc memory for n_overlap_bcd
    n_overlap_bcd = (uint8_t **) calloc(max_n_bin, sizeof(uint8_t*));
    if (NULL == n_overlap_bcd){
        fprintf(stderr, "ERROR! Failed to alloc memory for n_overlap_bcd!\n");
        exit(1);
    }
    for (int i = 0; i < max_n_bin; i++) {
        n_overlap_bcd[i] = (uint8_t *) calloc(max_n_bin, sizeof(uint8_t)); 
        if (NULL == n_overlap_bcd[i]){
            fprintf(stderr, "ERROR! Failed to alloc memory for n_overlap_bcd[%d]!\n", i);
            exit(1);
        }
    }
    
    all_frm_list = read_frm_list_from_bcd22_file(input_bcd22_file);
    fprintf(stderr, "number of frm: %lu\n", all_frm_list->size);

    for (tid1 = 0; tid1 < n_chr; tid1++){
        for (tid2 = 0; tid2 < n_chr; tid2++)
        {
            // if (tid1 != 11 || tid2 != 19) { continue; }
            fprintf(stderr, "tid1 = %d, tid2 = %d\n", tid1, tid2);
            fprintf(stderr, "clean-up memory\n");
            for (int i = 0; i < max_n_bin; i++)
            {
                for (int j = 0; j < max_n_bin; j++)
                {
                    n_overlap_bcd[i][j] = 0;
                }
            }
            fprintf(stderr, "grid overlap search started\n");
            two_chr_grid_overlap(tid1, tid2, chr_length_list, all_frm_list, bin_size, min_num_shared_barcodes, output_file, n_overlap_bcd);
            fprintf(stderr, "grid overlap search finished\n");
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
    FILE * out_fp;
    int min_num_shared_barcodes;

    if (argc < 5){
        usage(stderr);
        return 1;
    }

    input_bcd22_file = argv[1];
    output_file = argv[2];
    min_num_shared_barcodes = atoi(argv[3]);
    faidx_file = argv[4];

    if (min_num_shared_barcodes > UINT8_MAX-2){
        fprintf(stderr, "min_num_shared_barcodes cannot be larger than %d\n", UINT8_MAX-2);
        exit(1);
    }

    out_fp = fopen(output_file, "w");
    fclose(out_fp);

    whole_genome_grid_overlap(input_bcd22_file, output_file, (uint8_t) min_num_shared_barcodes, faidx_file);

    return 0;
}
