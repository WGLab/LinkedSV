# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>
# include <htslib/khash_str2int.h>

# include "tk.h"


BEDPE_CORE_LIST * init_bedpe_core_list(int capacity)
{
	BEDPE_CORE_LIST * new_bedpe_core_list;
	new_bedpe_core_list = (BEDPE_CORE_LIST *) calloc(1, sizeof(BEDPE_CORE_LIST));
	new_bedpe_core_list->size = 0;
	new_bedpe_core_list->capacity = capacity;
	new_bedpe_core_list->data_list = (BEDPE_CORE *) calloc(new_bedpe_core_list->capacity, sizeof(BEDPE_CORE));

    if (new_bedpe_core_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for BEDPE_CORE_LIST (capacity=%d)\n", new_bedpe_core_list->capacity);
        exit(1);
    }   

    return new_bedpe_core_list; 

}

int append_bedpe_core_list (BEDPE_CORE_LIST * bedpe_core_list,  const BEDPE_CORE * bedpe_core)
{
    if (bedpe_core_list->capacity == bedpe_core_list->size)
    {   
        bedpe_core_list->capacity  = (int)(bedpe_core_list->capacity * 1.5);
        bedpe_core_list->data_list = (BEDPE_CORE *) realloc (bedpe_core_list->data_list, bedpe_core_list->capacity*sizeof(BEDPE_CORE));
        if (bedpe_core_list->data_list == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for bedpe_core_list (capacity=%d)\n", bedpe_core_list->capacity);
            exit(1);
        }   
    }   
	memcpy(bedpe_core_list->data_list + bedpe_core_list->size, bedpe_core, sizeof(BEDPE_CORE));
    bedpe_core_list->size += 1;

    return 0;
}


CHR_INFO * get_chr_info(const char * faidx_file)
{
	int tid;
	FILE * faidx_fp;
	char * line;
	char * chrname;
	CHR_INFO * chr_info;

	line = (char *) calloc (LINE_MAX, sizeof(char));
	chrname = (char *) calloc (LINE_MAX, sizeof(char));

	chr_info = (CHR_INFO *) calloc (1, sizeof(CHR_INFO));
	chr_info->chrname_list = init_string_list (1000, LINE_MAX);
	chr_info->chrname_to_tid_hash = khash_str2int_init();
	chr_info->chr_length_list = get_chr_length_from_faidx_file(faidx_file);

	faidx_fp = fopen(faidx_file, "r");
	if (faidx_fp == NULL){
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", faidx_file);
		exit(1);
	}

	while (fgets(line, LINE_MAX, faidx_fp))
	{
		if (line[0] == '#') { continue; }
		sscanf(line, "%s\t%*s\n", chrname);
		append_string_list (chr_info->chrname_list, chrname);
		tid = chr_info->chrname_list->size - 1;
		khash_str2int_set(chr_info->chrname_to_tid_hash, chr_info->chrname_list->data_list[tid], tid);
	}

	fclose(faidx_fp);
	free(line);
	free(chrname);

	return chr_info;
}


int chr_to_tid(const char * chrname, CHR_INFO * chr_info)
{
	int tid;
	int ret;

	ret = khash_str2int_get(chr_info->chrname_to_tid_hash, chrname, &tid);
	if (ret == 0){
		return tid;
	}else{
		fprintf(stderr, "ERROR! %s is not found in hash map!\n", chrname);
		return -1;
	}
}

BEDPE_CORE_LIST * read_bedpe_core_file(const char * bedpe_core_file, const char * faidx_file)
{

	char * line;
	int32_t tid1, tid2;
	int32_t start1, end1, start2, end2;
	char * chr1, *chr2;

	BEDPE_CORE_LIST * bedpe_core_list;
	void * chrname_to_tid_hash;
	CHR_INFO * chr_info;

	gzFile bedpe_core_fp;	
	BEDPE_CORE * bedpe_core;

	line = (char *) calloc (LINE_MAX, sizeof(char));
	chr1 = (char *) calloc (LINE_MAX, sizeof(char));
	chr2 = (char *) calloc (LINE_MAX, sizeof(char));
	bedpe_core = (BEDPE_CORE *) calloc(1, sizeof(BEDPE_CORE));

	bedpe_core_list = init_bedpe_core_list(100);

	chr_info = get_chr_info(faidx_file);

	bedpe_core_fp = gzopen(bedpe_core_file, "r");
	if (Z_NULL == bedpe_core_fp) {
		fprintf(stderr, "ERROR! Failed to open file for reading: %s\n", bedpe_core_file);
		exit(1);
	}

	while (gzgets(bedpe_core_fp, line, LINE_MAX))
	{
		if (line[0] == '#') { continue; }
		sscanf(line, "%s\t%d\t%d\t%s\t%d\t%d\t%*s\n", chr1, &bedpe_core->start1, &bedpe_core->end1, chr2, &bedpe_core->start2, &bedpe_core->end2);
		if (bedpe_core->start1 < 0 || bedpe_core->start2 < 0){
			fprintf(stderr, "WARNING! Skipped invalid input. start position < 0. The record is: %s", line);
			continue;
		}
		if (bedpe_core->start1 > bedpe_core->end1 || bedpe_core->start2 > bedpe_core->end2){
			fprintf(stderr, "WARNING! Skiped invalid input. start position > end position. The record is: %s", line);
			exit(1);
		}
		bedpe_core->tid1 = chr_to_tid(chr1, chr_info);
		bedpe_core->tid2 = chr_to_tid(chr2, chr_info);
		append_bedpe_core_list(bedpe_core_list, bedpe_core);
	}

	gzclose(bedpe_core_fp);
	free(line);
	free(chr1);
	free(chr2);
	free(bedpe_core);

	return bedpe_core_list;

}

char * str_upper(const char * in_string) 
// return a new string in which all the letters are in upper case
{
    int l;
    int i;
    l = strlen(in_string);

    char * new_str;

    new_str = (char *) calloc(l+1, sizeof(char));

    for (i = 0; i < l; i++)
    {
        new_str[i] = toupper(in_string[i]);
    }
    return new_str;
}


STRING_LIST * init_string_list(int capacity, int32_t max_string_length)
{
    int i;

    STRING_LIST * new_string_list;
    new_string_list = (STRING_LIST *) calloc(1, sizeof(STRING_LIST));
    new_string_list->size = 0;
    new_string_list->capacity = capacity;
    new_string_list->max_string_length = max_string_length;
    new_string_list->data_list = (char **)calloc(new_string_list->capacity, sizeof(char *)); 

    if (new_string_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%d)\n", capacity, max_string_length);
        exit(1);
    }
    for (i = 0; i < new_string_list->capacity; i++){
        new_string_list->data_list[i] = (char *)calloc(max_string_length, sizeof(char));
        if (new_string_list->data_list[i] == NULL){
            fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%d)\n", capacity, max_string_length);
            exit(1);
        }
    }

    return new_string_list;

}

int append_string_list (STRING_LIST * string_list, char * new_string)
{
    int l;
    int i;
    l = strlen(new_string);

    if (l > string_list->max_string_length-1){
        fprintf(stderr, "ERROR! Cannot append new string to string list! Too long\n");
        exit(1);
    }
    if (string_list->capacity == string_list->size)
    {
        string_list->capacity  = string_list->capacity * 2;
        string_list->data_list = (char **) realloc (string_list->data_list, string_list->capacity*sizeof(char*));
        if (string_list->data_list[i] == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%d, max_string_length=%d)\n", string_list->capacity, string_list->max_string_length);
            exit(1);
        }
        for (i = string_list->size; i < string_list->capacity; i++)
        {
            string_list->data_list[i] = (char *)calloc(string_list->max_string_length, sizeof(char));
            if (string_list->data_list[i] == NULL){
                fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%d, max_string_length=%d)\n", string_list->capacity, string_list->max_string_length);
                exit(1);
            }
        }
    }
    strcpy(string_list->data_list[string_list->size], new_string);
    string_list->size += 1;

    return 0;

}

INT_2D_LIST * init_int_2d_list(int capacity1, int capacity2)
{
	INT_2D_LIST * new_int_2d_list;
	new_int_2d_list = (INT_2D_LIST *) calloc(1, sizeof(INT_2D_LIST));
	new_int_2d_list->size1 = 0;
	new_int_2d_list->size2 = 0;
	new_int_2d_list->capacity1 = capacity1;
	new_int_2d_list->capacity2 = capacity2;

	new_int_2d_list->data = (int32_t **) calloc(new_int_2d_list->capacity1, sizeof(int32_t *));
	if (new_int_2d_list->data == NULL) {
        fprintf(stderr, "ERROR! Failed to alloc memory for INT_2D_LIST (capacity1=%d, capacity2=%d)\n", capacity1, capacity2);
        exit(1);
	}
	for (int32_t i = 0; i < new_int_2d_list->capacity1; i++)
	{
		new_int_2d_list->data[i] = (int32_t *) calloc (new_int_2d_list->capacity2, sizeof(int32_t));
		if (new_int_2d_list->data[i] == NULL) {
			fprintf(stderr, "ERROR! Failed to alloc memory for INT_2D_LIST (capacity1=%d, capacity2=%d)\n", capacity1, capacity2);
			exit(1);
		}
	}

	return new_int_2d_list;
}

int free_int_2d_list(INT_2D_LIST * int_2d_list)
{
	for (int32_t i = 0; i < int_2d_list->capacity1; i++) {
		free(int_2d_list->data[i]);
	}
	free(int_2d_list);
	return 0;
}



INT_LIST * init_int_list (int capacity)
{
    INT_LIST * new_int_list;
    new_int_list = (INT_LIST *) calloc (1, sizeof(INT_LIST));
    new_int_list->size = 0;
    new_int_list->capacity = capacity;
    new_int_list->data_list = (int32_t *) calloc (new_int_list->capacity, sizeof(int32_t)); 

    if (new_int_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for int list (capacity=%d)\n", capacity);
        exit(1);
    }   

    return new_int_list;

}

inline INT_LIST * reset_int_list(INT_LIST * int_list) 
{
	//memset(int_list->data_list, 0, int_list->size * sizeof(int32_t) );
	for (int i = 0; i < int_list->size; i++) {
		int_list->data_list[i] = 0;
	}
	int_list->size = 0;
	return int_list;
}

inline INT_LIST * deep_reset_int_list(INT_LIST * int_list) 
{
	for (int i = 0; i < int_list->capacity; i++) {
		int_list->data_list[i] = 0;
	}
	int_list->size = 0;
	return int_list;
}

inline int free_int_list(INT_LIST * int_list)
{
    free(int_list->data_list);
    return 0;
}

inline int append_int_list (INT_LIST * int_list, int new_number)
{
    if (int_list->capacity == int_list->size)
    {   
        int_list->capacity  = (int)(int_list->capacity * 1.2);
        int_list->data_list = (int32_t *) realloc (int_list->data_list, int_list->capacity*sizeof(int32_t));
        if (int_list->data_list == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for int list (capacity=%d)\n", int_list->capacity);
            exit(1);
        }   
    }   
    int_list->data_list[int_list->size] = new_number;  
    int_list->size += 1;

    return 0;
}

int cmpfunc_int32(const void * a, const void * b)
{
	 return ( *(int32_t*)a - *(int32_t*)b );
}

int sort_int_list(INT_LIST * int_list)
{
	qsort(int_list->data_list, int_list->size, sizeof(int32_t), cmpfunc_int32);
	return 0;
}



NODE_LIST * init_node_list (int capacity)
{
    NODE_LIST * new_node_list;
    new_node_list = (NODE_LIST *) calloc (1, sizeof(NODE_LIST));
    new_node_list->size = 0;
    new_node_list->capacity = capacity;
    new_node_list->data_list = (NODE *) calloc (new_node_list->capacity, sizeof(NODE)); 

    if (new_node_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for node list (capacity=%d)\n", capacity);
        exit(1);
    }   

    return new_node_list;

}

int append_node_list (NODE_LIST * node_list, NODE new_node)
{
    if (node_list->capacity == node_list->size)
    {
        if ((int)(node_list->capacity * 1.2) > 100){
            node_list->capacity  = (int)(node_list->capacity * 1.2);
        }else{
            node_list->capacity  = (int)(node_list->capacity + 1000);
        }
        node_list->data_list = (NODE *) realloc (node_list->data_list, node_list->capacity*sizeof(NODE));
        if (node_list->data_list == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for int list (capacity=%d)\n", node_list->capacity);
            exit(1);
        }   
    }   
    node_list->data_list[node_list->size] = new_node;  
    node_list->size += 1;

    return 0;
}

INT_LIST * get_chr_length_from_faidx_file(const char * faidx_file)
{
    INT_LIST *chr_length_list;
    char * line = NULL;
    int chr_length = 0;
    FILE * faidx_fp;
    // char * token;
    faidx_fp = fopen(faidx_file, "r");
    if (NULL == faidx_fp){
        fprintf(stderr, "ERROR! Failed to open file: %s\n", faidx_file);
        exit(1);
    }   

    chr_length_list = init_int_list(100);
    line = (char *) calloc (LINE_MAX, sizeof(char)); 

    while (fgets(line, LINE_MAX, faidx_fp)) {
        sscanf(line, "%*s\t%d\t%*s\n", &chr_length);    
        append_int_list(chr_length_list, chr_length);
    }   

    fclose(faidx_fp);
    free(line);
    return chr_length_list;
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

FILE * open_file(const char * file, const char * mode)
{
	FILE * fp;

	fp = NULL;
	fp = fopen(file, mode);
	if (NULL == fp){
		fprintf(stderr, "ERROR! Failed to open file: %s\n", file);
		exit(1);
	}

	return fp;

}

