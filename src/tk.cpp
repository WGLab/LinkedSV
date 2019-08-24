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


int chrname2tid(const char * chrname, const CHR_INFO * chr_info)
{
    
    int tid;
	int ret;

	ret = khash_str2int_get(chr_info->chrname_to_tid_hash, chrname, &tid);
	if (ret == 0){
		return tid;
	}else{
		fprintf(stderr, "WARNING! %s is not found in hash map!\n", chrname);
		return -1;
	}
}

int chr_to_tid(const char * chrname, const CHR_INFO * chr_info)
{
	return chrname2tid(chrname, chr_info);
}

const char * tid2chrname(int tid, const CHR_INFO * chr_info)
{
    return chr_info->chrname_list->data_list[tid];
}

BEDPE_CORE_LIST * read_bedpe_core_file(const char * bedpe_core_file, const char * faidx_file)
{

	char * line;
	int32_t tid1, tid2;
	int32_t start1, end1, start2, end2;
	char * chr1, *chr2;

	BEDPE_CORE_LIST * bedpe_core_list;
	void * chrname_to_tid_hash;
	const CHR_INFO * chr_info;

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

INT_LIST * reset_int_list(INT_LIST * int_list) 
{
	//memset(int_list->data_list, 0, int_list->size * sizeof(int32_t) );
	for (int i = 0; i < int_list->size; i++) {
		int_list->data_list[i] = 0;
	}
	int_list->size = 0;
	return int_list;
}

INT_LIST * deep_reset_int_list(INT_LIST * int_list) 
{
	for (int i = 0; i < int_list->capacity; i++) {
		int_list->data_list[i] = 0;
	}
	int_list->size = 0;
	return int_list;
}

int free_int_list(INT_LIST * int_list)
{
    free(int_list->data_list);
    return 0;
}

int append_int_list (INT_LIST * int_list, int new_number)
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

uint32_t bcdseq2bcdint(const char * bcd_seq)
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

int bcdint2bcdseq (uint32_t bcd_int, char * bcd_seq) // bcd_seq must have at least 19 bytes spaces 
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


std::vector <std::string> split_cstring_with_delimiter (char * cstring, const char * delimiters)
{
    // this function will remove the '\n' and '\r' in the end of the string
    // the operation is directly done in the cstring
    std::vector <std::string> col_vector;
    char * pch;
    int str_len;
    str_len = strlen(cstring);
    int i;
    i = str_len-1;
    while (cstring[i] == '\n' || cstring[i] == '\0' || cstring[i] == '\r') // remove '\n' '\r' in the end 
    {   
        cstring[i] = '\0';
        i -= 1;
        if (i == -1)
        {
            break;
        }
    }

    pch = strtok (cstring, delimiters);
    while (pch != NULL)
    {   
        col_vector.push_back(pch);
        pch = strtok (NULL, delimiters);
    }   

    return col_vector;
}

/* cgranges invertal tree functions */

int get_interval_vector_from_bed_file(const std::string & input_bed_file, const CHR_INFO * chr_info, std::vector <Interval> & db_interval_vector)
{
    FILE * input_bed_fp;
    char * line, * ctg;
    Interval itv;

    line = new char [LINE_MAX];
    ctg = new char[LINE_MAX];
    input_bed_fp = fopen(input_bed_file.c_str(), "r");
    if (NULL == input_bed_fp)
    {
        fprintf(stderr, "Failed to open file: %s\n", input_bed_file.c_str());
        exit(1);
    }
    while (fgets(line, LINE_MAX, input_bed_fp))
    {
        sscanf (line, "%s\t%d\t%d\t%*s\n", ctg, &itv.start_pos, &itv.end_pos);
        itv.tid = chr_to_tid(ctg, chr_info);
        if (itv.tid < 0) { continue; }
        itv.ctg = chr_info->chrname_list->data_list[itv.tid];
        db_interval_vector.push_back(itv);
    }

    fclose(input_bed_fp);
    delete [] line;
    delete [] ctg;
    return 0;
}

cgranges_t * generate_cr_interval_tree(const std::vector <Interval> & db_interval_vector)
{
    cgranges_t * cr = cr_init();

    for (int i = 0; i < db_interval_vector.size(); i++)
    {
        cr_add(cr, db_interval_vector[i].ctg, db_interval_vector[i].start_pos, db_interval_vector[i].end_pos, i);
    }
    cr_index(cr); 

    // cr_label is the index in the db_interval_vector
    return cr; 
}

int search_overlap_from_cr_interval_tree(cgranges_t * cr, const Interval & input_itv, std::vector <size_t> & output_index_vector)
{
    // cr should be alreadly generated by generate_cr_interval_tree()

    int64_t i, num_overlap, *b = 0, max_b = 0;
    Interval ovl_itv;
    size_t index;

    num_overlap = cr_overlap (cr, input_itv.ctg, input_itv.start_pos, input_itv.end_pos, &b, &max_b); // overlap query; output array b[] can be reused
    output_index_vector.clear();
    for (i = 0; i < num_overlap; ++i) 
    {
        
        ovl_itv.ctg = input_itv.ctg;
        ovl_itv.tid = input_itv.tid;
        ovl_itv.start_pos = cr_start(cr, b[i]);
        ovl_itv.end_pos = cr_end(cr, b[i]);
        index = cr_label(cr, b[i]);
        output_index_vector.push_back(index);
    }
    free(b);
    return num_overlap; 
}



int calculate_distribution_from_count_vector(const std::vector<int> & input_count_vector, QuantileNumbers & quantile_numbers)
{
    const double invalid_value = std::numeric_limits<double>::lowest();

    double total_count; 
    double accu_count; 
    double fraction; 

    total_count = 0;
    for (auto & cnt : input_count_vector) {
        total_count += cnt;
    }
    accu_count = 0;
    for (int i = 0; i < input_count_vector.size(); i++)
    {
        accu_count += input_count_vector[i];
        fraction = accu_count / total_count; 
        for (int j = 0; j < 1001; j++)
        {
            if (quantile_numbers.q[j] == quantile_numbers.invalid_value && fraction >= quantile_numbers.b[j])
            {
                quantile_numbers.q[j] = i;
            }
        }
    }
    return 0;
}


int remove_file(const char * file_name)
{
    if (remove( file_name ) != 0)
    {
        fprintf(stderr, "Error deleting file: %s\n", file_name);
    }
    return 0;
}

