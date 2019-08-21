#ifndef TK_H 
#define TK_H

#define MAX_PATH_LENGTH 10240
#undef LINE_MAX
#define LINE_MAX 4096


#include <vector>
#include <string>
#include <limits>

#include <stdint.h>
#include "cgranges.h"


typedef struct {
	char ** data_list;
	int32_t size;
	int32_t capacity;
	int32_t max_string_length;
} STRING_LIST;

char * str_upper(const char * in_string);

STRING_LIST * init_string_list (int capacity, int32_t max_string_length);

int append_string_list (STRING_LIST * string_list, char * new_string);


typedef struct {
    int32_t * data_list;
    int32_t size;
    int32_t capacity;
} INT_LIST;

INT_LIST * init_int_list (int capacity);
int append_int_list (INT_LIST * int_list, int new_number);
int sort_int_list(INT_LIST * int_list);
INT_LIST * reset_int_list(INT_LIST * int_list);  // set the first int_list->size values to 0
INT_LIST * deep_reset_int_list(INT_LIST * int_list);  // set the first int_list->capacity values to 0
int free_int_list(INT_LIST * int_list);

typedef struct {
	int32_t ** data;
	int32_t size1;
	int32_t size2;
	int32_t capacity1;
	int32_t capacity2;
} INT_2D_LIST;

INT_2D_LIST * init_int_2d_list(int capacity1, int capacity2);
int free_int_2d_list(INT_2D_LIST * int_2d_list);

typedef struct {
	STRING_LIST * chrname_list;
	void * chrname_to_tid_hash;
	INT_LIST * chr_length_list; 
} CHR_INFO;

typedef struct {
	int64_t key1;
	int64_t key2;
	int64_t frag_id1;
	int64_t frag_id2;
}NODE;

typedef struct {
    NODE * data_list;
    int32_t size;
    int32_t capacity;
} NODE_LIST;

NODE_LIST * init_node_list (int capacity);

int append_node_list (NODE_LIST * node_list, NODE new_node);

INT_LIST * get_chr_length_from_faidx_file(const char * faidx_file);

typedef struct {
	int32_t tid1;
	int32_t start1;
	int32_t end1;
	int32_t tid2;
	int32_t start2;
	int32_t end2;
} BEDPE_CORE;

typedef struct {
	BEDPE_CORE * data_list;
	int32_t size;
	int32_t capacity;
} BEDPE_CORE_LIST;


BEDPE_CORE_LIST * init_bedpe_core_list(int capacity);
int append_bedpe_core_list (BEDPE_CORE_LIST * bedpe_core_list,  const BEDPE_CORE * bedpe_core);
BEDPE_CORE_LIST * read_bedpe_core_file(const char * bedpe_core_file, const char * faidx_file);



CHR_INFO * get_chr_info(const char * faidx_file); // example: chr_info = get_chr_info(faidx_file);
int chr_to_tid(const char * chrname, const CHR_INFO * chr_info); // example: tid1 = chr_to_tid(chr1, chr_info)
int chrname2tid(const char * chrname, const CHR_INFO * chr_info); // example: tid1 = chrname2tid(chr1, chr_info)
const char * tid2chrname(int tid, const CHR_INFO * chr_info); // example: chrname = tid2chrname(tid, chr_info);

uint32_t bcdseq2bcdint(const char * bcd_seq);
int bcdint2bcdseq (uint32_t bcd_int, char * bcd_seq);



FILE * open_file(const char * file, const char * mode);

std::vector <std::string> split_cstring_with_delimiter (char * cstring,  const char * delimiters);

class QuantileNumbers{
public:
    const double invalid_value = std::numeric_limits<double>::lowest();

    std::vector<double> q;
    std::vector<double> b;

    QuantileNumbers()
    {
        q.reserve(1001);
        b.reserve(1001);
        for (int i = 0; i < 1001; i++)
        {
            q.push_back(invalid_value);
            b.push_back( (double) i / 1000.0 ); 
        }
    }
    
    void print(FILE * out_fp, double shift_value) const
    {
        fprintf(out_fp, "0.1%% => %.4f\n", q[1] - shift_value);
        fprintf(out_fp, "1%% => %.4f\n", q[10] - shift_value); 
        fprintf(out_fp, "5%% => %.4f\n", q[50] - shift_value); 
        fprintf(out_fp, "10%% => %.4f\n", q[100] - shift_value); 
        fprintf(out_fp, "25%% => %.4f\n", q[250] - shift_value); 
        fprintf(out_fp, "50%% => %.4f\n", q[500] - shift_value); 
        fprintf(out_fp, "75%% => %.4f\n", q[750] - shift_value); 
        fprintf(out_fp, "90%% => %.4f\n", q[900] - shift_value); 
        fprintf(out_fp, "95%% => %.4f\n", q[950] - shift_value); 
        fprintf(out_fp, "99%% => %.4f\n", q[990] - shift_value); 
        fprintf(out_fp, "99.9%% => %.4f\n", q[999] - shift_value); 

    }
};

int calculate_distribution_from_count_vector(const std::vector<int> & input_count_vector, QuantileNumbers & quantile_numbers);

/* cgranges invertal tree functions  */


class Interval
{
public:
    char * ctg;
    int32_t tid ;
    int32_t start_pos;
    int32_t end_pos;
    void * data;
    Interval()
    {
        ctg = NULL;
        tid = -1;
        start_pos = -1;
        end_pos = -1;
        data = NULL; 
    };
};

int get_interval_vector_from_bed_file(const std::string & input_bed_file, const CHR_INFO * chr_info, std::vector <Interval> & db_interval_vector);
cgranges_t *generate_cr_interval_tree(const std::vector <Interval> & db_interval_vector);
int search_overlap_from_cr_interval_tree(cgranges_t * cr, const Interval & input_itv, std::vector <size_t> & output_index_vector);


#endif 
