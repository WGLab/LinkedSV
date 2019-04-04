#ifndef TK_H 
#define TK_H

#define MAX_READ_LENGTH 1048576
#define MAX_PATH_LENGTH 10240
#define LINE_MAX 4096

#include <stdint.h>

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
int chr_to_tid(const char * chrname, CHR_INFO * chr_info); // example: tid1 = chr_to_tid(chr1, chr_info)

uint32_t bcdseq2bcdint(const char * bcd_seq);
int bcdint2bcdseq (uint32_t bcd_int, char * bcd_seq);



FILE * open_file(const char * file, const char * mode);


#endif 
