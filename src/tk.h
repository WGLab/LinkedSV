#ifndef TK_H 
#define TK_H

#define MAX_READ_LENGTH 1048576
#define MAX_PATH_LENGTH 10240
#define LINE_MAX 4096

typedef struct {
	char ** data_list;
	size_t size;
	size_t capacity;
	size_t max_string_length;
} STRING_LIST;

char * str_upper(const char * in_string);

STRING_LIST * init_string_list (int capacity, size_t max_string_length);

int append_string_list (STRING_LIST * string_list, char * new_string);



typedef struct {
    int32_t * data_list;
    size_t size;
    size_t capacity;
} INT_LIST;

INT_LIST * init_int_list (int capacity);
int append_int_list (INT_LIST * int_list, int new_number);
int sort_int_list(INT_LIST * int_list);
INT_LIST * reset_int_list(INT_LIST * int_list);  // set the first int_list->size values to 0
INT_LIST * deep_reset_int_list(INT_LIST * int_list);  // set the first int_list->capacity values to 0
int free_int_list(INT_LIST * int_list);


typedef struct {
	int64_t key1;
	int64_t key2;
	int64_t frag_id1;
	int64_t frag_id2;
}NODE;

typedef struct {
    NODE * data_list;
    size_t size;
    size_t capacity;
} NODE_LIST;

NODE_LIST * init_node_list (int capacity);

int append_node_list (NODE_LIST * node_list, NODE new_node);

INT_LIST * get_chr_length_from_faidx_file(const char * faidx_file);
#endif 
