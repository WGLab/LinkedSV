# include <stdio.h>
# include <stdlib.h>
# include <zlib.h>
# include <string.h>
# include <ctype.h>
# include "tk.h"


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


STRING_LIST * init_string_list(int capacity, size_t max_string_length)
{
    int i;

    STRING_LIST * new_string_list;
    new_string_list = (STRING_LIST *) calloc(1, sizeof(STRING_LIST));
    new_string_list->size = 0;
    new_string_list->capacity = capacity;
    new_string_list->max_string_length = max_string_length;
    new_string_list->data_list = (char **)calloc(new_string_list->capacity, sizeof(char *)); 

    if (new_string_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%zu)\n", capacity, max_string_length);
        exit(1);
    }
    for (i = 0; i < new_string_list->capacity; i++){
        new_string_list->data_list[i] = (char *)calloc(max_string_length, sizeof(char));
        if (new_string_list->data_list[i] == NULL){
            fprintf(stderr, "ERROR! Failed to alloc memory for string list (capacity=%d, max_string_length=%zu)\n", capacity, max_string_length);
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
            fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%zu, max_string_length=%zu)\n", string_list->capacity, string_list->max_string_length);
            exit(1);
        }
        for (i = string_list->size; i < string_list->capacity; i++)
        {
            string_list->data_list[i] = (char *)calloc(string_list->max_string_length, sizeof(char));
            if (string_list->data_list[i] == NULL){
                fprintf(stderr, "ERROR! Failed to realloc memory for string list (capacity=%zu, max_string_length=%zu)\n", string_list->capacity, string_list->max_string_length);
                exit(1);
            }
        }
    }
    strcpy(string_list->data_list[string_list->size], new_string);
    string_list->size += 1;

    return 0;

}

inline INT_LIST * init_int_list (int capacity)
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
            fprintf(stderr, "ERROR! Failed to realloc memory for int list (capacity=%zu)\n", int_list->capacity);
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
            fprintf(stderr, "ERROR! Failed to realloc memory for int list (capacity=%zu)\n", node_list->capacity);
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

