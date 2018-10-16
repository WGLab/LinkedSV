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

INT_LIST * init_int_list (int capacity)
{
    INT_LIST * new_int_list;
    new_int_list = (INT_LIST *) calloc (1, sizeof(INT_LIST));
    new_int_list->size = 0;
    new_int_list->capacity = capacity;
    new_int_list->data_list = (int *) calloc (new_int_list->capacity, sizeof(int)); 

    if (new_int_list->data_list == NULL){
        fprintf(stderr, "ERROR! Failed to alloc memory for int list (capacity=%d)\n", capacity);
        exit(1);
    }   

    return new_int_list;

}

int append_int_list (INT_LIST * int_list, int new_number)
{
    if (int_list->capacity == int_list->size)
    {   
        int_list->capacity  = (int)(int_list->capacity * 1.2);
        int_list->data_list = (int*) realloc (int_list->data_list, int_list->capacity*sizeof(int));
        if (int_list->data_list == NULL){
            fprintf(stderr, "ERROR! Failed to realloc memory for int list (capacity=%d)\n", int_list->capacity);
            exit(1);
        }   
    }   
    int_list->data_list[int_list->size] = new_number;  
    int_list->size += 1;

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
