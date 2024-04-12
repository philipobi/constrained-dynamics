#pragma once
#include <globaldef.h>
#include <stdio.h>

void parse_array (FILE *p_file, sfloat **p_parsed_arr, int *p_n, int *p_m);

void print_array (sfloat *arr, int n, int m);