#pragma once
#include <globaldef.h>
#include <stdio.h>

void parse_array (FILE *p_file, sfloat **p_parsed_arr, int *p_n, int *p_m);

void print_array (sfloat *arr, int n, int m);

sfloat dot (const sfloat *p_vec1, const sfloat *p_vec2, const int n);

sfloat norm2 (const sfloat *p_vec, const int n);

void add_vec (const sfloat *p_vec1, const sfloat *p_vec2, const int n,
              sfloat *p_vec_result);

void mul_vec (const sfloat *p_vec, const sfloat a, const int n,
              sfloat *p_vec_result);

void sub_vec (const sfloat *p_vec1, const sfloat *p_vec2, const int n,
              sfloat *p_vec_result);