#pragma once
#include <globaldef.h>
#include <stdio.h>

void parse_array(FILE *p_file, sfloat **p_parsed_arr, int *p_n, int *p_m);

void print_array(const sfloat *arr, const int n, const int m);

void pprint_array(const sfloat *p_arr, const int n, const int m, const char *string);

sfloat dot(const sfloat *p_vec1, const sfloat *p_vec2, const int n);

sfloat norm2(const sfloat *p_vec, const int n);

void add_vec(const sfloat *p_x, const sfloat a, const sfloat *p_v, sfloat *p_y, const int n);

void add_vec_inplace(sfloat *p_x, const sfloat a, const sfloat *p_v, const int n);

void add_2vec_inplace(sfloat *p_x, const sfloat a, const sfloat *p_u, const sfloat b, const sfloat *p_v, const int n);

void copy_vec(const sfloat *p_src, sfloat *p_dest, const int n);

void reset_vec(sfloat *vec, const int n);
