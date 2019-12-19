/* 
This library may for example be used for the coupling of cell values between locally "identical" mesh zones. For example to make a custom rebuild of the porous media heat transfer 
algorithm of ANSYS Fluent. For example to make certain additions which are not 
accessible from the close sources application. A simple example is provided in 
example.c.

An extension to exchange properties between surface areas may come in future versions
of this library.

Tested with ANSYS Fluent 2019 R3 in 3DDP mode with 4 cores.

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
*/

#ifndef _NN_MATCHING_H

#define _NN_MATCHING_H 1

#include"udf.h"
#include "global.h"
#include "stdlib.h"

/* Settings: */
#define MAX_COMP_NODE_COUNT 64  /*Modify for HPC?!*/
#define NN_TOL 1E-4             /* NN Tolerance, can be small for identical meshes...*/

int checkGetCellCountInNode(Thread *ct, int *c_count_node);


void calc_start_idx_arr(int *start_idx_arr, int N);


int getnodecidx(int *start_idx_arr, cell_t c);


int getcoordarr_fromallnodes(
                                Thread *ct, 
                                real (**x_c_arr)[ND_ND],
                                int *no_c_in_ct, 
                                int *c_global_node_start_idcs_arr
                             );


void det_nn_matching_arrs( 
                            real (*x_arr1)[ND_ND],
                            real (*x_arr2)[ND_ND],
                            int N_arr1,
                            int N_arr2,
                            int **nn_matching_arr_1_to_2,
                            int **nn_matching_arr_2_to_1
                        );


void nn_matching_zones(
                        Thread* ts, 
                        Thread* tf, 
                        int *N_cs, 
                        int *N_cf ,
                        int *node_start_idx_arr_cts,
                        int *node_start_idx_arr_ctf,
                        int **matching_sf,
                        int **matching_fs
                        );


void set_cell_val_exchng_arr(
                            Thread *ct, 
                            int *node_start_idx_arr, 
                            real *cell_val_exchng_arr, 
                            int N,
                            real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*) // test
                            );


void write_cell_value_exch_to_UDMI(
                                    Thread *ct, 
                                    int *node_start_idx_arr, 
                                    int *matching_arr,
                                    int N_match,
                                    real *cell_value_exchng_arr, 
                                    int N_T,
                                    int iUDMI
                                );


void cell_value_exchange_udmi_cellzone(
                                        Thread* ts, 
                                        Thread* tf, 
                                        int *node_start_idx_arr_cts,
                                        int *node_start_idx_arr_ctf,
                                        real *CVAL_cs,
                                        real *CVAL_cf,
                                        int *matching_sf,
                                        int *matching_fs,
                                        int N_s,
                                        int N_f,
                                        real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*),
                                        int iUDMI
                                        );


#endif
