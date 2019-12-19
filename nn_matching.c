/* 
Description see "nn_matching.h"

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
*/

#include "nn_matching.h"


int checkGetCellCountInNode(Thread *ct, int *c_count_node)
{
/*  get cells count in current compute node and checks certain assumptions 
    for cell (c) numbering on different processors
    which should be valid since Fluent 2019 ... */

    cell_t c;
    int ci_min = INT_MAX;
    int ci_max = INT_MIN;
    (*c_count_node) = 0;

    begin_c_loop_int(c, ct) 
    {
        ci_min = MIN(ci_min, (*c_count_node));
        ci_max = MAX(ci_max, c);
        (*c_count_node) = (*c_count_node) + 1;
    }
    end_c_loop_int(c, ct)

    if ((*c_count_node) == 0 || (ci_min == 0 && ci_max == (*c_count_node) - 1))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


void calc_start_idx_arr(int *start_idx_arr, int N)
{
/*  Description
    [2,3,3] -> sizes
    [2,5,8] -> cum sumes
    [ 0, 1 ,2! ,3,4, 5!,6,7 -> consolidated array numbering
    [[0,1],[0,1,2],[0,1,2]] -> Arrays
    [0 ,2 , 5] -> this should be in start_idx_arr ... */
    int i;

    for(i=1; i<N; ++i)
    {
        start_idx_arr[i] += start_idx_arr[i-1];
    }

    for(i=(N-1); i>0;--i)
    {
        start_idx_arr[i] = start_idx_arr[i-1];
    }
     
    start_idx_arr[0] = 0;
}


int getnodecidx(int *start_idx_arr, cell_t c)
{
   return start_idx_arr[myid] + (int) c;
}


int getcoordarr_fromallnodes(
                                Thread *ct, 
                                real (**x_c_arr)[ND_ND],
                                int *no_c_in_ct, 
                                int *c_global_node_start_idcs_arr
                             )
{
int status = 0;
cell_t c;

real *rwork;
int iwork[MAX_COMP_NODE_COUNT]; 
int i_c;

real x[ND_ND];

real* x_c;
real* y_c;

#if RP_3D
real* z_c;
#endif

status = checkGetCellCountInNode (ct, no_c_in_ct);

if(myid >= MAX_COMP_NODE_COUNT)
{
    status = 0;
    Message("Error getcoordarr_fromallnodes(): Myid %i greater than max "
            "specified compute node count %i.\n", myid, MAX_COMP_NODE_COUNT);
}


if (status)
{
    c_global_node_start_idcs_arr[myid] = (int) (*no_c_in_ct); // error why?
    (*no_c_in_ct) = PRF_GISUM1((*no_c_in_ct));

    PRF_GISUM(c_global_node_start_idcs_arr, MAX_COMP_NODE_COUNT, iwork);

    calc_start_idx_arr(c_global_node_start_idcs_arr, compute_node_count);

    *x_c_arr = (real (*)[ND_ND]) calloc(ND_ND * (*no_c_in_ct), sizeof(real));

    x_c = (real *) calloc((*no_c_in_ct), sizeof(real));
    y_c = (real *) calloc((*no_c_in_ct), sizeof(real));
    z_c = (real *) calloc((*no_c_in_ct), sizeof(real));

    begin_c_loop_int(c, ct)
    {
        i_c = getnodecidx(c_global_node_start_idcs_arr, c);
        C_CENTROID(x, c, ct);

        if(i_c >= 0 && i_c < (*no_c_in_ct))
        {
            #if RP_3D
            ND_SET(x_c[i_c] , y_c[i_c], z_c[i_c], x[0], x[1], x[2]);
            #else
            ND_SET(x_c[i_c] , y_c[i_c], x[0], x[1]);
            #endif
        }
        else
        {
            Message("error i_c\n");
        }
        
    }
    end_c_loop_int(c, ct)

    rwork = (real *) calloc((*no_c_in_ct), sizeof(real)); 
    PRF_GRSUM(x_c, (*no_c_in_ct), rwork);
    PRF_GRSUM(y_c, (*no_c_in_ct), rwork);
    #if RP_3D
    PRF_GRSUM(z_c, (*no_c_in_ct), rwork);
    #endif

    free(rwork);

    for (i_c = 0; i_c<(*no_c_in_ct); ++i_c)
    {
        (*x_c_arr)[i_c][0] = x_c[i_c];
        (*x_c_arr)[i_c][1] = y_c[i_c];
        #if RP_3D
        (*x_c_arr)[i_c][2] = z_c[i_c];
        #endif
    }

    status = 1;
    Message0("Coordinates for cell thread ct: %i read!\n", THREAD_ID(ct));
}

return status;
}


void det_nn_matching_arrs( 
                            real (*x_arr1)[ND_ND],
                            real (*x_arr2)[ND_ND],
                            int N_arr1,
                            int N_arr2,
                            int **nn_matching_arr_1_to_2,
                            int **nn_matching_arr_2_to_1
                        )
{
int *iwork;

free(*nn_matching_arr_1_to_2);
(*nn_matching_arr_1_to_2) = (int *) calloc(N_arr1, sizeof(int));

free(*nn_matching_arr_2_to_1);
(*nn_matching_arr_2_to_1) = (int *) calloc(N_arr2, sizeof(int));

if (myid == 0)
{
/*  only serriel - parallelisation should be relativ easy 
    just distribute 1 coord array evenly to all processors and GRLOW after
    nnmin search
    but is not necessary right now ... */

    nn_idx_1_to_2(         
                    x_arr1,
                    x_arr2,
                    N_arr1,
                    N_arr2,
                    (*nn_matching_arr_1_to_2)
                 );

    nn_idx_1_to_2(         
                    x_arr2,
                    x_arr1,
                    N_arr2,
                    N_arr1,
                    (*nn_matching_arr_2_to_1)
                 );

}

iwork = (int *)calloc(N_arr1, sizeof(int));
PRF_GISUM((*nn_matching_arr_1_to_2), N_arr1, iwork);
free(iwork);

iwork = (int *)calloc(N_arr2, sizeof(int));
PRF_GISUM((*nn_matching_arr_2_to_1), N_arr2, iwork);
free(iwork);
}


void nn_idx_1_to_2(         
                    real (*x_arr1)[ND_ND],
                    real (*x_arr2)[ND_ND],
                    int N_arr1,
                    int N_arr2,
                    int *nn_matching_arr
                   )
{
/* Brute force NN search mapping indices of nearest values of arr2 to arr1
mapping 1 -> 2 */
int i, ii;
real distance;
real min_distance;

real x1[ND_ND];
real x2[ND_ND];
real x1_x2[ND_ND];

for(i = 0; i<N_arr1; ++i)
{
    min_distance = INT_MAX;
    NV_V(x1, = , x_arr1[i]);

    for(ii = 0; ii <N_arr2; ++ii)
    {
        NV_V(x2, = , x_arr2[ii]);

        NV_VV(x1_x2, = , x2, -, x1);

        if (NV_MAG2(x1_x2) < min_distance)
        {
            min_distance = NV_MAG2(x1_x2);
            nn_matching_arr[i] = ii;
        }
    }

    if ( min_distance > NN_TOL)
    {
        Message("Warning det_nn_matching_arr(): Matched Cell distance "
        "above NN_TOL of %lE", NN_TOL);
    }
}
}


void nn_matching_zones(
                        Thread* ts, 
                        Thread* tf, 
                        int *N_cs, 
                        int *N_cf ,
                        int *node_start_idx_arr_cts,
                        int *node_start_idx_arr_ctf,
                        int **matching_sf,
                        int **matching_fs
                       )
{
int * iwork;
int status = 0;

real (*x_cs)[ND_ND] = NULL;
real (*x_cf)[ND_ND] = NULL;

status = getcoordarr_fromallnodes(  
                                    ts, 
                                    &x_cs, 
                                    N_cs, 
                                    node_start_idx_arr_cts
                                );

if(status)
{
    status = getcoordarr_fromallnodes(  
                                        tf, 
                                        &x_cf, 
                                        N_cf, 
                                        node_start_idx_arr_ctf
                                    );
}
else
{
    free(x_cs);
}

if(status)
{
    Message0("Start NN Matching ...\n");
    det_nn_matching_arrs(
                        x_cs, 
                        x_cf,
                        (*N_cs),
                        (*N_cf),
                        &(*matching_sf),
                        &(*matching_fs)
                    );

    free(x_cf);
    free(x_cs);

    Message0("NN Matching ready!\n");
}
else
{
    free(x_cs);
    free(x_cf);
}

}



void set_cell_val_exchng_arr(
                            Thread *ct, 
                            int *node_start_idx_arr, 
                            real *cell_val_exchng_arr, 
                            int N,
                            real (*C_VAL_WRAPPER_FUN)(cell_t, Thread*) 
                            )
{
	cell_t c;
    int i_c;
	real* rwork;
    int error = 0;

    /* reset to zero (for debuging)*/
    for (i_c = 0; i_c < N; ++i_c)
	{
		cell_val_exchng_arr[i_c] = 0.0;
	}

    begin_c_loop_int(c, ct) 
	{
        i_c = getnodecidx(node_start_idx_arr, c);
        if (i_c < N)
        {
            cell_val_exchng_arr[i_c] = (*C_VAL_WRAPPER_FUN)(c,ct);
        }
        else
        {
            error = 1;

        }   
	}
	end_c_loop_int(c, ct)

    if(error)
    {
        Message("Error set_cell_val_exchng_arr(): Out of bounds error\n");
    }

    rwork = (real*)calloc(N, sizeof(real));
    PRF_GRSUM(cell_val_exchng_arr, N, rwork);
    free(rwork);
}



void write_cell_value_exch_to_UDMI(
                                    Thread *ct, 
                                    int *node_start_idx_arr, 
                                    int *matching_arr,
                                    int N_match,
                                    real *cell_value_exchng_arr, 
                                    int N_T,
                                    int iUDMI
                                )
{
    cell_t c;
    int error_i_match = 0;
    int error_i_T = 0;
    int i_c1;
    int i_c2;

    if (N_UDM  >= (iUDMI+1))
    {
        begin_c_loop_int(c, ct) 
        {
            i_c1 = getnodecidx(node_start_idx_arr, c);

            if(i_c1 < N_match)
            {
                i_c2 = matching_arr[i_c1];

                    if(i_c2 < N_T)
                    {
                        C_UDMI(c, ct, iUDMI) = cell_value_exchng_arr[i_c2];
                    }
                    else
                    {
                        error_i_T = 1;
                    }
            }
            else
            {
                error_i_match = 1;
            }
        }
        end_c_loop_int(c, ct)

        if (error_i_T)
        {
            Message("Error write_T_exch_to_UDMI(): Overflow cell_value_exchng_arr indexing\n");
        }
        if (error_i_match)
        {
            Message("Error write_T_exch_to_UDMI(): Overflow matching_arr indexing\n");
        }
    }
    else
    {
        Message0("Error in write_T_exch_to_UDMI(): Check set UDMI values "
                 "smaller than %i, check Fluent settings!\n" , (iUDMI+1));
    }
}



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
                                        )
{
    set_cell_val_exchng_arr(
                            ts, 
                            node_start_idx_arr_cts, 
                            CVAL_cs, 
                            N_s,
                            C_VAL_WRAPPER_FUN
                            );

    set_cell_val_exchng_arr(
                            tf, 
                            node_start_idx_arr_ctf, 
                            CVAL_cf, 
                            N_f,
                            C_VAL_WRAPPER_FUN
                           );

    write_cell_value_exch_to_UDMI(   
                                    tf, 
                                    node_start_idx_arr_ctf, 
                                    matching_fs,
                                    N_f,
                                    CVAL_cs, 
                                    N_s, 
                                    iUDMI
                                );

    write_cell_value_exch_to_UDMI(   
                                    ts, 
                                    node_start_idx_arr_cts, 
                                    matching_sf,
                                    N_s,
                                    CVAL_cf, 
                                    N_f, 
                                    iUDMI
                                );
}
