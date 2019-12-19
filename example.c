#include "udf.h"
#include "nn_matching.h"
/* example.c
         - simple heat transfer fluid-solid porous media (pm) coupling


To reconstruct this example:
- simply copy the a pm fluid zone and make the "copy" a solid zone
- setup solid and fluid id CELL_ID_SOLID/CELL_ID_FLUID
- setup interface area between solid and fluid (As_sf)
- setup heat transfer between solid and fluid (alpha_sf)


Example may be extended to multiple zones, moving zones, other exchange values, 
more complex alpha functions etc...

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


/* ---- SETTINGS */

#define EZNO 1 /* number of coupled solid - fluid regions */
#define CELL_ID_SOLID 11
#define CELL_ID_FLUID 63

real As_sf = 1; /* Contact Surface Solid / Fluid */
real alpha_sf = 40;  /* Alpha Solid / Fluid */

/* ---- END SETTINGS */


/* ---- GLOBAL COUPLING VARIABLES */

int g_N_cf[EZNO] = {0};    /* not sure if {0} works - on all compilers 
                            see https://stackoverflow.com/questions/14797810///why-are-my-structs-members-not-properly-initialised-using? */ 
int g_N_cs[EZNO] = {0};

int g_node_start_idx_arr_ctf[EZNO][MAX_COMP_NODE_COUNT] = {0};
int g_node_start_idx_arr_cts[EZNO][MAX_COMP_NODE_COUNT] = {0};

int **g_matching_fs = NULL;
int **g_matching_sf = NULL;

real **g_T_cf = NULL;
real **g_T_cs = NULL;
/* ---- END GLOBAL COUPLING VARIABLES */



/* Function Definitions */

void calloc_double_pointer_globals()
{
    int i,j;

    for(i = 0; i<MAX_COMP_NODE_COUNT; ++i)
    {
        for(j=0; j<EZNO; ++j)
        {
              g_node_start_idx_arr_ctf[j][i] = 0; /* just to be sure ... */
              g_node_start_idx_arr_cts[j][i] = 0;
        }
    }

	g_matching_fs = (int **) calloc(EZNO, sizeof(int *));
	g_matching_sf = (int **) calloc(EZNO, sizeof(int *));

	g_T_cf = (real **) calloc(EZNO, sizeof(real *));
	g_T_cs = (real **) calloc(EZNO, sizeof(real *));
}


void alloc_global_T_exch_arrays(int ezid)
{
    free(g_T_cf[ezid]);
    free(g_T_cs[ezid]);
    g_T_cf[ezid] = (real *)calloc(g_N_cf[ezid], sizeof(real));
    g_T_cs[ezid] = (real *)calloc(g_N_cs[ezid], sizeof(real));
}



void zero_udmis(int iUDMI)
{
#if !RP_HOST
Domain *domain = Get_Domain(1);
Thread *c_thread;
cell_t c;

if (N_UDM >= (iUDMI+1))
{
    thread_loop_c(c_thread, domain)
    {
        begin_c_loop(c, c_thread) 
        {
                C_UDMI(c,c_thread,iUDMI) = 0;
        }
    end_c_loop(c, c_thread) 
    } 
}
else
{
    Message0("Error in zero_udmis(): Check set UDMI values "
             "smaller than %i, check Fluent settings!\n" , (iUDMI+1));
}
#endif
}


real CELL_VAL_FUN_WRAPPER(cell_t c, Thread *ct)
{
   return C_T(c, ct);
}


/* End Functions */



/* Start - Fluent DEF Functions */

DEFINE_ON_DEMAND(NN_Matching_oD)
{
    #if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *ts = Lookup_Thread(d, CELL_ID_SOLID);
    Thread *tf = Lookup_Thread(d, CELL_ID_FLUID);

	zero_udmis(0);
	calloc_double_pointer_globals();

    nn_matching_zones(
                        ts, 
                        tf,
                        &g_N_cs[0],
                        &g_N_cf[0],
                        &g_node_start_idx_arr_cts[0][0],
                        &g_node_start_idx_arr_ctf[0][0],
                        &g_matching_sf[0],
                        &g_matching_fs[0]
                     );

    alloc_global_T_exch_arrays(0);

    #endif
}


DEFINE_ADJUST(temperature_exchange, d)
{
    #if !RP_HOST
    Thread* ts = Lookup_Thread(d, CELL_ID_SOLID);
	Thread* tf = Lookup_Thread(d, CELL_ID_FLUID);

    cell_value_exchange_udmi_cellzone(
                                        ts, 
                                        tf, 
                                        &g_node_start_idx_arr_cts[0][0],
                                        &g_node_start_idx_arr_ctf[0][0],
                                        &g_T_cs[0][0],
                                        &g_T_cf[0][0],
                                        &g_matching_sf[0][0],
                                        &g_matching_fs[0][0],
                                        g_N_cs[0],
                                        g_N_cf[0],
                                        &CELL_VAL_FUN_WRAPPER,
                                        0
                                      );
    #endif
}


DEFINE_SOURCE(heat_sink_link_solid_fluid_offgas, c, t, dS, eqn)
{
	real source = 0;

	if( C_UDMI(c, t, 0) > 273.15)
	{
		source = alpha_sf * As_sf * (C_UDMI(c, t, 0) - C_T(c, t));
		dS[eqn] = -alpha_sf * As_sf;
	}
	else
	{
		source = 0;
		dS[eqn] = 0;
	}
	

	return source;
}
