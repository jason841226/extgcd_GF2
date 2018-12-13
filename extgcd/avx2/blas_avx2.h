#ifndef _BLAS_AVX2_H_
#define _BLAS_AVX2_H_

#include <stdint.h>
#include <stdio.h>


#include "gf4591.h"


#ifdef  __cplusplus
extern  "C" {
#endif


///  if(condition) SWAP; i.e. a,b = (condition)? b,a : a,b;
void conditional_swap_avx2( uint8_t condition , void * a , void * b , unsigned num_byte );

/// r = (condition)? a : b;
void conditional_mov_avx2( void * _r , uint8_t condition , const void * a , const void * b , unsigned num_byte );


///////////////////


void gfv_add_avx2( sto_t * accu_b, const sto_t * a , unsigned n_ele );

void gfv_minus_avx2( sto_t * _a, const sto_t * a , unsigned n_ele );



////////////////////


void gfv_mul_scalar_avx2( sto_t *a, sto_t b, unsigned n_ele );

void gfv_madd_avx2( sto_t * accu_c, const sto_t * a , sto_t b, unsigned n_ele );



#ifdef  __cplusplus
}
#endif



#endif

