#ifndef _BLAS_H_
#define _BLAS_H_


#include "gf4591.h"

#include "blas_config.h"


#ifdef _BLAS_AVX2_

#include "blas_u32.h"
#include "blas_avx2.h"

#define conditional_swap conditional_swap_avx2
#define conditional_mov conditional_mov_avx2

#define gfv_add gfv_add_avx2
#define gfv_add_shiftleft_1unit _gfv_add_shiftleft_1unit
#define gfv_minus gfv_minus_avx2

#define gfv_mul_scalar gfv_mul_scalar_avx2
#define gfv_madd gfv_madd_avx2

#else

#include "blas_u32.h"

#define conditional_swap _conditional_swap
#define conditional_mov _conditional_mov

#define gfv_add _gfv_add
#define gfv_minus _gfv_minus

#define gfv_mul_scalar _gfv_mul_scalar
#define gfv_madd _gfv_madd

#endif



////////////////////////////////////////////////////////////


void gfv_set_zero( sto_t * a, unsigned n_ele );


/////////////////  polynomial multiplication //////////////

/// presuming: size_of(c) == len_a + len_b
void gfv_polymul( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b );

/// presuming: size_of(c) == len_a + len_b  /// XXX: should be size_of(c) == len_a + len_b - skiplen
void gfv_polymul_skipdeg( sto_t * c , unsigned skiplen, const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b );


#endif

