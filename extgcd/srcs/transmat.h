#ifndef _TRANSITION_MAT_H_
#define _TRANSITION_MAT_H_

#include <stdint.h>
#include <stdio.h>

#include "blas_config.h"

#include "gf4591.h"  /// for sto_t


#ifdef  __cplusplus
extern  "C" {
#endif




void transmat_prod_vec( sto_t * rf, sto_t * rg,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
	const sto_t * f, const sto_t * g, unsigned len_p );

void transmat_mul( sto_t * r00, sto_t * r01, sto_t * r10, sto_t * r11,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
	const sto_t * n00, const sto_t * n01, const sto_t * n10, const sto_t * n11, unsigned len_n );

////////////////////////////////

void deg1transmat_prod_vec( sto_t * rf, sto_t * rg,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11,
	const sto_t * f, const sto_t * g, unsigned len_p );


void deg1transmat_mul( sto_t * r00, sto_t * r01, sto_t * r10, sto_t * r11,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11,
	const sto_t * n00, const sto_t * n01, const sto_t * n10, const sto_t * n11, unsigned len_n );





#ifdef  __cplusplus
}
#endif



#endif  /// #define _BLAS_COMM_H_


