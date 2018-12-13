#ifndef _BLAS_COMM_H_
#define _BLAS_COMM_H_

#include <stdint.h>
#include <stdio.h>

#include "blas_config.h"

#include "gf4591.h"


#ifdef  __cplusplus
extern  "C" {
#endif



//////////////// input/output ////////////////////////////

void gfv_fdump(FILE * fp, const sto_t *v, unsigned n_ele , const char * );

unsigned gfv_is_eq( const sto_t * a , const sto_t * b , unsigned len );



#ifdef  __cplusplus
}
#endif



#endif  /// #define _BLAS_COMM_H_


