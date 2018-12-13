#ifndef _BLAS_U32_H_
#define _BLAS_U32_H_

#include <stdint.h>
#include <stdio.h>


#include "gf4591.h"


#ifdef  __cplusplus
extern  "C" {
#endif


/// if(condition) SWAP; i.e. a,b = (condition)? b,a : a,b;
static inline
void _conditional_swap( uint8_t condition , void * _a , void * _b , unsigned num_byte )
{
	uint8_t * a = (uint8_t *)_a;
	uint8_t * b = (uint8_t *)_b;
	uint8_t mask = 0-condition;
	for(unsigned i=0;i<num_byte;i++) {
		uint8_t sum = a[i]^b[i];
		a[i] = sum^(a[i]&mask)^(b[i]&(~mask));
		b[i] = sum^(b[i]&mask)^(a[i]&(~mask));
	}
}

/// ,or use x86 cmov instruction ?
/// r = (condition)? a: b;
static inline
void _conditional_mov( void * _r , uint8_t condition , const void * _a , const void * _b , unsigned num_byte )
{
	const uint8_t * a = (const uint8_t *)_a;
	const uint8_t * b = (const uint8_t *)_b;
	uint8_t * r = (uint8_t *)_r;
	uint8_t mask = 0-condition;
	for(unsigned i=0;i<num_byte;i++) {
		uint8_t sum = a[i]^b[i];
		r[i] = sum^(b[i]&mask)^(a[i]&(~mask));
	}
}



///////////////////


static inline
void _gfv_add( sto_t * accu_b, const sto_t * a , unsigned n_ele ) {
	for(unsigned i=0;i<n_ele;i++) accu_b[i] = gf_add( accu_b[i] , a[i] );
}

static inline
void _gfv_add_shiftleft_1unit( sto_t * accu_b, const sto_t * a , unsigned n_ele ) {
	//for(unsigned i=0;i<n_ele;i++) accu_b[i+1] = gf_add( accu_b[i+1] , a[i] );
	_gfv_add( accu_b+1 , a , n_ele );
}


static inline
void _gfv_minus( sto_t * _a, const sto_t * a , unsigned n_ele ) {
	for(unsigned i=0;i<n_ele;i++) _a[i] = gf_minus( a[i] );
}



////////////////////


static inline
void _gfv_mul_scalar( sto_t *a, sto_t b, unsigned n_ele ) {
	for(unsigned i=0;i<n_ele;i++) a[i] = gf_mul( a[i] , b );
}

static inline
void _gfv_madd( sto_t * accu_c, const sto_t * a , sto_t b, unsigned n_ele ) {
	for(unsigned i=0;i<n_ele;i++) {
		sto_t tmp = gf_mul( a[i] , b );
		accu_c[i] = gf_add( accu_c[i] , tmp );
	}
}





#ifdef  __cplusplus
}
#endif



#endif

