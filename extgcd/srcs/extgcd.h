#ifndef _EXTGCD_H_
#define _EXTGCD_H_

#include "stdint.h"

#include "gf4591.h"

#ifdef  __cplusplus
extern  "C" {
#endif





int16_t divsteps( sto_t * r00 , sto_t * r01 , sto_t * r10 , sto_t * r11 ,
	unsigned n , int16_t delta , const sto_t * f , const sto_t * g , unsigned plen );



#define _INIT_DELTA 1

/// reverse extended gcd
/// input:  f, g, and plen
/// output: degree of gcd (return value), gcd, s, and t.  gcd = s*f + t*g
/// presuming:
///     degree of f = plen - 1
///     degree of g = degree of f - _INIT_DELTA
///     storage of gcd, s, t, f, and g  = plen
int16_t rev_extgcd1( sto_t * gcd , sto_t * s , sto_t * t , const sto_t * f , const sto_t * g , unsigned plen );


/// reverse extended gcd
/// input:  f, g, and plen
/// output: degree of gcd (return value), gcd, s, and t.  gcd = s*f + t*g
/// presuming:
///     degree of f = plen - 1
///     degree of g = degree of f - _INIT_DELTA
///     storage of gcd, s, t, f, and g  = plen
int16_t rev_extgcd( sto_t * gcd , sto_t * s , sto_t * t , const sto_t * f , const sto_t * g , unsigned plen );


/// extended gcd
/// input:  f, g, and plen
/// output: degree of gcd (return value), gcd, s, and t.  gcd = s*f + t*g
/// presuming:
///     degree of f = plen - 1
///     degree of g = degree of f - _INIT_DELTA
///     storage of gcd, s, t, f, and g  = plen
int16_t ext_gcd( sto_t * r , sto_t * s , sto_t * t , const sto_t * f , const sto_t * g , unsigned len_f );




#ifdef  __cplusplus
}
#endif

#endif
