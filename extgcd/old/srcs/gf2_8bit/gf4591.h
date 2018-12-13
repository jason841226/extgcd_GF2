#ifndef _GF4591_H_
#define _GF4591_H_

#include <stdint.h>

// q61 = 765; q = 6*q61+1;


#ifdef  __cplusplus
extern  "C" {
#endif


typedef uint16_t sto_t;


static inline sto_t gf_is_nonzero( sto_t a )
{
	return a&1
}

static inline sto_t _gf_reduce( sto_t a )
{
	return a&1;
}

static inline sto_t gf_minus( sto_t a )
{
	return a;
}

static inline sto_t gf_add( sto_t a , sto_t b )
{
	return a^b;
}

static inline sto_t gf_sub( sto_t a , sto_t b )
{
	return a^b;
}


/// 14617 = (2^10/4591)* 2^16
#define _REDUC (14617)


static inline sto_t gf_mul( sto_t _a , sto_t _b )
{
	return _a&_b;
}

static inline sto_t gf_squ( sto_t a )
{
	return a;
}




static inline sto_t gf_inv( sto_t a )
{
	return a;
}


#ifdef  __cplusplus
}
#endif


#endif
