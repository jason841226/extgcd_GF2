#ifndef _GF4591_H_
#define _GF4591_H_

#include <stdint.h>

// q61 = 765; q = 6*q61+1;


#ifdef  __cplusplus
extern  "C" {
#endif


typedef uint8_t sto_t;


//TO DELETE
static inline sto_t gf_is_nonzero( sto_t a )
{
	return a&1;
}

static inline sto_t gf_is_nonzero2(sto_t a,unsigned i)
{
	return (a>>i)&1;
}
// static inline sto_t _gf_reduce( sto_t a )
// {
// 	return a&1;
// }

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

//b is scalar (b=0 or 1) a[i]=a[i]*b
static inline sto_t gf_mul( sto_t _a , sto_t _b )
{
	sto_t mask=~(_b-1);
	return _a&mask;
}

// //return 16-bit result
// static inline uint16_t gf_mul_8bit(sto_t a, sto_t b)
// {
// 	uint16_t aa,mask,_b=b,r=0;
// 	for(unsigned i=0;i<8;++i)
// 	{
// 		aa=(a>>i)&1;
// 		aa-=1;
// 		mask=~aa;
// 		r^=(_b<<i)&mask;
// 	}
// 	return r;
// }

// //return two 8-bit result
// static inline void gf_mul_8bit2(sto_t* c, sto_t a, sto_t b)
// {
// 	uint8_t aa, mask, g;
// 	c[0]=0; c[1]=0;
// 	for(unsigned i=0;i<8;++i)
// 	{
// 		aa=(a>>i)&1;
// 		aa-=1;
// 		mask=~aa;
// 		g=b&mask;
// 		c[0]^=g<<i;
// 		c[1]^=g>>(8-i);
// 	}
// }

//a[i]=a[i]^2
static inline sto_t gf_squ( sto_t a )
{
	return a;
}

//a[i]=a[i]^(-1)
static inline sto_t gf_inv( sto_t a )
{
	return a;
}


#ifdef  __cplusplus
}
#endif


#endif
