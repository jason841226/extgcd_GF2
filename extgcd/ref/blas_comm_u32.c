
#include <stdint.h>

#include <stdio.h>

#include "blas_config.h"

#include "gf4591.h"

#include "blas.h"


#include <assert.h>  /// have to be included at last.




void gfv_set_zero( sto_t * b, unsigned n_ele ) {
	for(unsigned i=0;i<n_ele;i++) b[i]=0;
}


/// school book multiplication.


// void gfv_polymul( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
// {
// 	gfv_set_zero( c , len_a + len_b );
// 	for(unsigned i=0;i<len_b;i++){
// 		gfv_madd( c + i , a , b[i] , len_a );
// 	}
// }

void gfv_polymul( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	sto_t aa, mask, g;
	for(unsigned i=0; i<len_a+len_b; ++i)
		c[i] = 0;
	for(unsigned i_a=0; i_a<len_a; ++i_a)
	{
		for(unsigned i_8=0; i_8<8; ++i_8)
		{
			aa = (a[i_a]>>i_8)&1;
			aa -= 1;
			mask = ~aa;
			for(unsigned i_b=0; i_b<len_b; ++i_b)
			{
				g = b[i_b]&mask;
				c[i_a+i_b] ^= g<<i_8;
				c[i_a+i_b+1] ^= g>>(8-i_8);
			}
		}
	}
}

void gfv_polymul_skipdeg( sto_t * c , unsigned skiplen, const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	gfv_polymul( c , a , len_a , b , len_b );
	unsigned i=0;
	for(;skiplen+i<len_a+len_b;i++) {
		c[i] = c[skiplen+i];
	}
	for(;i<len_a+len_b;i++) c[i] = 0;
}


