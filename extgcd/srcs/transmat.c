
#include <stdint.h>

#include <stdio.h>

#include "blas_config.h"

#include "gf4591.h"

#include "blas.h"

#include "blas_comm.h"

#include <assert.h>  /// have to be included at last.




//  m_l = len(mat[0][0])
//  f_m00 = poly_mul(mat[0][0],vec[0])
//  g_m01 = poly_mul(mat[0][1],vec[1])
//  f = list_add( f_m00 , g_m01 )
//  f_m10 = poly_mul(mat[1][0],vec[0])
//  g_m11 = poly_mul(mat[1][1],vec[1])
//  g = list_add( f_m10 , g_m11 )
//  ## XXX: assert 0 for truncated part.
//  rf = f[m_l-1:]
//  rg = list_pad_zero( g[m_l:] , len(vec[1]) )
//  return  rf , rg

#define MAX_POLYNOMIAL_LEN 1024

static unsigned num_mat_x_vec = 0;

void transmat_prod_vec( sto_t * rf, sto_t * rg,
        const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
        const sto_t * f, const sto_t * g, unsigned len_p )
{
//printf("matrix x vector: %dx%d (%d).\n", len_m , len_p , ++num_mat_x_vec );
	sto_t temp[MAX_POLYNOMIAL_LEN*2];

	gfv_polymul_skipdeg( temp , len_m-1, f , len_p , m00 , len_m );
	for(unsigned i=0;i<(len_p);i++) rf[i] = temp[i];
	gfv_polymul_skipdeg( temp , len_m-1, g , len_p , m01 , len_m );
	gfv_add( rf , temp , len_p );
	gfv_polymul_skipdeg( temp , len_m, f , len_p , m10 , len_m );
	for(unsigned i=0;i<(len_p);i++) rg[i] = temp[i];
	gfv_polymul_skipdeg( temp , len_m, g , len_p , m11 , len_m );
	gfv_add( rg , temp , len_p );
}


//  r_len = m1_len+m2_len
//  r00_u = poly_mul(mat1[0][0],mat2[0][0])
//  r00_d = poly_mul(mat1[0][1],mat2[1][0])
//  r00 = list_add(list_pad_front(r00_u,r_len),list_pad_zero(r00_d,r_len))
//  r01_u = poly_mul(mat1[0][0],mat2[0][1])
//  r01_d = poly_mul(mat1[0][1],mat2[1][1])
//  r01 = list_add(list_pad_front(r01_u,r_len),list_pad_zero(r01_d,r_len))
//  r10_u = poly_mul(mat1[1][0],mat2[0][0])
//  r10_d = poly_mul(mat1[1][1],mat2[1][0])
//  r10 = list_add(list_pad_front(r10_u,r_len),list_pad_zero(r10_d,r_len))

static unsigned num_mat_x_mat = 0;

void transmat_mul( sto_t * r00, sto_t * r01, sto_t * r10, sto_t * r11,
        const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
        const sto_t * n00, const sto_t * n01, const sto_t * n10, const sto_t * n11, unsigned len_n )
{
//printf("matrix x matrix: %dx%d (%d).\n", len_m , len_n , ++num_mat_x_mat );
	sto_t temp[MAX_POLYNOMIAL_LEN*2];
	gfv_polymul( temp , m00 , len_m , n00 , len_n );
	gfv_polymul( r00 , m01 , len_m , n10 , len_n );
	//gfv_add( r00 + 1 , temp , len_m+len_n-1 );
	gfv_add_shiftleft_1unit( r00 , temp , len_m+len_n-1 );
	gfv_polymul( temp , m00 , len_m , n01 , len_n );
	gfv_polymul( r01 , m01 , len_m , n11 , len_n );
	//gfv_add( r01 + 1 , temp , len_m+len_n-1 );
	gfv_add_shiftleft_1unit( r01 , temp , len_m+len_n-1 );

	gfv_polymul( temp , m10 , len_m , n00 , len_n );
	gfv_polymul( r10 , m11 , len_m , n10 , len_n );
	//gfv_add( r10 + 1 , temp , len_m+len_n-1 );
	gfv_add_shiftleft_1unit( r10 , temp , len_m+len_n-1 );
	gfv_polymul( temp , m10 , len_m , n01 , len_n );
	gfv_polymul( r11 , m11 , len_m , n11 , len_n );
	//gfv_add( r11 + 1 , temp , len_m+len_n-1 );
	gfv_add_shiftleft_1unit( r11 , temp , len_m+len_n-1 );
}




////////////////////////////////////////////////////////////////////



void deg1transmat_prod_vec( sto_t * rf, sto_t * rg,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11,
	const sto_t * f, const sto_t * g, unsigned len_p )
{
	conditional_mov( rf , m00[0] , f , g , sizeof(sto_t)*len_p );
	rf[len_p] = 0;

	for(unsigned i=0;i<len_p-1;i++) rg[i] = f[1+i];
	rg[len_p] = 0;
	gfv_mul_scalar( rg , m10[0] , len_p-1 );
	gfv_madd( rg , g+1 , m11[0] , len_p-1 );
	rg[len_p] = 0;
}



void deg1transmat_mul( sto_t * r00, sto_t * r01, sto_t * r10, sto_t * r11,
	const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11,
	const sto_t * n00, const sto_t * n01, const sto_t * n10, const sto_t * n11, unsigned len_n )
{
	/// adjust the alignment for elements of the matrix.
	r10[0] = 0; for(unsigned i=0;i<len_n;i++) r10[i+1] = n00[i];
	r11[len_n] = 0; /// don't wanna do: r11 <- n10
	conditional_mov( r00 , m00[0] , r10 , n10 , sizeof(sto_t)*len_n );
	conditional_mov( r00+len_n , m00[0] , r10+len_n , r11+len_n , sizeof(sto_t) );

	r10[0] = 0; for(unsigned i=0;i<len_n;i++) r10[i+1] = n01[i];
	conditional_mov( r01 , m00[0] , r10 , n11 , sizeof(sto_t)*len_n );
	conditional_mov( r01+len_n , m00[0] , r10+len_n , r11+len_n , sizeof(sto_t) );

	r10[0] = 0; for(unsigned i=0;i<len_n;i++) r10[i+1] = n00[i];
	gfv_mul_scalar( r10+1 , m10[0] , len_n );
	gfv_madd( r10 , n10 , m11[0] , len_n );

	r11[0] = 0; for(unsigned i=0;i<len_n;i++) r11[i+1] = n01[i];
	gfv_mul_scalar( r11+1 , m10[0] , len_n );
	gfv_madd( r11 , n11 , m11[0] , len_n );

}

