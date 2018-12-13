


#include "extgcd.h"
#include "transmat.h"
#include "blas.h"


#define _MAX_LEN 2048

#define SWP(a,b) do { tmp=a; a=b; b=tmp; } while(0)


static void rev_poly( sto_t * rev_a , const sto_t * a, unsigned len, unsigned len_array )
{
	for(unsigned i=0;i<len;i++) rev_a[i] = a[len-1-i];
	for(unsigned i=len;i<len_array;i++) rev_a[i] = 0;
}


int16_t ext_gcd( sto_t * r , sto_t * s , sto_t * t , const sto_t * f , const sto_t * g , unsigned len_f )
{
/// 1. reverse f and g
	sto_t rev_f[_MAX_LEN];
	sto_t rev_g[_MAX_LEN];
	rev_poly( rev_f , f , len_f , len_f );
	rev_poly( rev_g , g , len_f - _INIT_DELTA , len_f );

/// 2. call rev_gcd
	sto_t rev_r[_MAX_LEN];
	sto_t rev_s[_MAX_LEN];
	sto_t rev_t[_MAX_LEN];
	int16_t deg = rev_extgcd( rev_r , rev_s , rev_t , rev_f , rev_g , len_f );

/// 3. reverse r, s, and t.
/// XXX: presuming deg always 0
/// XXX: The following code violates time constancy.
	rev_poly( r , rev_r , 1 , len_f );
	rev_poly( s , rev_s , len_f-1 - _INIT_DELTA , len_f );
	rev_poly( t , rev_t , len_f-1 , len_f );

	return deg;
}


static int16_t div1step_8bit(int16_t delta , sto_t f , sto_t g )
{
	int16_t a00,a01,a10,a11,t00,t01,t10,t11,delta2=0;
	unsigned condition = (delta>0)&(g&1);
	if(condition)
		{a00=0;a01=1;a10=g&1;a11=f&1;}
	else
		{a00=1;a01=0;a10=f&1;a11=g&1;}

	for(unsigned i=1;i<16;++i)
	{
		condition = (delta>0)&((g>>i)&1);
		if(condition)
			{t00=0;t01=1;t10=g&1;t11=f&1;}
		else
			{t00=1;t01=0;t10=f&1;t11=g&1;}


	}
}


//  example:
//  d=0, f=[3290], g=[2890]
//  n=1 --> d=1, r00=[1], r01=[0], r10=[1693], r11=[3290]

///
///  generating 2x2 transition matrix.
///
static int16_t div1step( sto_t * r00 , sto_t * r01 , sto_t * r10 , sto_t * r11 ,
        int16_t delta , const sto_t * f , const sto_t * g )
{
	unsigned condition = (delta>0)&(gf_is_nonzero(g[0]))&1;
	r00[0] = 1;
	r01[0] = 0;
	conditional_swap( condition , r00 , r01 , sizeof(sto_t) );
	sto_t _g0 = gf_minus(g[0]);
	conditional_mov( r10 , condition , g , &_g0 , sizeof(sto_t) );
	sto_t _f0 = gf_minus(f[0]);
	conditional_mov( r11 , condition , &_f0 , f , sizeof(sto_t) );
///
	int16_t _delta = -delta;
	int16_t r_d;
	conditional_mov( &r_d ,condition , &_delta , &delta , sizeof(int16_t) );
	return r_d+1;
}




///
///  Generating transition matirx by multiplying 2x2 matrix.
///    n steps.
///
static int16_t divsteps_1( sto_t * r00 , sto_t * r01 , sto_t * r10 , sto_t * r11 ,
	unsigned n , int16_t delta , const sto_t * f , const sto_t * g , unsigned plen )
{
	sto_t md1_00[1];
	sto_t md1_01[1];
	sto_t md1_10[1];
	sto_t md1_11[1];
	sto_t _m00[_MAX_LEN];
	sto_t _m01[_MAX_LEN];
	sto_t _m10[_MAX_LEN];
	sto_t _m11[_MAX_LEN];

	sto_t _p1_f[_MAX_LEN];
	sto_t _p1_g[_MAX_LEN];
	sto_t _p2_f[_MAX_LEN];
	sto_t _p2_g[_MAX_LEN];
	sto_t * p1_f = &_p1_f[0];
	sto_t * p1_g = &_p1_g[0];
	sto_t * p2_f = &_p2_f[0];
	sto_t * p2_g = &_p2_g[0];

	delta = div1step( md1_00 , md1_01 , md1_10 , md1_11 , delta , f , g );
	deg1transmat_prod_vec( p1_f , p1_g , md1_00 , md1_01 , md1_10 , md1_11 , f , g , (n<plen)?n:plen );

	sto_t * m00 = &_m00[0];
	sto_t * m01 = &_m01[0];
	sto_t * m10 = &_m10[0];
	sto_t * m11 = &_m11[0];
	m00[0] = md1_00[0];
	m01[0] = md1_01[0];
	m10[0] = md1_10[0];
	m11[0] = md1_11[0];
	sto_t * n00 = r00;
	sto_t * n01 = r01;
	sto_t * n10 = r10;
	sto_t * n11 = r11;
	sto_t * tmp;

	for(unsigned i=1;i<n;i++) {
		delta = div1step( md1_00 , md1_01 , md1_10 , md1_11 , delta , p1_f , p1_g );
		deg1transmat_prod_vec( p2_f , p2_g , md1_00 , md1_01 , md1_10 , md1_11 , p1_f , p1_g , ((n-i)<plen)?(n-i):plen );
		deg1transmat_mul( n00 , n01 , n10 , n11 , md1_00 , md1_01 , md1_10 , md1_11 , m00 , m01 , m10 , m11 , i );

		SWP( p1_f , p2_f );
		SWP( p1_g , p2_g );
		SWP( n00 , m00 );
		SWP( n01 , m01 );
		SWP( n10 , m10 );
		SWP( n11 , m11 );
	}

	for(unsigned i=0;i<n;i++) r00[i] = m00[i];
	for(unsigned i=0;i<n;i++) r01[i] = m01[i];
	for(unsigned i=0;i<n;i++) r10[i] = m10[i];
	for(unsigned i=0;i<n;i++) r11[i] = m11[i];
	return delta;
}


///
///  Generating transition matirx with recursive calls
///
int16_t divsteps( sto_t * r00 , sto_t * r01 , sto_t * r10 , sto_t * r11 ,
        unsigned n , int16_t delta , const sto_t * f , const sto_t * g , unsigned plen )
{
	if( 0 == n ) { r00[0]=1; r01[0]=0; r10[0]=0; r11[0]=1; return delta; }
	if( 1 == n ) return div1step( r00 , r01 , r10 , r11 , delta , f , g );

	sto_t m00[_MAX_LEN];
	sto_t m01[_MAX_LEN];
	sto_t m10[_MAX_LEN];
	sto_t m11[_MAX_LEN];

	sto_t n00[_MAX_LEN];
	sto_t n01[_MAX_LEN];
	sto_t n10[_MAX_LEN];
	sto_t n11[_MAX_LEN];

	sto_t p1_f[_MAX_LEN];
	sto_t p1_g[_MAX_LEN];

	unsigned j=n>>1;
	unsigned n_plen = (n<plen)?n:plen;
	delta = divsteps( m00 , m01 , m10 , m11 , j , delta , f , g , n_plen );

	transmat_prod_vec( p1_f , p1_g , m00 , m01 , m10 , m11 , j , f , g , n_plen );

	delta = divsteps( n00 , n01 , n10 , n11 , n-j , delta , p1_f , p1_g , n_plen );
	transmat_mul( r00 , r01 , r10 , r11 , n00 , n01 , n10 , n11 , n-j , m00 , m01 , m10 , m11 , j );

	return delta;
}




/////////////////////////////////////////////////////////////////



/// input:  f, g, and plen
/// output: degree of gcd (return value), gcd, s, and t.  gcd = s*f + t*g
/// presuming:
///     degree of f = plen - 1
///     degree of g = degree of f - _INIT_DELTA
///     storage of gcd, s, t, f, and g  = plen

int16_t rev_extgcd1( sto_t * gcd , sto_t * r00 , sto_t * r01 , const sto_t * f , const sto_t * g , unsigned plen )
{
	if( 0 == plen ) return 0;

	sto_t md1_00[1];
	sto_t md1_01[1];
	sto_t md1_10[1];
	sto_t md1_11[1];
	sto_t _m00[_MAX_LEN];
	sto_t _m01[_MAX_LEN];
	sto_t _m10[_MAX_LEN];
	sto_t _m11[_MAX_LEN];
	sto_t _n00[_MAX_LEN];
	sto_t _n01[_MAX_LEN];
	sto_t _n10[_MAX_LEN];
	sto_t _n11[_MAX_LEN];
	sto_t * m00 = &_m00[0];
	sto_t * m01 = &_m01[0];
	sto_t * m10 = &_m10[0];
	sto_t * m11 = &_m11[0];
	sto_t * n00 = &_n00[0];
	sto_t * n01 = &_n01[0];
	sto_t * n10 = &_n10[0];
	sto_t * n11 = &_n11[0];
	sto_t * tmp;

	sto_t _p1_f[_MAX_LEN];
	sto_t _p1_g[_MAX_LEN];
	sto_t * p1_f = &_p1_f[0];
	sto_t * p1_g = &_p1_g[0];
	sto_t * p2_f = r00;
	sto_t * p2_g = r01;

	int16_t delta = div1step( m00 , m01 , m10 , m11 , _INIT_DELTA , f , g );
	deg1transmat_prod_vec( p1_f , p1_g , m00 , m01 , m10 , m11 , f , g , plen );

	if( 1 < plen ) {
	unsigned n = 2*(plen-1)-_INIT_DELTA;
	for(unsigned i=1;i<n;i++) {
		delta = div1step( md1_00 , md1_01 , md1_10 , md1_11 , delta , p1_f , p1_g );
		deg1transmat_prod_vec( p2_f , p2_g , md1_00 , md1_01 , md1_10 , md1_11 , p1_f , p1_g , plen );
		deg1transmat_mul( n00 , n01 , n10 , n11 , md1_00 , md1_01 , md1_10 , md1_11 , m00 , m01 , m10 , m11 , i );

		SWP( p1_f , p2_f );
		SWP( p1_g , p2_g );
		SWP( n00 , m00 );
		SWP( n01 , m01 );
		SWP( n10 , m10 );
		SWP( n11 , m11 );
	}
	}

	for(unsigned i=0;i<plen;i++) gcd[i] = p1_f[i];
	sto_t inv_hc = gf_inv( gcd[0] );
	gfv_mul_scalar( gcd , inv_hc , plen );
	for(unsigned i=0;i<plen;i++) r00[i] = m00[i];
	gfv_mul_scalar( r00 , inv_hc , plen );
	for(unsigned i=0;i<plen;i++) r01[i] = m01[i];
	gfv_mul_scalar( r01 , inv_hc , plen );
	return delta>>1;
}



inline unsigned _get_msb( unsigned j ) {
	unsigned remain = (j &(j-1));
	while( remain ) {
		j = remain;
		remain = (j&(j-1));
	}
	return j;
}


/// input:  f, g, and plen
/// output: degree of gcd (return value), gcd, s, and t.  gcd = s*f + t*g
/// presuming:
///     degree of f = plen - 1
///     degree of g = degree of f - _INIT_DELTA
///     storage of gcd, s, t, f, and g  = plen

/// XXX:
int16_t rev_extgcd( sto_t * gcd , sto_t * s , sto_t * t , const sto_t * f , const sto_t * g , unsigned plen )
{
	//if( 1 >= plen ) return rev_extgcd( gcd , r00 , r01 , f , g , plen );

	sto_t _m00[_MAX_LEN];
	sto_t _m01[_MAX_LEN];
	sto_t _m10[_MAX_LEN];
	sto_t _m11[_MAX_LEN];
	sto_t * m00 = &(_m00[0]);
	sto_t * m01 = &(_m01[0]);
	sto_t * m10 = &(_m10[0]);
	sto_t * m11 = &(_m11[0]);

	if(1 >= plen ) {
		int16_t r = div1step( s , t , m00 , m01 , 0 , f , g );
		conditional_mov( gcd , s[0] , f , g , sizeof(sto_t) );
		return r;
	}

	sto_t n00[_MAX_LEN];
	sto_t n01[_MAX_LEN];
	sto_t n10[_MAX_LEN];
	sto_t n11[_MAX_LEN];

	sto_t _r00[_MAX_LEN];
	sto_t _r01[_MAX_LEN];
	sto_t _r10[_MAX_LEN];
	sto_t _r11[_MAX_LEN];
	sto_t * r00 = &(_r00[0]);
	sto_t * r01 = &(_r01[0]);
	sto_t * r10 = &(_r10[0]);
	sto_t * r11 = &(_r11[0]);

	sto_t _p1_f[_MAX_LEN];
	sto_t _p1_g[_MAX_LEN];
	sto_t _p2_f[_MAX_LEN];
	sto_t _p2_g[_MAX_LEN];
	sto_t * p1_f = &(_p1_f[0]);
	sto_t * p1_g = &(_p1_g[0]);
	sto_t * p2_f = &(_p2_f[0]);
	sto_t * p2_g = &(_p2_g[0]);
	sto_t * tmp;

	unsigned n=2*(plen-1)-_INIT_DELTA;
printf("plen=%d, n=%d.\n", plen , n );
	unsigned delta = _INIT_DELTA;
	unsigned j= (n&(n-1))? _get_msb(n) : n>>1;

printf("j=%d, n-j=%d.\n" , j , n-j );
	delta = divsteps( m00 , m01 , m10 , m11 , j , delta , f , g , plen );
	transmat_prod_vec( p1_f , p1_g , m00 , m01 , m10 , m11 , j , f , g , plen );
	unsigned m_len = j;
	while( n-j ) {
		n = n-j;
		j = (n&(n-1))?  _get_msb(n) : n ;
printf("j=%d, n-j=%d.\n" , j , n-j );
		delta = divsteps( n00 , n01 , n10 , n11 , j , delta , p1_f , p1_g , plen );
		transmat_prod_vec( p2_f , p2_g , n00 , n01 , n10 , n11 , j , p1_f , p1_g , plen );
		transmat_mul( r00 , r01 , r10 , r11 , n00 , n01 , n10 , n11 , j , m00 , m01 , m10 , m11 , m_len );
		m_len += j;

		SWP( p1_f , p2_f );
		SWP( p1_g , p2_g );
		SWP( r00 , m00 );
		SWP( r01 , m01 );
		SWP( r10 , m10 );
		SWP( r11 , m11 );
	}

	/// output
	for(unsigned i=0;i<plen;i++) gcd[i] = p1_f[i];
	for(unsigned i=0;i<plen;i++) s[i] = m00[i];
	for(unsigned i=0;i<plen;i++) t[i] = m01[i];
	/// normalize
	sto_t inv_hc = gf_inv( gcd[0] );
	gfv_mul_scalar( gcd , inv_hc , plen );
	gfv_mul_scalar( s , inv_hc , plen );
	gfv_mul_scalar( t , inv_hc , plen );

	return delta>>1;
}
