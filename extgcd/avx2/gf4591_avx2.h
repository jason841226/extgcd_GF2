#ifndef _GF4591_AVX2_H_
#define _GF4591_AVX2_H_

#include <stdint.h>

// q61 = 765; q = 6*q61+1;


/// SSE2
#include "emmintrin.h"
/// SSSE3
#include "tmmintrin.h"
/// AVX2
#include "immintrin.h"


#include "gf4591.h"

//typedef uint16_t sto_t;


static const uint16_t __mask_4591[16] __attribute__((aligned(32))) = { 4591,4591, 4591,4591, 4591,4591, 4591,4591, 4591,4591, 4591,4591, 4591,4591, 4591,4591 };
static const __m256i * __ymm_4591 = (__m256i*) &__mask_4591[0];

static inline __m256i _gf_reduce_avx2( __m256i a )
{
	__m256i m4591 = _mm256_load_si256(__ymm_4591);

	__m256i _a = _mm256_sub_epi16( a , m4591 );
	__m256i a_minus = _mm256_cmpgt_epi16( _mm256_setzero_si256() , _a );
	return _mm256_add_epi16( _a , m4591&a_minus );
}

static inline __m256i gf_minus_avx2( __m256i a )
{
	__m256i a_0s = _mm256_cmpeq_epi16( a , _mm256_setzero_si256() );
	return _mm256_sub_epi16( _mm256_andnot_si256( a_0s,_mm256_load_si256(__ymm_4591) ) , a );
}

static inline __m256i gf_add_avx2( __m256i a , __m256i b )
{
	__m256i r = _mm256_add_epi16(a,b);
	return _gf_reduce_avx2( r );
}

static inline __m256i gf_sub_avx2( __m256i a , __m256i b )
{
	__m256i r = _mm256_sub_epi16( a , b );
	__m256i r_minus = _mm256_cmpgt_epi16( _mm256_setzero_si256(), r );
	return _mm256_add_epi16( r_minus&_mm256_load_si256(__ymm_4591) , r );
}


/// 58469 = (2^12/4591)* 2^16
//#define _REDUC (58469)

static const uint16_t __mask_58469[16] __attribute__((aligned(32))) = { 58469,58469, 58469,58469, 58469,58469, 58469,58469, 58469,58469, 58469,58469, 58469,58469, 58469,58469 };
static const __m256i * __ymm_58469 = (__m256i*) &__mask_58469[0];

static const unsigned char __mask_0000ffff[32] __attribute__((aligned(32))) = {0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0, 0xff,0xff,0,0 };
static const __m256i * __ymm_0000ffff = (__m256i*) &__mask_0000ffff[0];


static inline __m256i _estimate_4591x( __m256i c0 , __m256i c1 , __m256i i0000ffff )
{
	/// estimating:  e   s.t. e*4591 ~ c , and e*4591 <= c

	__m256i i58469 = _mm256_load_si256(__ymm_58469);
	//uint32_t d = (c>>10);
	__m256i d = _mm256_srli_epi32(c0,12)^_mm256_andnot_si256(i0000ffff,_mm256_slli_epi32(c1,4));
	//uint32_t e = (d*_REDUC)>>16;
	__m256i e = _mm256_mulhi_epu16(d, i58469 ); /// use unsigned MUL since _REDUC > 32768

	return e;
}


static inline __m256i gf_mul_avx2( __m256i a , __m256i b )
{
	//return ((a*b)%4591);

	__m256i i0000ffff = _mm256_load_si256(__ymm_0000ffff);
	__m256i c0 = _mm256_madd_epi16( a & i0000ffff , b );
	__m256i c1 = _mm256_madd_epi16( _mm256_andnot_si256(i0000ffff,a) , b );

	__m256i e = _estimate_4591x( c0 , c1 , i0000ffff );

	__m256i i4591 = _mm256_load_si256(__ymm_4591);
	/// reduce
	__m256i c = _mm256_mullo_epi16( a , b );
	//uint32_t f = e*4591;
	__m256i f = _mm256_mullo_epi16( e , i4591 );
	//uint32_t _r = c-f;
	__m256i r = _mm256_sub_epi16( c , f );

	return _gf_reduce_avx2(r);
}


static inline __m256i _packlo_16_from_32( __m256i c0 , __m256i c1 , __m256i i0000ffff )
{
	return (c0&i0000ffff)^_mm256_andnot_si256(i0000ffff,_mm256_slli_si256(c1,2));
}


///  (a0,a1), 0, a0, a1
static const uint8_t __shuffle_2x2_a[32] __attribute__((aligned(32))) = { 0,1, 2,3, -1,-1, -1,-1, 0,1, -1,-1, 2,3, -1,-1,    0,1, 2,3, -1,-1, -1,-1, 0,1, -1,-1, 2,3, -1,-1 };
static const __m256i * __ymm_shuffle_2x2_a = (__m256i*) & __shuffle_2x2_a[0];

///  (b1,b0), 0, b0, b1
static const uint8_t __shuffle_2x2_b[32] __attribute__((aligned(32))) = { 2,3, 0,1, -1,-1, -1,-1, 0,1, -1,-1, 2,3, -1,-1,    2,3, 0,1, -1,-1, -1,-1, 0,1, -1,-1, 2,3, -1,-1 };
static const __m256i * __ymm_shuffle_2x2_b = (__m256i*) & __shuffle_2x2_b[0];

static const uint32_t __mask_64[8] __attribute__((aligned(32))) = {-1,-1,0,0,  -1,-1,0,0};
static const __m256i * __ymm_m64 = (__m256i*) &__mask_64[0];


/// XXX: have to check all results value < 4591
static inline __m256i gf_polymul_2x2_avx2( __m256i a , __m256i b )
{
	__m256i aa = _mm256_shuffle_epi8( a , _mm256_load_si256(__ymm_shuffle_2x2_a) );
	__m256i bb = _mm256_shuffle_epi8( b , _mm256_load_si256(__ymm_shuffle_2x2_b) );

	__m256i c1 = _mm256_madd_epi16( aa , bb );
	__m256i c0 = _mm256_srli_si256( c1 , 8 );
	c1 &= _mm256_load_si256( __ymm_m64 );

	__m256i i0000ffff = _mm256_load_si256(__ymm_0000ffff);
	__m256i e = _estimate_4591x( c0 , c1 , i0000ffff );

	__m256i i4591 = _mm256_load_si256(__ymm_4591);
	/// reduce
	__m256i c = _packlo_16_from_32( c0 , c1 , i0000ffff );
	//uint32_t f = e*4591;
	__m256i f = _mm256_mullo_epi16( e , i4591 );
	//uint32_t _r = c-f;
	__m256i r = _mm256_sub_epi16( c , f );
	r = _gf_reduce_avx2(r);
	return _gf_reduce_avx2(r);
}


/// XXX: have to check all results value < 4591
static inline void gf_polymul_16x2_avx2( __m256i * d0 , __m256i * d1 , __m256i a , __m256i b0 , __m256i b1 )
{
	__m256i i0000ffff = _mm256_load_si256(__ymm_0000ffff);
	/// multiplication
	__m256i a0_b0 = _mm256_madd_epi16( a , b0 & i0000ffff );
	__m256i ab1 = _mm256_madd_epi16( a , _mm256_unpacklo_epi16( b1 , b0 )  );
	__m256i a1_b1 = _mm256_madd_epi16( a , _mm256_andnot_si256( i0000ffff, b1 ) );

	/// data arrangement
	__m256i tmp = _mm256_permute2x128_si256( a1_b1 , a1_b1 , 0x08 );
	__m256i ab2_low = _mm256_alignr_epi8( a1_b1 , tmp , 12  );  /// leftshift 4 bytes
	a0_b0 = _mm256_add_epi32( a0_b0 , ab2_low );
	__m256i ab2_high = _mm256_permute4x64_epi64( _mm256_srli_si256(a1_b1,12) , 0xfe );

	/// reduce
	__m256i i4591 = _mm256_load_si256(__ymm_4591);
	__m256i e0 = _estimate_4591x( a0_b0 , ab1 , i0000ffff );
	__m256i c0 = _packlo_16_from_32( a0_b0 , ab1 , i0000ffff );
	__m256i e1 = _estimate_4591x( ab2_high , _mm256_setzero_si256() , i0000ffff );
	__m256i c1 = ab2_high&i0000ffff;
	//uint32_t f = e*4591;
	__m256i f0 = _mm256_mullo_epi16( e0 , i4591 );
	__m256i f1 = _mm256_mullo_epi16( e1 , i4591 );
	//uint32_t _r = c-f;
	__m256i r0 = _mm256_sub_epi16( c0 , f0 );
	__m256i r1 = _mm256_sub_epi16( c1 , f1 );
	r0 = _gf_reduce_avx2(r0);
	r0 = _gf_reduce_avx2(r0);
	r1 = _gf_reduce_avx2(r1);

	_mm256_store_si256( d0 , r0 );
	_mm256_store_si256( d1 , r1 );
}




///  a0, (a0,a1) (a1,a2) a3
static const uint8_t __shuffle_4x4_a0[32] __attribute__((aligned(32))) = { 0,1,-1,-1, 0,1,2,3, 2,3,4,5, 6,7,-1,-1,    0,1,-1,-1, 0,1,2,3, 2,3,4,5, 6,7,-1,-1 };
static const __m256i * __ymm_shuffle_4x4_a0 = (__m256i*) & __shuffle_4x4_a0[0];

///  0,  a2,  a3,  0
static const uint8_t __shuffle_4x4_a1[32] __attribute__((aligned(32))) = { -1,-1,-1,-1, 4,5,-1,-1, 6,7,-1,-1, -1,-1,-1,-1,    -1,-1,-1,-1, 4,5,-1,-1, 6,7,-1,-1, -1,-1,-1,-1 };
static const __m256i * __ymm_shuffle_4x4_a1 = (__m256i*) & __shuffle_4x4_a1[0];

///  (a0,a1), (a0,a1),  (a2,a3),  (a2,a3)
static const uint8_t __shuffle_4x4_a2[32] __attribute__((aligned(32))) = { 0,1,2,3, 0,1,2,3, 4,5,6,7, 4,5,6,7,    0,1,2,3, 0,1,2,3, 4,5,6,7, 4,5,6,7 };
static const __m256i * __ymm_shuffle_4x4_a2 = (__m256i*) & __shuffle_4x4_a2[0];

///  b0,  (b2,b1), (b3,b2), b3
static const uint8_t __shuffle_4x4_b0[32] __attribute__((aligned(32))) = { 0,1,-1,-1, 4,5,2,3, 6,7,4,5, 6,7,-1,-1,    0,1,-1,-1, 4,5,2,3, 6,7,4,5, 6,7,-1,-1, };
static const __m256i * __ymm_shuffle_4x4_b0 = (__m256i*) & __shuffle_4x4_b0[0];

///  0,  b0,  b1,  0
static const uint8_t __shuffle_4x4_b1[32] __attribute__((aligned(32))) = { -1,-1,-1,-1, 0,1,-1,-1, 2,3,-1,-1, -1,-1,-1,-1,    -1,-1,-1,-1, 0,1,-1,-1, 2,3,-1,-1, -1,-1,-1,-1 };
static const __m256i * __ymm_shuffle_4x4_b1 = (__m256i*) & __shuffle_4x4_b1[0];

///  (b1,b0),  (b3,b2), (b3,b2),  (b1,b0)
static const uint8_t __shuffle_4x4_b2[32] __attribute__((aligned(32))) = { 2,3,0,1, 6,7,4,5, 6,7,4,5, 2,3,0,1,    2,3,0,1, 6,7,4,5, 6,7,4,5, 2,3,0,1 };
static const __m256i * __ymm_shuffle_4x4_b2 = (__m256i*) & __shuffle_4x4_b2[0];


static const uint32_t __mask_96[8] __attribute__((aligned(32))) = {-1,-1,-1,0,  -1,-1,-1,0};
static const __m256i * __ymm_m96 = (__m256i*) &__mask_96[0];



/// XXX: have to check all results value < 4591
static inline __m256i gf_polymul_4x4_avx2( __m256i a , __m256i b )
{
	__m256i a0 = _mm256_shuffle_epi8( a , _mm256_load_si256(__ymm_shuffle_4x4_a0) );
	__m256i a1 = _mm256_shuffle_epi8( a , _mm256_load_si256(__ymm_shuffle_4x4_a1) );
	__m256i a2 = _mm256_shuffle_epi8( a , _mm256_load_si256(__ymm_shuffle_4x4_a2) );

	__m256i b0 = _mm256_shuffle_epi8( b , _mm256_load_si256(__ymm_shuffle_4x4_b0) );
	__m256i b1 = _mm256_shuffle_epi8( b , _mm256_load_si256(__ymm_shuffle_4x4_b1) );
	__m256i b2 = _mm256_shuffle_epi8( b , _mm256_load_si256(__ymm_shuffle_4x4_b2) );

	__m256i c0 = _mm256_add_epi32( _mm256_madd_epi16( a0 , b0 ) , _mm256_madd_epi16( a1 , b1 ) );
	__m256i c1x = _mm256_madd_epi16( a2 , b2 );
	__m256i c1l = c1x & _mm256_load_si256( __ymm_m96 );
	__m256i c1 = _mm256_add_epi32( c1l , _mm256_srli_si256( c1x^c1l , 8 ) );

	__m256i i0000ffff = _mm256_load_si256(__ymm_0000ffff);
	__m256i e = _estimate_4591x( c0 , c1 , i0000ffff );

	__m256i i4591 = _mm256_load_si256(__ymm_4591);
	/// reduce
	__m256i c = _packlo_16_from_32( c0 , c1 , i0000ffff );
	//uint32_t f = e*4591;
	__m256i f = _mm256_mullo_epi16( e , i4591 );
	//uint32_t _r = c-f;
	__m256i r = _mm256_sub_epi16( c , f );
	r = _gf_reduce_avx2(r);
	/// 2nd time
	return _gf_reduce_avx2(r);
}


/// XXX: have to check all results value < 4591
static inline void gf_polymul_16x4_avx2( __m256i * d0 , __m256i * d1 , __m256i a , __m256i b0 , __m256i b1, __m256i b2, __m256i b3 )
{

	__m256i i0000ffff = _mm256_load_si256(__ymm_0000ffff);
	/// multiplication
	__m256i a0xb0 = _mm256_madd_epi16( a , b0 & i0000ffff );  /// ab0
	__m256i a0xb1_a1xb0 = _mm256_madd_epi16( a , _mm256_unpacklo_epi16( b1 , b0 )  );  /// ab1
	__m256i a0xb2_a1xb1 = _mm256_madd_epi16( a , _mm256_unpacklo_epi16( b2 , b1 )  );  /// ab2
	__m256i a0xb3_a1xb2 = _mm256_madd_epi16( a , _mm256_unpacklo_epi16( b3 , b2 )  );  /// ab3
	__m256i a1xb3 = _mm256_madd_epi16( a , _mm256_andnot_si256( i0000ffff, b3 ) );   /// ab4

	/// data arrangement
	__m256i tmp1 = _mm256_permute2x128_si256( a0xb2_a1xb1 , a0xb2_a1xb1 , 0x08 );
	__m256i ab2_low = _mm256_alignr_epi8( a0xb2_a1xb1 , tmp1 , 12  );  /// leftshift 4 bytes
	__m256i ab0 = _mm256_add_epi32( a0xb0 , ab2_low );
	__m256i ab2_high = _mm256_permute4x64_epi64( _mm256_srli_si256(a0xb2_a1xb1,12) , 0xfe );

	__m256i tmp2 = _mm256_permute2x128_si256( a0xb3_a1xb2 , a0xb3_a1xb2 , 0x08 );
	__m256i ab3_low = _mm256_alignr_epi8( a0xb3_a1xb2 , tmp2 , 12  );  /// leftshift 4 bytes
	__m256i ab1 = _mm256_add_epi32( a0xb1_a1xb0 , ab3_low );
	__m256i ab3_high = _mm256_permute4x64_epi64( _mm256_srli_si256(a0xb3_a1xb2,12) , 0xfe );

	__m256i tmp3 = _mm256_permute2x128_si256( a1xb3 , a1xb3 , 0x08 );
	__m256i ab4_low = _mm256_alignr_epi8( a1xb3 , tmp3 , 8  );  /// leftshift 8 bytes
	ab0 = _mm256_add_epi32( ab0 , ab4_low );
	ab2_high = _mm256_add_epi32( ab2_high , _mm256_permute4x64_epi64( _mm256_srli_si256(a1xb3,8) , 0xfe ) );

	/// reduce
	__m256i i4591 = _mm256_load_si256(__ymm_4591);
	__m256i e0 = _estimate_4591x( ab0 , ab1 , i0000ffff );
	__m256i c0 = _packlo_16_from_32( ab0 , ab1 , i0000ffff );

	__m256i e1 = _estimate_4591x( ab2_high , ab3_high , i0000ffff );
	__m256i c1 = _packlo_16_from_32( ab2_high , ab3_high , i0000ffff );
	//uint32_t f = e*4591;
	__m256i f0 = _mm256_mullo_epi16( e0 , i4591 );
	__m256i f1 = _mm256_mullo_epi16( e1 , i4591 );
	//uint32_t _r = c-f;
	__m256i r0 = _mm256_sub_epi16( c0 , f0 );
	__m256i r1 = _mm256_sub_epi16( c1 , f1 );
	r0 = _gf_reduce_avx2(r0);
	r1 = _gf_reduce_avx2(r1);
	/// 2nd time
	r0 = _gf_reduce_avx2(r0);
	r1 = _gf_reduce_avx2(r1);

	_mm256_store_si256( d0 , r0 );
	_mm256_store_si256( d1 , r1 );
}







static const uint32_t __mask_h128[8] __attribute__((aligned(32))) = {0,0,0,0, -1,-1,-1,-1 };
static const __m256i * __ymm_h128 = (__m256i*) &__mask_h128[0];

static inline __m256i gf_polymul_8x8_avx2( __m256i a , __m256i b )
{
	__m256i a1_a0 = _mm256_permute4x64_epi64(a,0x98);  ///  imm8 = 10,01,10,00
	__m256i a0_a1 = _mm256_permute4x64_epi64(a,0x89);  ///  imm8 = 10,00,10,01
	__m256i b0_b0 = _mm256_permute4x64_epi64(b,0x88);  ///  imm8 = 10,00,10,00
	__m256i b1_b1 = _mm256_permute4x64_epi64(b,0x99);  ///  imm8 = 10,01,10,01

	__m256i c1_c0 = gf_polymul_4x4_avx2( a1_a0 , b0_b0 );
	__m256i c1_c2 = gf_polymul_4x4_avx2( a0_a1 , b1_b1 );

	__m256i c2_c0 = _mm256_permute2x128_si256( c1_c0 , c1_c2 , 0x20 );  /// imm8 = 00,10|00,00
	__m256i c1_xx = gf_add_avx2( c1_c0 , c1_c2 );
	__m256i c1_00 = c1_xx & _mm256_load_si256( __ymm_h128 );
	__m256i _c1_0 = _mm256_permute4x64_epi64( c1_00 , 0x39 ); /// imm8 = 00,11,10,01
	return gf_add_avx2( _c1_0 , c2_c0 );
}


static const uint8_t __shuffle_broadcast_u16[32] __attribute__((aligned(32))) = { 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1 };
static const __m256i * __ymm_broadcast_u16 = (__m256i*) & __shuffle_broadcast_u16[0];

static inline void gf_polymul_16x8_avx2( __m256i * d0 , __m256i * d1 , __m256i a , __m256i b )
{
	__m256i bb = _mm256_permute4x64_epi64(b,0x44);  ///  imm8 = 01,00,01,00
	__m256i b0 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b1 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b2 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b3 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i c0,c1;
	gf_polymul_16x4_avx2( &c0 , &c1 , a , b0, b1, b2, b3 );

	__m256i b4 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b5 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b6 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	bb = _mm256_srli_si256( bb , 2 );
	__m256i b7 = _mm256_shuffle_epi8( bb , _mm256_load_si256(__ymm_broadcast_u16) );
	__m256i c40, c41;
	gf_polymul_16x4_avx2( &c40 , &c41 , a , b4, b5, b6, b7 );

	/// data arrangement
	__m256i tmp1 = _mm256_permute2x128_si256( c40 , c40 , 0x08 );
	__m256i c40_low = _mm256_alignr_epi8( c40 , tmp1 , 8  );  /// leftshift 8 bytes
	__m256i tmp2 = _mm256_permute2x128_si256( c40 , c40 , 0x83 );
	__m256i c40_high = _mm256_alignr_epi8( c41 , tmp2 , 8 );   /// right shift 8 bytes

	c0 = gf_add_avx2( c0 , c40_low );
	c1 = gf_add_avx2( c1 , c40_high );

	_mm256_store_si256( d0 , c0 );
	_mm256_store_si256( d1 , c1 );
}


static inline void gf_polymul_16x16_avx2( __m256i * r0 , __m256i * r1 , __m256i a , __m256i b )
{
#if 1
	__m256i c0,c1;
	gf_polymul_16x8_avx2( &c0 , &c1 , a , b );
	__m256i b1 = _mm256_permute2x128_si256( b , b , 0x81 );
	__m256i c08,c18;
	gf_polymul_16x8_avx2( &c08 , &c18 , a , b1 );

	/// data arrangement
	__m256i c08_rev = _mm256_permute2x128_si256( c08 , c08 , 0x01 );
	__m256i c18_l128_to_h128 = _mm256_permute2x128_si256( c18 , c18 , 0x08 );
	c18 = c18_l128_to_h128 | _mm256_andnot_si256( _mm256_load_si256( __ymm_h128 ) , c08_rev );

	c0 = gf_add_avx2( c0 , c08_rev & _mm256_load_si256( __ymm_h128 ) );
	c1 = gf_add_avx2( c1 , c18 );

	_mm256_store_si256( r0 , c0 );
	_mm256_store_si256( r1 , c1 );
#else
	__m256i a0 = _mm256_permute2x128_si256( a , a , 0x80 );
	__m256i a1 = _mm256_permute2x128_si256( a , a , 0x82 );

	__m256i b0 = _mm256_permute2x128_si256( b , b , 0x80 );
	__m256i b1 = _mm256_permute2x128_si256( b , b , 0x82 );

	__m256i c0 = gf_polymul_8x8_avx2( a0 , b0 );
	__m256i c2 = gf_polymul_8x8_avx2( a1 , b1 );
	__m256i c1__ = gf_polymul_8x8_avx2( gf_add_avx2(a0,a1) , gf_add_avx2(b0,b1) );
	__m256i c1_ = gf_sub_avx2( c1__ , c0 );
	__m256i c1 = gf_sub_avx2( c1_ , c2 );

	__m256i c1l_to_h = _mm256_permute2x128_si256( c1 , c1 , 0x08 );
	__m256i c1h_to_l = _mm256_permute2x128_si256( c1 , c1 , 0x81 );

	c0 = gf_add_avx2( c0 , c1l_to_h );
	c2 = gf_add_avx2( c2 , c1h_to_l );

	_mm256_store_si256( r0 , c0 );
	_mm256_store_si256( r1 , c2 );
#endif
}



static inline __m256i gf_squ_avx2( __m256i a )
{
	return gf_mul_avx2(a,a);
}




static inline __m256i gf_inv_avx2( __m256i a )
{
	/// 4591 = 0x11ef <-- 1,0001,1110,1111
	/// 4549 : 1,0001,1110,1101
	__m256i x10 = gf_squ_avx2( a );
	__m256i x11 = gf_mul_avx2( x10 , a );
	__m256i x100 = gf_squ_avx2( x10 );
	__m256i x1000 = gf_squ_avx2( x100 );
	__m256i r = x1000;
	for(unsigned i=0;i<2;i++) r = gf_squ_avx2(r);
	r = gf_mul_avx2( r , x11 );
	for(unsigned i=0;i<2;i++) r = gf_squ_avx2(r);
	r = gf_mul_avx2( r , x11 );
	for(unsigned i=0;i<3;i++) r = gf_squ_avx2(r);
	r = gf_mul_avx2( r , x11 );
	for(unsigned i=0;i<2;i++) r = gf_squ_avx2(r);
	r = gf_mul_avx2( r , a );

	return r;
}



#endif
