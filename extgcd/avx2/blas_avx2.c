
#include <stdint.h>
#include <stdio.h>


#include "gf4591.h"

#include "gf4591_avx2.h"

#include "blas.h"




static inline
__m256i _load_ymm( const sto_t * a , unsigned n_ele ) {
	//n_ele &= 0xf; // make vale in the range [0,15]
	uint16_t temp[16] __attribute__((aligned(32)));
	for(unsigned i=0;i<n_ele;i++) temp[i]=a[i];
	return _mm256_load_si256( (const __m256i*) &temp[0] );
}


static inline
__m256i _load_ymm_padzero( const sto_t * a , unsigned n_ele ) {
	//n_ele &= 0xf; // make vale in the range [0,15]
	uint16_t temp[16] __attribute__((aligned(32))) = {0};
	for(unsigned i=0;i<n_ele;i++) temp[i]=a[i];
	return _mm256_load_si256( (const __m256i*) &temp[0] );
}


static inline
void _store_ymm( sto_t * r , __m256i a , unsigned n_ele ) {
	uint16_t temp[16] __attribute__((aligned(32)));
	_mm256_store_si256( (__m256i*) &temp[0] , a );
	for(unsigned i=0;i<n_ele;i++) r[i]=temp[i];
}

///////////////////


static inline
__m256i _load_ymm_bytes( const void * _a , unsigned n_byte ) {
	const uint8_t * a = (const uint8_t *)_a;
	//n_byte &= 0x1f;
	uint8_t temp[32] __attribute__((aligned(32)));
	for(unsigned i=0;i<n_byte;i++) temp[i]=a[i];
	return _mm256_load_si256( (const __m256i*) &temp[0] );
}

static inline
void _store_ymm_bytes( void * _r , __m256i a , unsigned n_byte ) {
	uint8_t * r = (uint8_t *)_r;
	uint8_t temp[32] __attribute__((aligned(32)));
	_mm256_store_si256( (__m256i*) &temp[0] , a );
	for(unsigned i=0;i<(n_byte);i++) r[i]=temp[i];
}


void conditional_swap_avx2( uint8_t condition , void * _a , void * _b , unsigned num_byte )
{
	__m256i * a = (__m256i *)_a;
	__m256i * b = (__m256i *)_b;
	__m256i mask = _mm256_set1_epi16( 0-((short)condition) );
	unsigned num_ymm = num_byte>>5;
	unsigned rem = num_byte&0x1f;
	for(unsigned i=num_ymm;i!=0;i--){
		__m256i ai = _mm256_loadu_si256(a);
		__m256i bi = _mm256_loadu_si256(b);
		__m256i sum = ai^bi;
		_mm256_storeu_si256( a , sum^(ai&mask)^_mm256_andnot_si256(mask,bi) );
		_mm256_storeu_si256( b , sum^(bi&mask)^_mm256_andnot_si256(mask,ai) );
		a++;
		b++;
	}
	if( rem ) {
		__m256i ai = _load_ymm_bytes(a,rem);
		__m256i bi = _load_ymm_bytes(b,rem);
		__m256i sum = ai^bi;
		_store_ymm_bytes( a , sum^(ai&mask)^_mm256_andnot_si256(mask,bi) , rem );
		_store_ymm_bytes( b , sum^(bi&mask)^_mm256_andnot_si256(mask,ai) , rem );
	}
}


void conditional_mov_avx2( void * _r , uint8_t condition , const void * _a , const void * _b , unsigned num_byte )
{
	const __m256i * a = (const __m256i *)_a;
	const __m256i * b = (const __m256i *)_b;
	__m256i * r = (__m256i *)_r;
	__m256i mask = _mm256_set1_epi16( 0-((short)condition) );
	unsigned num_ymm = num_byte>>5;
	unsigned rem = num_byte&0x1f;
	for(unsigned i=num_ymm;i!=0;i--){
		__m256i ai = _mm256_loadu_si256(a);
		__m256i bi = _mm256_loadu_si256(b);
		__m256i sum = ai^bi;
		_mm256_storeu_si256( r , sum^(bi&mask)^_mm256_andnot_si256(mask,ai) );
		a++;
		b++;
		r++;
	}
	if( rem ) {
		__m256i ai = _load_ymm_bytes(a,rem);
		__m256i bi = _load_ymm_bytes(b,rem);
		__m256i sum = ai^bi;
		_store_ymm_bytes( r , sum^(bi&mask)^_mm256_andnot_si256(mask,ai) , rem );
	}
}


///////////////////


void gfv_add_avx2( sto_t * accu_b, const sto_t * a , unsigned n_ele ) {
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i inp2 = _mm256_loadu_si256( (const __m256i*)accu_b );
		__m256i out = gf_add_avx2( inp1 , inp2 );
		_mm256_storeu_si256( (__m256i*)accu_b , out );
		a += 16;
		accu_b += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm( a , rem );
		__m256i inp2 = _load_ymm( accu_b , rem );
		__m256i out = gf_add_avx2( inp1 , inp2 );
		_store_ymm( accu_b , out , rem );
	}
}


void gfv_minus_avx2( sto_t * _a, const sto_t * a , unsigned n_ele ) {
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i out = gf_minus_avx2( inp1 );
		_mm256_storeu_si256( (__m256i*)_a , out );
		a += 16;
		_a += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm( a , rem );
		__m256i out = gf_minus_avx2( inp1 );
		_store_ymm( _a , out , rem );
	}
}



////////////////////


void gfv_mul_scalar_avx2( sto_t *a, sto_t b, unsigned n_ele ) {
	__m256i ymm_b =_mm256_set1_epi16( b );
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i out = gf_mul_avx2( inp1 , ymm_b );
		_mm256_storeu_si256( (__m256i*)a , out );
		a += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm( a , rem );
		__m256i out = gf_mul_avx2( inp1 , ymm_b );
		_store_ymm( a , out , rem );
	}
}


void gfv_mul_scalar_avx2_2( sto_t * c , const sto_t *a, sto_t b, unsigned n_ele ) {
	__m256i ymm_b =_mm256_set1_epi16( b );
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i out = gf_mul_avx2( inp1 , ymm_b );
		_mm256_storeu_si256( (__m256i*)c , out );
		a += 16;
		c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm( a , rem );
		__m256i out = gf_mul_avx2( inp1 , ymm_b );
		_store_ymm( c , out , rem );
	}
}


void gfv_madd_avx2( sto_t * accu_c, const sto_t * a , sto_t b, unsigned n_ele ) {
	__m256i ymm_b =_mm256_set1_epi16( b );
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i inp2 = _mm256_loadu_si256( (const __m256i*)accu_c );
		__m256i out0 = gf_mul_avx2( inp1 , ymm_b );
		__m256i out1 = gf_add_avx2( out0 , inp2 );
		_mm256_storeu_si256( (__m256i*)accu_c , out1 );
		a += 16;
		accu_c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm( a , rem );
		__m256i inp2 = _load_ymm( accu_c , rem );
		__m256i out0 = gf_mul_avx2( inp1 , ymm_b );
		__m256i out1 = gf_add_avx2( out0 , inp2 );
		_store_ymm( accu_c , out1 , rem );
	}
}





////////////////////  functions for polynomials ///////////////////////////////





void gfv_set_zero( sto_t * b, unsigned n_ele ) {
	//for(unsigned i=0;i<n_ele;i++) b[i]=0;
	unsigned num_ymm = n_ele>>4;
	unsigned rem = n_ele&0xf;
	__m256i zero = _mm256_setzero_si256();
	for(unsigned i=0;i<num_ymm;i++) {
		_mm256_storeu_si256( (__m256i*)b , zero );
		b += 16;
	}
	if( rem ) {
		_store_ymm( b , zero , rem );
	}
}




static void gfv_polymul_2x2( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	__m256i a0 = _load_ymm_padzero( a , len_a );
	__m256i b0 = _load_ymm_padzero( b , len_b );
	__m256i c0 = gf_polymul_2x2_avx2( a0 , b0 );
//	XXX: check which one is faster
//	__m256i c0,c1;
//	__m256i ymm_b =_mm256_set1_epi32( *(const int32_t*)(&b[0]) );
//	gf_polymul_16x2_avx2( &c0 , &c1 , a0 , ymm_b );
	_store_ymm( c , c0 , len_a + len_b );
}

static void gfv_polymadd_2x2( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	__m256i a0 = _load_ymm_padzero( a , len_a );
	__m256i b0 = _load_ymm_padzero( b , len_b );
	__m256i _c = _load_ymm_padzero( c , len_a + len_b );
	__m256i c0 = gf_polymul_2x2_avx2( a0 , b0 );
	c0 = gf_add_avx2( c0 , _c );
	_store_ymm( c , c0 , len_a + len_b );
}

static void gfv_polymul_nx2( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b )
{
#if 0
	gfv_set_zero( c , len_a + 2 );
	for(unsigned i=0;i<len_a;i+=2) {
		unsigned la = (i+2 <= len_a)? 2 : 1;
		gfv_polymadd_2x2( c + i , a + i , la , b , 2 );
	}
#else
	__m256i b0 = _mm256_set1_epi16( b[0] );
	__m256i b1 = _mm256_set1_epi16( b[1] );
	unsigned num_ymm = len_a>>4;
	unsigned rem = len_a&0xf;

	__m256i tail = _mm256_setzero_si256();
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i r0,r1;
		gf_polymul_16x2_avx2( &r0 , &r1 , inp1 , b0, b1 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		_mm256_storeu_si256( (__m256i*) c , r0 );
		a += 16;
		c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm_padzero( a , rem );

		__m256i r0,r1;
		gf_polymul_16x2_avx2( &r0 , &r1 , inp1 , b0, b1 );
		r0 = gf_add_avx2( r0 , tail );
		_store_ymm( c , r0 , rem + 1 );
		tail = _mm256_setzero_si256(); /// here. no extra tail required.
		_store_ymm( c + rem + 1 , tail , 1 );
	}
	else { _store_ymm( c , tail , 2 ); }
#endif
}



static void gfv_polymul_4x4( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	__m256i a0 = _load_ymm_padzero( a , len_a );
	__m256i b0 = _load_ymm_padzero( b , len_b );
	__m256i c0 = gf_polymul_4x4_avx2( a0 , b0 );
	_store_ymm( c , c0 , len_a + len_b );
}


static void gfv_polymul_nx4( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b )
{
	__m256i b0 = _mm256_set1_epi16( b[0] );
	__m256i b1 = _mm256_set1_epi16( b[1] );
	__m256i b2 = _mm256_set1_epi16( b[2] );
	__m256i b3 = _mm256_set1_epi16( b[3] );
	unsigned num_ymm = len_a>>4;
	unsigned rem = len_a&0xf;

	__m256i tail = _mm256_setzero_si256();
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i r0,r1;
		gf_polymul_16x4_avx2( &r0 , &r1 , inp1 , b0, b1, b2, b3 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		_mm256_storeu_si256( (__m256i*) c , r0 );
		a += 16;
		c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm_padzero( a , rem );

		__m256i r0,r1;
		gf_polymul_16x4_avx2( &r0 , &r1 , inp1 , b0, b1, b2, b3 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		if( rem + 4 <= 16 ) _store_ymm( c , r0 , rem + 4 );
		else {
			_mm256_storeu_si256( (__m256i*) c , r0 );
			_store_ymm( c+16 , r1 , rem + 4 - 16 );
		}
	}
	else { _store_ymm( c , tail , 4 ); }
}





static void gfv_polymul_8x8( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	__m256i a0 = _load_ymm_padzero( a , len_a );
	__m256i b0 = _load_ymm_padzero( b , len_b );
	__m256i c0 = gf_polymul_8x8_avx2( a0 , b0 );
	_store_ymm( c , c0 , len_a + len_b );
}


static void gfv_polymul_nx8( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b )
{
	__m256i b0 = _load_ymm_padzero( b , 8 );

	unsigned num_ymm = len_a>>4;
	unsigned rem = len_a&0xf;

	__m256i tail = _mm256_setzero_si256();
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i r0,r1;
		gf_polymul_16x8_avx2( &r0 , &r1 , inp1 , b0 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		_mm256_storeu_si256( (__m256i*) c , r0 );
		a += 16;
		c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm_padzero( a , rem );

		__m256i r0,r1;
		gf_polymul_16x8_avx2( &r0 , &r1 , inp1 , b0 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		if( rem + 8 <= 16 ) _store_ymm( c , r0 , rem + 4 );
		else {
			_mm256_storeu_si256( (__m256i*) c , r0 );
			_store_ymm( c+16 , r1 , rem + 8 - 16 );
		}
	}
	else { _store_ymm( c , tail , 8 ); }
}


/// un-tested code.
static void gfv_polymul_16x16( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	__m256i a0 = _load_ymm_padzero( a , len_a );
	__m256i b0 = _load_ymm_padzero( b , len_b );
	__m256i c0,c1;
	gf_polymul_16x16_avx2( &c0 , &c1 , a0 , b0 );
	_mm256_storeu_si256( (__m256i*)c , c0 );
	_store_ymm( c+8 , c1 , len_a+len_b-16 );
}

static void gfv_polymul_nx16( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b )
{
	__m256i b0 = _mm256_loadu_si256( (const __m256i*) b );
	unsigned num_ymm = len_a>>4;
	unsigned rem = len_a&0xf;

	__m256i tail = _mm256_setzero_si256();
	for(unsigned i=0;i<num_ymm;i++) {
		__m256i inp1 = _mm256_loadu_si256( (const __m256i*)a );
		__m256i r0,r1;
		gf_polymul_16x16_avx2( &r0 , &r1 , inp1 , b0 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		_mm256_storeu_si256( (__m256i*) c , r0 );
		a += 16;
		c += 16;
	}
	if( rem ) {
		__m256i inp1 = _load_ymm_padzero( a , rem );

		__m256i r0,r1;
		gf_polymul_16x16_avx2( &r0 , &r1 , inp1 , b0 );
		r0 = gf_add_avx2( r0 , tail );
		tail = r1;
		_mm256_storeu_si256( (__m256i*) c , r0 );
		_store_ymm( c+16 , tail , rem );
	}
	else { _mm256_storeu_si256( (__m256i*)c , tail ); }
}




#define _PROF_MUL_

#ifdef _PROF_MUL_
static unsigned _x_762 = 0;
static unsigned _x_n[10] = {0};

void gfv_polymul_report()
{
	for(unsigned i=0;i<10;i++) {
		printf("x %d: %d\n", 1<<i , _x_n[i] );
	}
	printf("x 762: %d\n", _x_762 );
}
#endif

void gfv_polymul( sto_t * c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	if( len_a < len_b ) {
		len_b ^= len_a; len_a = len_b^len_a; len_b ^= len_a;
		//b ^= a; a = b^a; b ^= a;
		const sto_t * ptmp = a; a = b; b = ptmp;
	}
#ifdef _PROF_MUL_
	if( (len_b&(len_b-1))  ) {
		if( 762 == len_b ) {
			//printf("%d x 762 : %d.\n", len_a , ++x_762 );
			_x_762++;
		} else {
			printf("unit: %d.\n", len_b );
			exit(-1);
		}
	} else {
		_x_n[ __builtin_ctz(len_b) ] ++;
	}
#endif
	//printf("poly mul: %d x %d.\n", len_a , len_b );
	if( 1 == len_b ) {  gfv_mul_scalar_avx2_2( c, a, b[0], len_a ); c[len_a] = 0; return; }
	if( 2 == len_b ) {
		if( 2 >= len_a ) { gfv_polymul_2x2( c , a , len_a , b , len_b ); return; }
		else { gfv_polymul_nx2( c , a , len_a , b ); return; }
	}
	if( 4 == len_b ) {
		if( 4 >= len_a ) { gfv_polymul_4x4( c , a , len_a , b , len_b ); return; }
		else { gfv_polymul_nx4( c , a , len_a , b ); return; }
	}
	if( 8 == len_b ) {
		if( 8 >= len_a ) { gfv_polymul_8x8( c , a , len_a , b , len_b ); return; }
		else { gfv_polymul_nx8( c , a , len_a , b ); return; }
	}
	if( 16 == len_b ) { gfv_polymul_nx16( c , a , len_a , b ); return; }

	gfv_set_zero( c , len_a + len_b );
	for(unsigned i=0;i<len_b;i++){
		gfv_madd( c + i , a , b[i] , len_a );
	}
}

void gfv_polymul_skipdeg( sto_t * c , unsigned skiplen, const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b )
{
	//printf("skipped(%d) poly mul: %d x %d.\n", skiplen , len_a , len_b );
	gfv_polymul( c , a , len_a , b , len_b );
	unsigned i=0;
	for(;skiplen+i<len_a+len_b;i++) {
		c[i] = c[skiplen+i];
	}
	for(;i<len_a+len_b;i++) c[i] = 0;
}
