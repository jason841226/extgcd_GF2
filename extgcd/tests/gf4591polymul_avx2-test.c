#include "gf4591_avx2.h"

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "string.h"

#include "blas_comm.h"

#define TEST_RUN 5000

void polymul( sto_t* c , const sto_t * a , unsigned len_a , const sto_t * b , unsigned len_b ) {
	unsigned len_c = len_a + len_b;
	for(unsigned i=0;i<len_c;i++) c[i] = 0;
	for(unsigned i=0;i<len_b;i++) {
		for(unsigned j=0;j<len_a;j++) {
			sto_t ab = gf_mul( a[j] , b[i] );
			c[i+j] = gf_add( c[i+j] , ab );
		}
	}
}



int main()
{

	union Data {
		sto_t u16[32];
		uint32_t u32[16];
		__m256i u256[2];
	} a , b , c, d;
	sto_t ca[32],cb[32],cc[32];
	int fail = 0;

	srand(time(NULL));

/////////////////////////////////////////////////////

	printf("\nTester random polymul 2x2:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		b.u16[0] = rand()%4591;
		b.u16[1] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			b.u16[0] = 4590;
			b.u16[1] = 4590;
		}

		d.u256[0] = gf_polymul_2x2_avx2( a.u256[0] , b.u256[0] );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 2*sizeof(sto_t) );
		polymul( cc , ca , 2 , cb , 2 );

		if( ! gfv_is_eq( cc , d.u16 , 4 ) ) {
			printf("2x2 fail\n");
			gfv_fdump( stdout , a.u16 , 2 , "a : " );
			gfv_fdump( stdout , b.u16 , 2 , "xb: " );
			gfv_fdump( stdout , d.u16 , 4 , "(x) " );
			gfv_fdump( stdout , cc , 4 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


//////////////////////////////////////////////////

	printf("\nTester random polymul 16x2:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		b.u16[0] = rand()%4591;
		b.u16[1] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			b.u16[0] = 4590;
			b.u16[1] = 4590;
		}
		__m256i b0 = _mm256_set1_epi16( b.u16[0] );
		__m256i b1 = _mm256_set1_epi16( b.u16[1] );
		gf_polymul_16x2_avx2( & c.u256[0] , & c.u256[1] , a.u256[0] , b0, b1 );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 2*sizeof(sto_t) );
		polymul( cc , ca , 16 , cb , 2 );

		if( ! gfv_is_eq( cc , c.u16 , 18 ) ) {
			printf("16x2 fail\n");
			gfv_fdump( stdout , a.u16 , 16 , "a : " );
			gfv_fdump( stdout , b.u16 , 2 , "xb: " );
			gfv_fdump( stdout , c.u16 , 18 , "(x) " );
			gfv_fdump( stdout , cc , 18 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");



//////////////////////////////////////////////////

	printf("\nTester random polymul 4x4:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		b.u16[0] = rand()%4591;
		b.u16[1] = rand()%4591;
		b.u16[2] = rand()%4591;
		b.u16[3] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			b.u16[0] = 4590;
			b.u16[1] = 4590;
			b.u16[2] = 4590;
			b.u16[3] = 4590;
		}

		d.u256[0] = gf_polymul_4x4_avx2( a.u256[0] , b.u256[0] );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 4*sizeof(sto_t) );
		polymul( cc , ca , 4 , cb , 4 );

		if( ! gfv_is_eq( cc , d.u16 , 8 ) ) {
			printf("4x4 fail\n");
			gfv_fdump( stdout , a.u16 , 4 , "a : " );
			gfv_fdump( stdout , b.u16 , 4 , "xb: " );
			gfv_fdump( stdout , d.u16 , 8 , "(x) " );
			gfv_fdump( stdout , cc , 8 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


//////////////////////////////////////////////////


	printf("\nTester random polymul 16x4:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		b.u16[0] = rand()%4591;
		b.u16[1] = rand()%4591;
		b.u16[2] = rand()%4591;
		b.u16[3] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			b.u16[0] = 4590;
			b.u16[1] = 4590;
			b.u16[2] = 4590;
			b.u16[3] = 4590;
		}
		__m256i b0 = _mm256_set1_epi16( b.u16[0] );
		__m256i b1 = _mm256_set1_epi16( b.u16[1] );
		__m256i b2 = _mm256_set1_epi16( b.u16[2] );
		__m256i b3 = _mm256_set1_epi16( b.u16[3] );
		gf_polymul_16x4_avx2( & c.u256[0] , & c.u256[1] , a.u256[0] , b0, b1, b2, b3 );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 4*sizeof(sto_t) );
		polymul( cc , ca , 16 , cb , 4 );

		if( ! gfv_is_eq( cc , c.u16 , 20 ) ) {
			printf("16x4 fail\n");
			gfv_fdump( stdout , a.u16 , 16 , "a : " );
			gfv_fdump( stdout , b.u16 , 4 , "xb: " );
			gfv_fdump( stdout , c.u16 , 20 , "(x) " );
			gfv_fdump( stdout , cc , 20 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");



//////////////////////////////////////////////////

	printf("\nTester random polymul 8x8:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		for(unsigned j=0;j<8;j++) b.u16[j] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			for(unsigned j=0;j<8;j++) b.u16[j] = 4590;
		}

		d.u256[0] = gf_polymul_8x8_avx2( a.u256[0] , b.u256[0] );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 8*sizeof(sto_t) );
		polymul( cc , ca , 8 , cb , 8 );

		if( ! gfv_is_eq( cc , d.u16 , 16 ) ) {
			printf("8x8 fail\n");
			gfv_fdump( stdout , a.u16 , 8 , "a : " );
			gfv_fdump( stdout , b.u16 , 8 , "xb: " );
			gfv_fdump( stdout , d.u16 , 16 , "(x) " );
			gfv_fdump( stdout , cc , 16 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


//////////////////////////////////////////////////


	printf("\nTester random polymul 16x8:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		for(unsigned j=0;j<8;j++) b.u16[j] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			for(unsigned j=0;j<8;j++) b.u16[j] = 4590;
		}
		gf_polymul_16x8_avx2( & c.u256[0] , & c.u256[1] , a.u256[0] , b.u256[0] );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 8*sizeof(sto_t) );
		polymul( cc , ca , 16 , cb , 8 );

		if( ! gfv_is_eq( cc , c.u16 , 24 ) ) {
			printf("16x8 fail\n");
			gfv_fdump( stdout , a.u16 , 16 , "a : " );
			gfv_fdump( stdout , b.u16 , 8 , "xb: " );
			gfv_fdump( stdout , c.u16 , 24 , "(x) " );
			gfv_fdump( stdout , cc , 24 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


///////////////////////////////////////////////////

	printf("\nTester random polymul 16x16:\n\n");
	fail = 0;
	for(unsigned i=0;i<TEST_RUN;i++) {
		for(unsigned j=0;j<16;j++) a.u16[j] = rand()%4591;
		for(unsigned j=0;j<16;j++) b.u16[j] = rand()%4591;
		if( 0 == i ) {
			for(unsigned j=0;j<16;j++) a.u16[j] = 4590;
			for(unsigned j=0;j<16;j++) b.u16[j] = 4590;
		}
		gf_polymul_16x16_avx2( & c.u256[0] , & c.u256[1] , a.u256[0] , b.u256[0] );

		memcpy( ca , &a.u16[0] , 16*sizeof(sto_t) );
		memcpy( cb , &b.u16[0] , 16*sizeof(sto_t) );
		polymul( cc , ca , 16 , cb , 16 );

		if( ! gfv_is_eq( cc , c.u16 , 32 ) ) {
			printf("16x16 fail\n");
			gfv_fdump( stdout , a.u16 , 16 , "a : " );
			gfv_fdump( stdout , b.u16 , 16 , "xb: " );
			gfv_fdump( stdout , c.u16 , 32 , "(x) " );
			gfv_fdump( stdout , cc , 32 , "(o) " );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");



	return 0;
}
