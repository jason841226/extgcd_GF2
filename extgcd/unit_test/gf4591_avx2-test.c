#include "gf4591_avx2.h"

#include "stdio.h"
#include "stdlib.h"

#define TEST_RUN 5000

int main()
{

	union Data {
		sto_t u16;
		sto_t u16x2[2];
		__m256i u256;
	} a , b , c;
	sto_t ca,cb,cc;
	int fail = 0;

	printf("\nTester random addition:\n\n");
	fail = 0;
	for(unsigned i=1;i<TEST_RUN;i++) {
		a.u16 = rand()%4591;
		b.u16 = rand()%4591;
		c.u256 = gf_add_avx2( a.u256 , b.u256 );

		ca = a.u16;
		cb = b.u16;
		cc = gf_add( ca , cb );

		if( cc != c.u16 ) {
			printf("(x) a=%d + b=%d --> c=%d.\n", a.u16 , b.u16 , c.u16 );
			printf("(o) a=%d + b=%d --> c=%d.\n\n", ca , cb , cc );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");

	printf("\nTester SQU:\n\n");
	fail = 0;
	for(unsigned i=0;i<4591;i++) {
		a.u16 = i;
		c.u256 = gf_squ_avx2( a.u256 );

		ca = i;
		cc = gf_squ( ca );

		if( cc != c.u16 ) {
			printf("(x) a=%d x b=%d --> c=%d.\n", a.u16 , a.u16 , c.u16 );
			printf("(o) a=%d x b=%d --> c=%d.\n\n", ca , ca , cc );
			fail++;
			//break;
		}
	}
	printf("(%d) %s\n\n", fail , (fail)?"FAILD.":"PASSED.");


	printf("\nTester for muls:\n\n");
	fail = 0;
	for(unsigned i=0;i<4591;i++) {
		for(unsigned j=i;j<4591;j++) {
		a.u16 = i;
		a.u16x2[1] = i;
		b.u16 = j;
		b.u16x2[1] = j;
		c.u256 = gf_mul_avx2( a.u256 , b.u256 );

		ca = a.u16;
		cb = b.u16;
		cc = gf_mul( ca , cb );

		if( ( cc != c.u16 )&&( c.u16 != c.u16x2[1] ) ) {
			printf("(x) a=%d x b=%d --> c=%d.\n", a.u16 , b.u16 , c.u16 );
			printf("(x) a=%d x b=%d --> c=%d.\n", a.u16x2[1] , b.u16x2[1] , c.u16x2[1] );
			printf("(o) a=%d x b=%d --> c=%d.\n\n", ca , cb , cc );
			fail = -1;
			break;
		}
		}
		if( 0 != fail ) break;
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");



	printf("\nTester 1== gf(a)*gf_inv(a) for a in (0,4591):\n\n");
	fail = 0;
	for(unsigned i=1;i<4591;i++) {
		a.u16 = i;
		b.u256 = gf_inv_avx2( a.u256 );
		c.u256 = gf_mul_avx2( a.u256 , b.u256 );

		ca = i;
		cb = gf_inv( ca );
		cc = gf_mul( ca , cb );

		if( 1 != c.u16 ) {
			printf("(x) a=%d , a^-1=%d, a*a^-1=%d.\n", a.u16 , b.u16 , c.u16 );
			printf("(o) a=%d , a^-1=%d, a*a^-1=%d.\n\n", ca , cb , cc );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");

	return 0;
}
