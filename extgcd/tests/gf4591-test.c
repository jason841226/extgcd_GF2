#include "gf4591.h"

#include "stdio.h"


int main()
{

	printf("\nTester 1== gf(a)*gf_inv(a) for a in (0,4591):\n\n");

	int fail = 0;
	for(unsigned i=1;i<4591;i++) {
		sto_t a = i;
		sto_t a_1 = gf_inv(a);
		sto_t c = gf_mul( a , a_1 );
		if( 1 != c ) {
			printf("a=%d , a^-1=%d, a*a^-1=%d.\n", a , a_1 , c );
			fail = -1;
			break;
		}
	}
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");

	return 0;
}
