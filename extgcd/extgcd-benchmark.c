
#include "stdio.h"

#include "blas.h"

#include "blas_comm.h"

#include "extgcd.h"

#include "benchmark.h"

#include "stdlib.h"




#define TEST_RUN 1

#define LEN 1024
#define TEST_LEN 762

sto_t f[LEN];
sto_t g[LEN];

sto_t s[LEN];
sto_t t[LEN];
sto_t r[LEN];

////////////////////////////////////////////////////////////

sto_t alloc_f[LEN];
sto_t alloc_g[LEN];
unsigned alloc_len = LEN;
sto_t alloc_r[LEN];
sto_t alloc_s[LEN];
sto_t alloc_t[LEN];

////////////////////////////////////////////////////////////



unsigned is_eq_vec( const sto_t * a , const sto_t * b , unsigned len )
{
	for(unsigned i=0;i<len;i++) { if(a[i]!=b[i]) return 0; }
	return 1;
}


void set_case( sto_t * f , sto_t * g , unsigned * len )
{
	//printf("set case:\n");
	len[0] = TEST_LEN;
	for(unsigned i=0;i<len[0];i++) f[i] = rand()%4591;
	for(unsigned i=0;i<len[0];i++) g[i] = rand()%4591;
	//gfv_fdump( stdout , f , len[0] , "f[]: " );
	//gfv_fdump( stdout , g , len[0] , "g[]: " );
}


int check_case( const sto_t * r, int deg , const sto_t * s , const sto_t * t )
{
	//printf("check case:\n");
	int c=0;
	//if( 0 > deg ) { printf("deg: %d. fail.\n", deg ); c=-1; }
	//if( !(deg==cdeg) ) { printf("deg fail.%d:%d.\n", deg , cdeg ); c=-1; }
	//if( !is_eq_vec(r,cr,clen) ) { printf("gcd[] fail.\n"); gfv_fdump(stdout,r,clen,"r[]: "); gfv_fdump(stdout,cr,clen,"cr]: "); c=-1; }
	//if( !is_eq_vec(s,cs,clen) ) { printf("s[] fail.\n"); gfv_fdump(stdout,s,clen,"s[]: "); gfv_fdump(stdout,cs,clen,"cs]: "); c=-1; }
	//if( !is_eq_vec(t,ct,clen) ) { printf("t[] fail.\n"); gfv_fdump(stdout,t,clen,"t[]: "); gfv_fdump(stdout,ct,clen,"ct]: "); c=-1; }

	return c;
}

////////////////////////////////////////////////////////////


int main()
{

	unsigned len;
	int16_t deg;
	int err;

	struct benchmark bm;
	char msg[256];

	printf("\nBenchmarking rev_extgcd1():\n");
	bm_init(&bm);
	for(unsigned i=0;i<TEST_RUN;i++) {
		set_case( f , g , &len );
BENCHMARK(bm,{
		deg = rev_extgcd1( r , s , t , f , g , len );
});
		err = check_case( r , deg , s , t );
	}
	printf("#err: %d.\n", err );

	bm_dump(msg,256,&bm);
	printf("result: %s\n\n", msg );

/////////////////////////////////////////////////


	printf("\nBenchmarking rev_extgcd():\n");
	bm_init(&bm);
	for(unsigned i=0;i<TEST_RUN;i++) {
		set_case( f , g , &len );
BENCHMARK(bm,{
		deg = rev_extgcd( r , s , t , f , g , len );
});
		err = check_case( r , deg , s , t );
	}
	printf("#err: %d.\n", err );

	bm_dump(msg,256,&bm);
	printf("result: %s\n\n", msg );

	gfv_polymul_report();

	return 0;
}
