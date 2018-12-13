
#include "stdio.h"

#include "blas.h"

#include "blas_comm.h"

#include "extgcd.h"


#define LEN 16

sto_t f[LEN];
sto_t g[LEN];

sto_t s[LEN];
sto_t t[LEN];
sto_t r[LEN];

////////////////////////////////////////////////////////////

sto_t case1_f[4] = {2184,1395,531,0};
sto_t case1_g[4] = {1980,2720,2730,0};
unsigned case1_len = 3;
sto_t case1_r[4] = {1,0,0,0};
sto_t case1_s[4] = {3606,1360,0,0};
sto_t case1_t[4] = {530,541,0,0};
int case1_deg = 0;

////////////////////////////////////////////////////////////

sto_t case2_f[4] = {2562,1556,267,2213};
sto_t case2_g[4] = {3268,408,3397,4272};
unsigned case2_len = 4;
sto_t case2_r[4] = {1,923,0,0};
sto_t case2_s[4] = {0,0,2535,3130};
sto_t case2_t[4] = {0,0,2067,4527};
int case2_deg = 2;

////////////////////////////////////////////////////////////

unsigned is_eq_vec( const sto_t * a , const sto_t * b , unsigned len )
{
	for(unsigned i=0;i<len;i++) { if(a[i]!=b[i]) return 0; }
	return 1;
}

void set_case( sto_t * f , sto_t * g , unsigned * len , const sto_t * cf, const sto_t * cg , unsigned clen )
{
	printf("set case:\n");
	len[0] = clen;
	for(unsigned i=0;i<clen;i++) f[i] = cf[i];
	for(unsigned i=0;i<clen;i++) g[i] = cg[i];
	gfv_fdump( stdout , f , len[0] , "f[]: " );
	gfv_fdump( stdout , g , len[0] , "g[]: " );
}

int check_case( const sto_t * r, int deg , const sto_t * s , const sto_t * t ,
	const sto_t *cr, int cdeg , const sto_t * cs , const sto_t * ct , unsigned clen )
{
	printf("check case:\n");

	int c=0;
	if( 0 > deg ) { printf("deg: %d. fail.\n", deg ); c=-1; }
	if( !(deg==cdeg) ) { printf("deg fail.%d:%d.\n", deg , cdeg ); c=-1; }
	if( !is_eq_vec(r,cr,clen) ) { printf("gcd[] fail.\n"); gfv_fdump(stdout,r,clen,"r[]: "); gfv_fdump(stdout,cr,clen,"cr]: "); c=-1; }
	if( !is_eq_vec(s,cs,clen) ) { printf("s[] fail.\n"); gfv_fdump(stdout,s,clen,"s[]: "); gfv_fdump(stdout,cs,clen,"cs]: "); c=-1; }
	if( !is_eq_vec(t,ct,clen) ) { printf("t[] fail.\n"); gfv_fdump(stdout,t,clen,"t[]: "); gfv_fdump(stdout,ct,clen,"ct]: "); c=-1; }

	return c;
}

////////////////////////////////////////////////////////////


int main()
{

	unsigned len;
	int16_t deg;
	int err;

	printf("Testing rev_extgcd1().\n");

	set_case( f , g , &len , case1_f , case1_g , case1_len );
	deg = rev_extgcd1( r , s , t , f , g , len );
	err = check_case( r , deg , s , t , case1_r , case1_deg , case1_s , case1_t , len );
	printf("%s.\n\n", (0==err)?"PASS":"FAIL" );

	set_case( f , g , &len , case2_f , case2_g , case2_len );
	deg = rev_extgcd1( r , s , t , f , g , len );
	err = check_case( r , deg , s , t , case2_r , case2_deg , case2_s , case2_t , len );
	printf("%s.\n\n", (0==err)?"PASS":"FAIL" );

/////////////////////////////////////////////////

	printf("Testing rev_extgcd().\n");

	set_case( f , g , &len , case1_f , case1_g , case1_len );
	deg = rev_extgcd( r , s , t , f , g , len );
	err = check_case( r , deg , s , t , case1_r , case1_deg , case1_s , case1_t , len );
	printf("%s.\n\n", (0==err)?"PASS":"FAIL" );

	set_case( f , g , &len , case2_f , case2_g , case2_len );
	deg = rev_extgcd( r , s , t , f , g , len );
	err = check_case( r , deg , s , t , case2_r , case2_deg , case2_s , case2_t , len );
	printf("%s.\n\n", (0==err)?"PASS":"FAIL" );

	return 0;
}
