
#include "stdio.h"

#include "blas.h"

#include "blas_comm.h"

#include "extgcd.h"


#define LEN 32

sto_t f[LEN];
sto_t g[LEN];

sto_t s[LEN];
sto_t t[LEN];
sto_t r[LEN];

////////////////////////////////////////////////////////////

//sto_t case1_f[18] = {  1393,   2753, 1213, 1932, 1404,    718,4525, 3553, 1396,    1618, 1682, 2300, 1312,    3526, 1794, 1172, 2446,   0};
//sto_t case1_g[18] = { 2712,4559, 634, 102,     1090,350, 207, 4020,     1716,3847, 3572, 46,     1696, 2194, 1482, 2351,      0, 0 };
sto_t case1_f[18] = { 1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 1, 1, 1,   1, 0 };
sto_t case1_g[18] = { 1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 1, 0, 0,   0, 0 };
unsigned case1_len = 17;
sto_t case1_r[18] = { 1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0 };
sto_t case1_s[18] = { 1, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0 };
//sto_t case1_s[18] = { 400, 256, 1812,      2396,    2358, 828, 186,      4146, 492, 3886, 620,      1171, 1936, 2050, 3830,     0    ,0 ,0 };
sto_t case1_t[18] = { 1, 1, 1, 1,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0, 0, 0,   0, 0 };
//sto_t case1_t[18] = { 988, 70, 425, 2434,      3049,   2360, 3206, 1147,      3010, 3948, 2488, 4364,      379, 1810, 1838, 616,      0, 0 };
int case1_deg = 0;

////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////


unsigned is_eq_vec( const sto_t * a , const sto_t * b , unsigned len )
{
	for(unsigned i=0;i<len;i++) { if(a[i]!=b[i]) return 0; }
	return 1;
}



void set_case( sto_t * f , sto_t * g , unsigned * len_f , const sto_t * cf , const sto_t * cg , unsigned clen )
{
	printf("set case:\n");
	len_f[0] = clen;
	for(unsigned i=0;i<clen;i++) f[i] = cf[i];
	for(unsigned i=0;i<clen;i++) g[i] = cg[i];
	gfv_fdump( stdout , f , len_f[0] , "f[]: " );
	gfv_fdump( stdout , g , len_f[0] , "g[]: " );
}

int check_case( const sto_t * r, int deg , const sto_t * s , const sto_t * t ,
	const sto_t * cr , int cdeg , const sto_t * cs , const sto_t *ct , unsigned len_f )
{
	printf("check case:\n");

	int c=0;
	if( 0 > deg ) { printf("deg: %d. fail.\n", deg ); c=-1; }
	if( !(deg==cdeg) ) { printf("deg fail.%d:%d.\n", deg , cdeg ); c=-1; }

	unsigned rlen = cdeg +1;
	if( !is_eq_vec(r,cr,rlen) ) { printf("gcd[] fail.\n"); c=-1; }
	gfv_fdump(stdout,r,rlen,"r[]: "); gfv_fdump(stdout,cr,rlen,"cr]: ");
	unsigned clen = len_f;
	if( !is_eq_vec(s,cs,clen) ) { printf("s[] fail.\n"); c=-1; }
	gfv_fdump(stdout,s,clen,"s[]: "); gfv_fdump(stdout,cs,clen,"cs]: ");
	if( !is_eq_vec(t,ct,clen) ) { printf("t[] fail.\n"); c=-1; }
	gfv_fdump(stdout,t,clen,"t[]: "); gfv_fdump(stdout,ct,clen,"ct]: ");

	return c;
}


////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////


int main()
{

	unsigned len;
	int16_t deg;
	int err;

	printf("Testing ext_gcd().\n");



	set_case( f , g , &len , case1_f , case1_g , case1_len );
	deg = ext_gcd( r , s , t , f , g , len );
	err = check_case( r , deg , s , t , case1_r , case1_deg , case1_s , case1_t , case1_len );
	printf("%s.\n\n", (0==err)?"PASS":"FAIL" );

	return 0;
}
