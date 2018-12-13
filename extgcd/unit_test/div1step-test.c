#include "extgcd.h"

#include "stdio.h"

#include "blas.h"
#include "blas_comm.h"


#define LEN 4


sto_t f[LEN];
sto_t g[LEN];

sto_t r00[LEN];
sto_t r01[LEN];
sto_t r10[LEN];
sto_t r11[LEN];

sto_t cr00[LEN];
sto_t cr01[LEN];
sto_t cr10[LEN];
sto_t cr11[LEN];


//int divsteps( sto_t * r00 , sto_t * r01 , sto_t * r10 , sto_t * r11 ,
//        unsigned n , int delta , const sto_t * f , const sto_t * g );

////////////////////////////////////////////////

//  case 1:
//  d=0, f=[3290], g=[2898]
//  n=1 --> d=1, r00=[1], r01=[0], r10=[1693], r11=[3290]
//
//

void set_case1( unsigned *n , int * delta , sto_t * f , sto_t * g )
{
	f[0] = 3290;
	g[0] = 2898;
	delta[0] = 0;
	n[0] = 1;

	unsigned plen = 1;
	printf("case 1: set:\n");
	printf("n=%d, delta=%d, plen=%d\n", n[0], delta[0], plen );
	gfv_fdump( stdout , f , plen , "f[]: " );
	gfv_fdump( stdout , g , plen , "g[]: " );
}

unsigned check_case1( unsigned delta , const sto_t * r00, const sto_t * r01, const sto_t * r10, const sto_t * r11 )
{
	unsigned mlen = 1;
	printf("case 1: check:\n");
	printf("delta=%d, mlen=%d\n", delta, mlen );
	gfv_fdump( stdout , r00 , mlen , "r00[]: " );
	gfv_fdump( stdout , r01 , mlen , "r01[]: " );
	gfv_fdump( stdout , r10 , mlen , "r10[]: " );
	gfv_fdump( stdout , r11 , mlen , "r11[]: " );

	if( delta != 1) { printf("delta fails.\n"); return -1; }
	if( r00[0] != 1 ) { printf("r00 fails.\n"); return -1; }
	if( r01[0] != 0 ) { printf("r01 fails.\n"); return -1; }
	if( r10[0] != 1693 ) { printf("r10 fails.\n"); return -1; }
	if( r11[0] != 3290 ) { printf("r11 fails.\n"); return -1; }

	return 0;
}

////////////////////////////////////////


//  case 2:
//  delta =1, f = [2466,3302,2193,1475], g=[2898,4092,1754,0]
//  n=1 --> delta = 0, m = [ [0,1] , [2898, 2125] ]
//  n=1 --> delta = 1, x [ [1,0] , [2893, 2898] ] ->  [ [[0],[0,1]] , [[1465,0],[1719,2893]] ]
//

void set_case2( unsigned *n , int * delta , sto_t * f , sto_t * g )
{
	f[0] = 2466; f[1] = 3302; f[2] = 2193; f[3] = 1475;
	g[0] = 2898; g[1] = 4092; g[2] = 1754; g[3] = 0;
	delta[0] = 1;
	n[0] = 2;
	unsigned plen = 4;

	printf("case 2: set:\n");
	printf("n=%d, delta=%d, plen=%d\n", n[0], delta[0], plen );
	gfv_fdump( stdout , f , plen , "f[]: " );
	gfv_fdump( stdout , g , plen , "g[]: " );
}

unsigned check_case2( unsigned delta , const sto_t * r00, const sto_t * r01, const sto_t * r10, const sto_t * r11 )
{
	unsigned mlen = 2;
	printf("case 2: check:\n");
	printf("delta=%d, mlen=%d\n", delta, mlen );
	gfv_fdump( stdout , r00 , mlen , "r00[]: " );
	gfv_fdump( stdout , r01 , mlen , "r01[]: " );
	gfv_fdump( stdout , r10 , mlen , "r10[]: " );
	gfv_fdump( stdout , r11 , mlen , "r11[]: " );

	if( delta != 1) { printf("delta fails.\n"); return -1; }
	if( !((r00[0]==0)&&(r00[1]==0)) ) { printf("r00 fails.\n"); return -1; }
	if( !((r01[0]==0)&&(r01[1]==1)) ) { printf("r01 fails.\n"); return -1; }
	if( !((r10[0]==1465)&&(r10[1]==0)) ) { printf("r10 fails.\n"); return -1; }
	if( !((r11[0]==1719)&&(r11[1]==2893)) ) { printf("r11 fails.\n"); return -1; }

	return 0;
}

////////////////////////////////////////////////////

//  case 3:
//  delta =1, f = [2466,3302,2193,1475], g=[2898,4092,1754,0]
//  n=1 --> delta = 0, m = [ [0,1] , [2898, 2125] ]
//  n=1 --> delta = 1, x [ [1,0] , [2893, 2898] ] ->  [ [[0],[0,1]] , [[1465,0],[1719,2893]] ]
//  n=1 --> delta = 0, x [ [0,1] , [442,  1693] ] ->  [ [[1465,0,0],[1719,2893,0]] , [[1105,0,0],[4164,3843,442]] ]

void set_case3( unsigned *n , int * delta , sto_t * f , sto_t * g )
{
	f[0] = 2466; f[1] = 3302; f[2] = 2193; f[3] = 1475;
	g[0] = 2898; g[1] = 4092; g[2] = 1754; g[3] = 0;
	delta[0] = 1;
	n[0] = 3;
	unsigned plen = 4;

	printf("case 3: set:\n");
	printf("n=%d, delta=%d, plen=%d\n", n[0], delta[0], plen );
	gfv_fdump( stdout , f , plen , "f[]: " );
	gfv_fdump( stdout , g , plen , "g[]: " );
}

unsigned check_case3( unsigned delta , const sto_t * r00, const sto_t * r01, const sto_t * r10, const sto_t * r11 )
{
	unsigned mlen = 3;
	printf("case 3: check:\n");
	printf("delta=%d, mlen=%d\n", delta, mlen );
	gfv_fdump( stdout , r00 , mlen , "r00[]: " );
	gfv_fdump( stdout , r01 , mlen , "r01[]: " );
	gfv_fdump( stdout , r10 , mlen , "r10[]: " );
	gfv_fdump( stdout , r11 , mlen , "r11[]: " );

	if( delta != 0) { printf("delta fails.\n"); return -1; }
	if( !((r00[0]==1465)&&(r00[1]==0)&&(r00[2]==0)) ) { printf("r00 fails.\n"); return -1; }
	if( !((r01[0]==1719)&&(r01[1]==2893)&&(r01[2]==0)) ) { printf("r01 fails.\n"); return -1; }
	if( !((r10[0]==1105)&&(r10[1]==0)&&(r10[2]==0)) ) { printf("r10 fails.\n"); return -1; }
	if( !((r11[0]==4164)&&(r11[1]==3843)&&(r11[2]==442)) ) { printf("r11 fails.\n"); return -1; }

	return 0;
}


/////////////////////////////////////////

int main()
{

	unsigned n;
	int delta0;
	int delta;
	unsigned fail;

	printf("\nTester for divsteps:\n\n");

	set_case1( &n , &delta0 , f , g );
	delta = divsteps( r00 , r01 , r10 , r11 , n , delta0 , f , g , LEN );
	fail = check_case1( delta , r00 , r01 , r10 , r11 );
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


	set_case2( &n , &delta0 , f , g );
	delta = divsteps( r00 , r01 , r10 , r11 , n , delta0 , f , g , LEN );
	fail = check_case2( delta , r00 , r01 , r10 , r11 );
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");

	set_case3( &n , &delta0 , f , g );
	delta = divsteps( r00 , r01 , r10 , r11 , n , delta0 , f , g , LEN );
	fail = check_case3( delta , r00 , r01 , r10 , r11 );
	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


	return 0;
}
