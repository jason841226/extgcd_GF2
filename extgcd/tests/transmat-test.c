
#include "stdio.h"

#include "blas.h"
#include "blas_comm.h"
#include "transmat.h"


#define LEN 4

sto_t f[LEN];
sto_t g[LEN];

sto_t m00[LEN];
sto_t m01[LEN];
sto_t m10[LEN];
sto_t m11[LEN];

sto_t rf[LEN*2];
sto_t rg[LEN*2];



//void transmat_prod_vec( sto_t * rf, sto_t * rg,
//        const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
//        const sto_t * f, const sto_t * g, unsigned len_p );
//  case 1:
//  m = [ [0,1] , [2898, 2125] ] , f = [2466,3302,2193,1475], g=[2898,4092,1754,0]
//  rf=[2890,4092,1754,0] , rg = [1698,728,329,0]
//
//

void set_case1( sto_t * m00, sto_t * m01, sto_t * m10, sto_t * m11, unsigned * len_m,
        sto_t * f, sto_t * g, unsigned * len_p )
{
	m00[0] = 0;
	m01[0] = 1;
	m10[0] = 82;
	m11[0] = 162;
	len_m[0] = 1;
	f[0] = 162; f[1] = 230; f[2] = 145; f[3] = 195;
	g[0] =  82; g[1] = 252; g[2] = 218; g[3] = 0;
	len_p[0] = 4;
}

unsigned check_case1( const sto_t * rf, const sto_t * rg )
{
	if( rf[0] != 82) { printf("f[0] fails.\n"); return -1; }
	if( rf[1] != 252) { printf("f[1] fails.\n"); return -1; }
	if( rf[2] != 218) { printf("f[2] fails.\n"); return -1; }
	if( rf[3] != 0)    { printf("f[3] fails.\n"); return -1; }

	if( rg[0] != 84) { printf("g[0] fails.\n"); return -1; }
	if( rg[1] != 208)  { printf("g[1] fails.\n"); return -1; }
	if( rg[2] != 45)  { printf("g[2] fails.\n"); return -1; }
	if( rg[3] != 61)    { printf("g[3] fails.\n"); return -1; }

	return 0;
}



///////////////////////////////////////

sto_t n00[LEN];
sto_t n01[LEN];
sto_t n10[LEN];
sto_t n11[LEN];

sto_t r00[LEN];
sto_t r01[LEN];
sto_t r10[LEN];
sto_t r11[LEN];

//void transmat_mul( sto_t * r00, sto_t * r01, sto_t * r10, sto_t * r11,
//        const sto_t * m00, const sto_t * m01, const sto_t * m10, const sto_t * m11, unsigned len_m,
//        const sto_t * n00, const sto_t * n01, const sto_t * n10, const sto_t * n11, unsigned len_n )
//  case 1:
//  n = [ [1,0] , [2893,2898] ]  x  m = [ [0,1] , [2898, 2125] ]
//  r = [ [[0],[0,1]] , [[1465,0],[1719,2893]] ]

void set_mul_case1( sto_t * m00, sto_t * m01, sto_t * m10, sto_t * m11, unsigned * len_m,
        sto_t * n00, sto_t * n01, sto_t * n10, sto_t * n11, unsigned * len_n )
{
	m00[0] = 0;
	m01[0] = 1;
	m10[0] = 2898;
	m11[0] = 2125;
	len_m[0] = 1;

	n00[0] = 1;
	n01[0] = 0;
	n10[0] = 2893;
	n11[0] = 2898;
	len_n[0] = 1;
}

unsigned check_mul_case1( const sto_t * r00, const sto_t * r01, const sto_t * r10, const sto_t * r11 )
{
	if( r00[0] != 0) { printf("r00[0] fails.\n"); return -1; }
	if( r00[1] != 0) { printf("r00[1] fails.\n"); return -1; }

	if( r01[0] != 0) { printf("r01[0] fails.\n"); return -1; }
	if( r01[1] != 1) { printf("r01[1] fails.\n"); return -1; }

	if( r10[0] != 1465) { printf("r10[0] fails.\n"); return -1; }
	if( r10[1] != 0) { printf("r10[1] fails.\n"); return -1; }

	if( r11[0] != 1719) { printf("r11[0] fails.\n"); return -1; }
	if( r11[1] != 2893) { printf("r11[1] fails.\n"); return -1; }
	return 0;
}


//  case 2:
//  m = [[0,1] , [442,1693] ] x r = [ [[0],[0,1]] , [[1465,0],[1719,2893]] ]
//  n = [ [[1465,0,0],[1719,2893,0]] , [[1105,0,0],[4164,3843,442]] ]


void set_mul_case2( sto_t * m00, sto_t * m01, sto_t * m10, sto_t * m11, unsigned * len_m )
{
	m00[0] = 0;
	m01[0] = 1;
	m10[0] = 442;
	m11[0] = 1693;
	len_m[0] = 1;
}

unsigned check_mul_case2( const sto_t * r00, const sto_t * r01, const sto_t * r10, const sto_t * r11 )
{
	if( r00[0] != 1465) { printf("r00[0] fails.\n"); return -1; }
	if( r00[1] != 0) { printf("r00[1] fails.\n"); return -1; }
	if( r00[2] != 0) { printf("r00[2] fails.\n"); return -1; }

	if( r01[0] != 1719) { printf("r01[0] fails.\n"); return -1; }
	if( r01[1] != 2893) { printf("r01[1] fails.\n"); return -1; }
	if( r01[2] != 0) { printf("r01[2] fails.\n"); return -1; }

	if( r10[0] != 1105) { printf("r10[0] fails.\n"); return -1; }
	if( r10[1] != 0) { printf("r10[1] fails.\n"); return -1; }
	if( r10[2] != 0) { printf("r10[2] fails.\n"); return -1; }

	if( r11[0] != 4164) { printf("r11[0] fails.\n"); return -1; }
	if( r11[1] != 3843) { printf("r11[1] fails.\n"); return -1; }
	if( r11[2] != 442) { printf("r11[2] fails.\n"); return -1; }
	return 0;
}



//////////////////////////////////////////////////////////////////////



int main()
{

	unsigned len_m,len_p;

	set_case1( m00 , m01 , m10 , m11 , &len_m , f , g , &len_p );

	transmat_prod_vec( rf , rg , m00 , m01 , m10 , m11 , len_m , f , g , len_p );

	printf("\nTester for transmat_mul():\n\n");
	printf("test case1:\n");

	gfv_fdump( stdout , m00 , len_m , "m[00]: " );
	gfv_fdump( stdout , m01 , len_m , "m[01]: " );
	gfv_fdump( stdout , m10 , len_m , "m[10]: " );
	gfv_fdump( stdout , m11 , len_m , "m[11]: " );
	gfv_fdump( stdout , f , len_p , "f[]: " );
	gfv_fdump( stdout , g , len_p , "g[]: " );
	gfv_fdump( stdout , rf , sizeof(rf)/sizeof(sto_t) , "rf[]: " );
	gfv_fdump( stdout , rg , sizeof(rg)/sizeof(sto_t) , "rg[]: " );

	unsigned fail = check_case1( rf , rg );

	printf("%s\n", (fail)?"FAILD.":"PASSED.");

	printf("\n\n");

///////////////////////////////////////

	unsigned len_n;

	set_mul_case1( m00 , m01 , m10 , m11 , &len_m , n00 , n01 , n10 , n11 , &len_n );

	transmat_mul( r00 , r01 , r10 , r11 , n00 , n01 , n10 , n11 , len_n , m00 , m01 , m10 , m11 , len_m );

	printf("Tester for transmat_prod_vec():\n\n");
	printf("test mul case1: n[] x m[] :\n");

	gfv_fdump( stdout , m00 , len_m , "m[00]: " );
	gfv_fdump( stdout , m01 , len_m , "m[01]: " );
	gfv_fdump( stdout , m10 , len_m , "m[10]: " );
	gfv_fdump( stdout , m11 , len_m , "m[11]: " );

	gfv_fdump( stdout , n00 , len_n , "n[00]: " );
	gfv_fdump( stdout , n01 , len_n , "n[01]: " );
	gfv_fdump( stdout , n10 , len_n , "n[10]: " );
	gfv_fdump( stdout , n11 , len_n , "n[11]: " );

	gfv_fdump( stdout , r00 , len_m+len_n , "r[00]: " );
	gfv_fdump( stdout , r01 , len_m+len_n , "r[01]: " );
	gfv_fdump( stdout , r10 , len_m+len_n , "r[10]: " );
	gfv_fdump( stdout , r11 , len_m+len_n , "r[11]: " );

	fail = check_mul_case1( r00 , r01 , r10 , r11 );

	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


	set_mul_case2( m00 , m01 , m10 , m11 , &len_m );

	transmat_mul( n00 , n01 , n10 , n11 , m00 , m01 , m10 , m11 , len_m , r00 , r01 , r10 , r11 , 2 );

	printf("test mul cas2: m[] x r[]\n");

	fail = check_mul_case2( n00 , n01 , n10 , n11 );

	printf("%s\n\n", (fail)?"FAILD.":"PASSED.");


	return 0;
}
