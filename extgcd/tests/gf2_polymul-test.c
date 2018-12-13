#include <stdio.h>
#include "blas.h"
#include "blas_comm.h"

sto_t case1_a[4]={0b00000110,0b11100000,0b10101101,0b10111100};
sto_t case1_b[4]={0b10001010,0b10101101,0b10101010,0b01101101};
sto_t case1_c[8]={0b00111100,0b00101101,0b01011011,0b01011110,0b00101011,0b10100000,0b01011000,0b00111111};


sto_t case2_a[4]={162,230,145,195};
sto_t case2_b[1]={82};
sto_t case2_c[5]={228,15,68,90,61};

sto_t case3_a[4]={82,252,218,0};
sto_t case3_b[1]={162};
sto_t case3_c[5]={228,91,148,119,0};



sto_t c[8]={0,0,0,0,0,0,0,0};
int main()
{
	gfv_polymul(c,case2_a,4,case2_b,1);
	for(int i=0;i<5;i++)
		// printf("%d\n", c[i]);
		if(c[i]!=case2_c[i])
		{
			printf("gf2v_polymul fail\n");
			break;
		}
	printf("PASS\n");
	return 0;
}