
#include <stdint.h>

#include <stdio.h>

#include "blas_config.h"

#include "gf4591.h"

#include "blas.h"


#include <assert.h>  /// have to be included at last.






//////////////// input/output ////////////////////////////

void gfv_fdump(FILE * fp, const sto_t *v, unsigned n_ele, const char * str ) {
	if( NULL != str ) fprintf(fp, "%s" , str );
	fprintf(fp,"[%2d][",n_ele);
	for(unsigned i=0;i<n_ele;i++) { fprintf(fp,"%4d,",v[i]); if(7==(i%8)) fprintf(fp," ");}
	fprintf(fp,"]\n");
}



unsigned gfv_is_eq( const sto_t * a , const sto_t * b , unsigned len )
{
	sto_t r = 0;
	for(unsigned i=0;i<len;i++) r |= a[i]^b[i];
	return (r==0)?1:0;
}
