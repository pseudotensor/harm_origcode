/* 
	produces an "r8" file.  
*/

#include "decs.h"

void image(FILE *fp)
{
	float iq,liq,a,lmax,lmin,min,max ;
	int i,j,k ;

	/* density mapping is logarithmic, in 255 steps
	   between e^lmax and e^lmin */
	k = RHO ;
	max = -1.e9 ;
	min = 1.e9 ;
	ZLOOP {
		if(p[i][j][k] > max) max = p[i][j][k] ;
		if(p[i][j][k] < min) min = p[i][j][k] ;
	}
	fprintf(stderr,"min,max: %g %g\n",min,max) ;

	lmax = log(max) ;
	lmin = log(min) ;
	a = 256./(lmax - lmin) ;

	for(j=NH-1;j>=0;j--) 
	for(i=0;i<NR;i++) {
		iq = log(p[i][j][k]) ;
		/* liq = a*log(iq) + b ; */
		liq = a*(iq - lmin) ;
		if(liq > 255.) liq = 255. ;
		if(liq < 0.) liq = 0. ;
		fprintf(fp,"%c",(char)((int)liq)) ;
	}
}
