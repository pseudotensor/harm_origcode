
/* transforms u^i to kerr-schild from boyer-lindquist */

#include "decs.h"

void bltoks(double *pr,int i,int j)
{
	double ucon[NDIM],tmp[NDIM] ;
	double trans[NDIM][NDIM] ;
	void ucon_calc(double *pr, double *ucon) ;
	void blgset(int i, int j); 
	int k ;

	blgset(i,j) ;
	ucon_calc(pr,ucon) ;

	/* make transform matrix */
	DLOOP trans[j][k] = 0. ;
	DLOOPA trans[j][j] = 1. ;
	trans[0][1] = 2.*rcurr/(rcurr*rcurr - 2.*rcurr + a*a) ;
	trans[3][1] = a/(rcurr*rcurr - 2.*rcurr + a*a) ;

	/* transform ucon; solve for v */
	DLOOPA tmp[j] = 0. ;
	DLOOP tmp[j] += trans[j][k]*ucon[k] ;
	DLOOPA ucon[j] = tmp[j] ;

	pr[UR] = ucon[RR] ;
	pr[UH] = ucon[TH] ;
	pr[UP] = ucon[PH] ;

	/* done! */
}

