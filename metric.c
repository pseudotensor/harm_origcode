
#include "decs.h"

/* calculate Kerr metric related quantities,
   in Kerr-Schild coordinates */

/* 
   this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   where i = {1,2,3,4} corresponds to {t,r,theta,phi}
*/

/* NOTE: parameter hides global variable */
void conn_func(double r, double th, double conn[][NDIM][NDIM])
{
	int i,j,k,l ;
	double sth,s2,s4,cth,c2,c4,c6,a2,a4,a6,r2,r3,r6,rho2,rho4 ;
	double tmp[NDIM][NDIM][NDIM] ;

	sth = fabs(sin(th)) ;
	// if(sth < SMALL) sth = SMALL ;
	s2 = sth*sth ;
	s4 = s2*s2 ;
	
	cth = cos(th) ;
	c2 = cth*cth ;
	c4 = c2*c2 ;
	c6 = c4*c2 ;

	a2 = a*a ;
	a4 = a2*a2 ;
	a6 = a4*a2 ;

	r2 = r*r ;
	r3 = r2*r ;
	r6 = r3*r3 ;

	rho2 = r2 + a2*c2 ;
	rho4 = rho2*rho2 ;

	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++) conn[i][j][k] = 0. ;

	/* first calculate out g_{ij,k} */
	conn[TT][TT][RR] = -2.*(r2 - a2*c2)/rho4 ;
	conn[TT][TT][TH] = 4.*a2*r*cth*sth/rho4 ;

	conn[TT][RR][RR] = -2.*(r2 - a2*c2)/rho4 ;
	conn[TT][RR][TH] = 4.*a2*r*cth*sth/rho4 ;

	conn[TT][PH][RR] = 2.*a*(r2 - a2*c2)*s2/rho4 ;
	conn[TT][PH][TH] = -4.*a*r*(a2 + r2)*sth*cth/rho4 ;

	conn[RR][TT][RR] = -2.*(r2 - a2*c2)/rho4 ;
	conn[RR][TT][TH] = 4.*a2*r*cth*sth/rho4 ;

	conn[RR][RR][RR] = -2.*(r2 - a2*c2)/rho4 ;
	conn[RR][RR][TH] = 4.*a2*r*cth*sth/rho4 ;

	conn[RR][PH][RR] = -2.*a*(-r2 + a2*c2)*s2/rho4 ;
	conn[RR][PH][TH] = 2.*a*cth*sth*(-1. - 2.*r/rho2 - 2.*a2*r*s2/rho4) ;

	conn[TH][TH][RR] = 2.*r ;
	conn[TH][TH][TH] = -2.*a2*cth*sth ;

	conn[PH][TT][RR] = 2.*a*(r2 - a2*c2)*s2/rho4 ;
	conn[PH][TT][TH] = -4.*a*r*(a2 + r2)*sth*cth/rho4 ;

	conn[PH][RR][RR] = 2.*a*(r2 - a2*c2)*s2/rho4 ;
	conn[PH][RR][TH] = 2.*a*cth*sth*(-1. - 2.*r/rho2 - 2.*a2*r*s2/rho4) ;

	conn[PH][PH][RR] = 2.*s2*(r + (-a2*r2 + a4*c2)*s2/rho4) ;
	conn[PH][PH][TH] = 2.*cth*sth*(r6 + a6*c6 +
			a2*r3*(4. + r)*s2 +
			2.*a4*r*s4 + c4*(3.*a4*r2 + a6*s2) +
			a2*r*c2*(3.*r3 + 2.*a2*(2. + r)*s2))/rho4 ;

	/* now rearrange to find \Gamma_{ijk} */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++) 
		tmp[i][j][k] = 0.5*(conn[j][i][k] + conn[k][i][j] - conn[k][j][i]) ;

	/* now raise index */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)  {
		conn[i][j][k] = 0. ;
		for(l=0;l<NDIM;l++) conn[i][j][k] += Gcon[i][l]*tmp[l][j][k] ;
	}
}

double gdet_func(double r, double th)
{
	double a2,r2,c2,sth ;

	a2 = a*a ;
	r2 = r*r ;
	sth = fabs(sin(th)) ;
	// if(sth < SMALL) sth = SMALL ;
	c2 = cos(th) ;
	c2 *= c2 ;
	return( 
		fabs( (r2 + a2*c2)*sth )
	) ;
}

void gcov_func(double r, double th, double gcov[][NDIM])
{
	int j,k ;
	double sth,cth,s2,rho2 ;

	DLOOP gcov[j][k] = 0. ;

	cth = cos(th) ;
	sth = fabs(sin(th)) ;
	if(sth < SMALL) sth = SMALL ;
	s2 = sth*sth ;
	rho2 = r*r + a*a*cth*cth ;
	
	gcov[TT][TT] = -1. + 2.*r/rho2 ;
	gcov[TT][RR] = 2.*r/rho2 ;
	gcov[TT][PH] = -2.*a*r*s2/rho2 ;

	gcov[RR][TT] = gcov[TT][RR] ;
	gcov[RR][RR] = 1. + 2.*r/rho2 ;
	gcov[RR][PH] = -a*s2*(1. + 2.*r/rho2) ;

	gcov[TH][TH] = rho2 ;

	gcov[PH][TT] = gcov[TT][PH] ;
	gcov[PH][RR] = gcov[RR][PH] ;
	gcov[PH][PH] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)) ;

}

void gcon_func(double r, double th, double gcon[][NDIM])
{
	int j,k ;
	double sth,cth,s2,rho2 ;

	DLOOP gcon[j][k] = 0. ;

	cth = cos(th) ;
	sth = fabs(sin(th)) ;
	if(sth < SMALL) sth = SMALL ;
	s2 = sth*sth ;
	rho2 = r*r + a*a*cth*cth ;

	gcon[TT][TT] = -1. - 2.*r/rho2 ;
	gcon[TT][RR] = 2.*r/rho2 ;

	gcon[RR][TT] = gcon[TT][RR] ;
	gcon[RR][RR] = 1. - 2.*r/rho2 + a*a*s2/rho2 ;
	gcon[RR][PH] = a/rho2 ;

	gcon[TH][TH] = 1./rho2 ;

	gcon[PH][RR] = gcon[RR][PH] ;
	gcon[PH][PH] = 1./(s2*rho2) ;

}

double delta(int i, int j)
{
	if(i == j) return(1.) ;
	else return(0.) ;
}

/* Minkowski metric; signature +2 */
double mink(int i, int j)
{
	if(i == j) {
		if(i == 0) return(-1.) ;
		else return(1.) ;
	}
	else return(0.) ;
}

/* Boyer-Lindquist ("bl") metric functions */

void blgset(int i, int j)
{
	double bl_gdet_func(double r, double th) ;
	void bl_gcov_func(double r, double th, double gcov[][NDIM]) ;
	void bl_gcon_func(double r, double th, double gcon[][NDIM]) ;
	static double gtmp,gcovtmp[NDIM][NDIM],gcontmp[NDIM][NDIM] ;

	rcurr = Rin + (i + 0.5)*dr ;
	hcurr = 0. + (j + 0.5)*dh ;
	icurr = i ;
	jcurr = j ;
	gtmp = bl_gdet_func(rcurr,hcurr) ;
	bl_gcov_func(rcurr,hcurr,gcovtmp) ;
	bl_gcon_func(rcurr,hcurr,gcontmp) ;

	g = gtmp ;
	Gcov = gcovtmp ;
	Gcon = gcontmp ;
}

double bl_gdet_func(double r, double th)
{
	double a2,r2 ;

	a2 = a*a ;
	r2 = r*r ;
	return( 
		r*r*fabs(sin(th))*(1. + 0.5*(a2/r2)*(1. + cos(2.*th)))
	) ;
}

void bl_gcov_func(double r, double th, double gcov[][NDIM])
{
	int j,k ;
	double sth,cth,s2,a2,r2,DD,mu ;

	DLOOP gcov[j][k] = 0. ;

	sth = fabs(sin(th)) ;
	s2 = sth*sth ;
	cth = cos(th) ;
	a2 = a*a ;
	r2 = r*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;
	
	gcov[TT][TT] = -(1. - 2./(r*mu)) ;
	gcov[TT][PH] = -2.*a*s2/(r*mu) ;
	gcov[PH][TT] = gcov[TT][PH] ;
	gcov[RR][RR] = mu/DD ;
	gcov[TH][TH] = r2*mu ;
	gcov[PH][PH] = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)) ;

}

void bl_gcon_func(double r, double th, double gcon[][NDIM])
{
	int j,k ;
	double sth,cth,a2,r2,r3,DD,mu ;

	DLOOP gcon[j][k] = 0. ;

	sth = fabs(sin(th)) ;
	cth = cos(th) ;
	a2 = a*a ;
	r2 = r*r ;
	r3 = r2*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;

	gcon[TT][TT] = -1. - 2.*(1. + a2/r2)/(r*DD*mu) ;
	gcon[TT][PH] = -2.*a/(r3*DD*mu) ;
	gcon[PH][TT] = gcon[TT][PH] ;
	gcon[RR][RR] = DD/mu ;
	gcon[TH][TH] = 1./(r2*mu) ;
	gcon[PH][PH] = (1. - 2./(r*mu))/(r2*sth*sth*DD) ;

}

