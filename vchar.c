
#include "decs.h"


/*
 * calculate components of magnetosonic velocity 
 * corresponding to primitive variables p 
 *
 * cfg 7-10-01
 * 
 */

#define SAFE 0.2
#define TOL  1.e-8 
#define NTRY 100 

double A[NDIM],B[NDIM],c1,c2,c3,c4,c5,b[NDIM],cms2 ;

double tmpu,tmprho ;

double vrchar(double *pr)
{
	double vchar(double *pr, int js) ;

	return( vchar(pr,RR) ) ;
}

double vhchar(double *pr)
{
	double vchar(double *pr, int js) ;

	return( vchar(pr,TH) ) ;
}

double vchar(double *pr, int js)
{
	extern double chk_disp(double v) ;
	double ucon[NDIM], ecov[NDIM][NDIM], econ[NDIM][NDIM] ;
	double zbrent( double (*func)(double), double v1, double v2, double tol) ;
	void ucon_calc(double *pr, double *ucon) ;
	void b_calc(double *pr, double *ucon, double *b) ;
	void make_co_to_comov(double *ucon, double ecov[][NDIM], double econ[][NDIM]) ;
	void transform(double *vec, double t[][NDIM]) ;
	void coeff_set(double rho, double u) ;
	double v1,v2,vr,err1,err2 ;
	double dv,v0 ;
	double tmp,v ;
	double vr1,vr2 ;
	int j,k,l,m ;

	double r1,h1 ;

	DLOOPA A[j] = 0. ;
	A[js] = 1. ;

	DLOOPA B[j] = 0. ;
	B[TT] = 1. ;

	if(g < SMALL) return(SMALL) ;	/* for regions along poles */

	/* find all components of b, u */
	ucon_calc(pr, ucon) ;
	v0 = ucon[js]/ucon[TT] ;
	b_calc(pr, ucon, b) ;

	/* solve for transform to comoving frame */
	make_co_to_comov(ucon, ecov, econ) ;

	/* transform b to comoving frame */
	transform(b, ecov) ;

	/* find A,B coefficient: K = A + B v */
        transform(A, econ) ;
	transform(B, econ) ;

	/*
	if((jcurr == 0) || (jcurr == NH)) {
		fprintf(stderr,"%d %d\n",icurr,jcurr) ;
		DLOOP fprintf(stderr,"%d %d %g %g\n",j,k,Tetr_con[j][k],Tetr_cov[j][k]) ;
	}
	*/

	/* set up coefficients in dispersion relation */
	coeff_set(pr[0],pr[1]) ;

	/* root finding bit */
	/* bracket */
        dv = SAFE*sqrt(cms2)/sqrt(fabs(Gcov[js][js]) + 1.) ; /* guess; not totally general */
	v1 = -v0 ;
	dv = copysign(dv,v1) ;
	v2 = -v0 + dv ;
        err1 = chk_disp(v1) ;
        err2 = chk_disp(v2) ;

	
        /* based on NR's zbrac */
        for (j=0;j<NTRY;j++) {
                if (err1*err2 < 0.0) break ;
		v1 = v2 ;
		err1 = err2 ;
		v2 += dv ;
		err2 = chk_disp(v2) ;
        }

	if(j == NTRY) {		/* brac failed-- just use speed of light */
		fprintf(stderr,"brac fail %d %d %d\n",icurr,jcurr,pcurr) ;
		PLOOP fprintf(stderr,"%d %g\n",k,pr[k]) ;
		if(js == RR) vr = (-Gcon[TT][RR] + sqrt(Gcon[TT][RR]*Gcon[TT][RR] 
					- Gcon[RR][RR]*Gcon[TT][TT]))/Gcon[TT][TT] ;
		else vr = (-Gcon[TT][TH] + sqrt(Gcon[TT][TH]*Gcon[TT][TH] 
					- Gcon[TH][TH]*Gcon[TT][TT]))/Gcon[TT][TT] ;
	}
	else {			/* brac succeeded-- find root using Brent's method */
		vr =  zbrent(chk_disp, v1, v2, TOL)  ;
	}

	return(fabs(vr)) ;
}

/* make the matrix that transforms from the coordinate
   basis to the comoving orthonormal tetrad */

void make_co_to_comov(double *ucon, double ecov[][NDIM], double econ[][NDIM])
{
	static double boost[NDIM][NDIM] ;
	static double tmp[NDIM] ;
	double v2,f ;
	int j,k,l ;
	void transform(double *ucon, double t[][NDIM]) ;
	double mink(int i, int j) ;

	/** transform ucon to tetrad basis **/
	transform(ucon, Tetr_cov) ;

	/** make boost array **/
	/* space-space components */
	SLOOPA tmp[j] = -ucon[j]/ucon[TT] ;
	v2 = tmp[1]*tmp[1] + tmp[2]*tmp[2] + tmp[3]*tmp[3] ;
	f = (ucon[TT] - 1.)/v2 ;
	SLOOP boost[j][k] = (j == k) ? 1. + tmp[j]*tmp[k]*f : tmp[j]*tmp[k]*f ;

	/* space-time components */
	SLOOPA boost[0][j] = -ucon[j] ;
	SLOOPA boost[j][0] = -ucon[j] ;

	/* time-time component */
	boost[0][0] = ucon[TT] ;

	/** multiply boost by tetrad transform **/
	/* make covariant transform */
	DLOOP {
		ecov[j][k] = 0. ;
		for(l=0;l<NDIM;l++) ecov[j][k] += boost[j][l]*Tetr_cov[l][k] ;
	}
	/* make contravariant transform */
	DLOOP {	 /* re-use boost array as tmp array */
		boost[j][k] = 0. ;
		for(l=0;l<NDIM;l++) boost[j][k] += ecov[j][l]*Gcon[l][k] ;
	}
	DLOOP {
		econ[j][k] = 0. ;
		for(l=0;l<NDIM;l++) econ[j][k] += mink(j,l)*boost[l][k] ;
	}

	return ;
}

void transform(double *v, double t[][NDIM]) 
{
	double tmp[NDIM] ;
	int i,j ;

	for(i=0;i<NDIM;i++) tmp[i] = v[i] ;
	for(i=0;i<NDIM;i++) {
		v[i] = 0. ;
		for(j=0;j<NDIM;j++) v[i] += t[i][j]*tmp[j] ;
	}

	return ;
}


/* evaluate whether dispersion relation is satisfied for 
   covariant wavevector (v,1,0,0) or (v,0,1,0) in coordinate basis */
double chk_disp(double v) 
{
	double w,kx,ky,kz,k2,kdB,kdB2,disp,err ;

	w = A[0] + v*B[0] ;
	kx = A[1] + v*B[1] ;
	ky = A[2] + v*B[2] ;
	kz = A[3] + v*B[3] ;

	k2 = kx*kx + ky*ky + kz*kz ;
	kdB = kx*b[1] + ky*b[2] + kz*b[3] ;
	kdB2 = kdB*kdB ;

	disp = c1*k2 + c2*kdB2 + sqrt(fabs(c3*k2*k2 + c4*k2*kdB2 + c5*kdB2*kdB2)) ;

	err = w*w - disp ;
 
	return(err) ;
}

/* set coefficients for dispersion relation so that they do
   not need to be re-evaluated while v is being solved for */
void coeff_set(double rho, double u)
{
	double B2,EF,EE,va2,cs2 ;

	B2 = b[1]*b[1] + b[2]*b[2] + b[3]*b[3] ;

	EF = rho + gam*u ;
	EE = B2 + EF ;

	va2 = B2/EE ;
	cs2 = gam*(gam - 1.)*u/EF ;

	tmprho = rho ;
	tmpu = u ;

	cms2 = cs2 + va2 - cs2*va2 ;
	if(cms2 < 0.) {
		fprintf(stderr,"crap.\n") ;
		fprintf(stderr,"%g %g %g %g\n",cs2,EE,EF,va2) ;
		exit(1) ;
	}

	/* here are the DR coefficients; see srspeed.ma for
	   an explicit calculation */
	c1 = 0.5*cms2 ;
	c2 = 0.5*cs2/EE ;
	c3 = 0.25*cms2*cms2 ;
	c4 = (0.5*cs2*cms2 - cs2)/EE ;
	c5 = 0.25*cs2*cs2/(EE*EE) ;

}
