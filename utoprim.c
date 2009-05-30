
#include "decs.h"

double U_target[NPR] ;

#define NMAX	200

/* pr *MUST* contain initial guess */

void Utoprim(double *U, double *pr)
{
	double tolx,tolf ;
	int ntrial,k,test ;
	int mnewt(int ntrail, double *p, int n, double tolx, double tolf) ;
	void primtoU(double *p, double *U) ;
	double pr0[NPR],U0[NPR] ;

	ntrial = 100 ;
	/*
	tolx = 1.e-15;
	tolf = 1.e-15 ;
	*/
	tolx = 1.e-5;
	tolf = 1.e-5;

	if(U[0] < 0.) {
		fprintf(stderr,"uh oh, negative particle density\n") ;
		fail = 1 ;
		return ;
	}

	PLOOP U_target[k] = U[k] ;
	PLOOP pr0[k] = pr[k] ;

	for(k=BR;k<=BP;k++) pr[k] = U[k]/g ;	/* solution is known */
	test = mnewt(ntrial,pr-1, NPR-3, tolx, tolf) ;

	if(test == 0) {
		fprintf(stderr,"mnewt failure: %d %d\n",icurr,jcurr) ;
		primtoU(pr0,U0) ;
		PLOOP fprintf(stderr,"%4d %15.8g %15.8g %15.8g %15.8g\n",
				k,pr[k],pr0[k],U_target[k],U0[k]) ;
		PLOOP pr[k] = pr0[k] ;
		fail = 1 ;
	}

	return ;
}

/* auxiliary function required by mnewt */
void usrfun(double *pr,int n,double *beta,double **alpha)
{
	static double U_curr[NPR] ;
	int k ;
	void primtoU(double *p, double *U), dudp_calc(double *p, double **a);


	primtoU(pr+1,U_curr) ;
	for(k=0;k<NPR-3;k++) beta[k+1] = U_curr[k] - U_target[k] ;
	dudp_calc(pr+1,alpha) ;
}

