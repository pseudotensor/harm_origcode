

/* in "tens" it is assumed that the *first* index is up, and
 * the *second* index is down */


#include "decs.h"

/* calculate fluxes in radial direction */
void primtorflux(double *pr, double *f)
{
	double tens[NDIM][NDIM],b[NDIM],ucon[NDIM] ;
	void stress_calc(double *pr, double tens[][NDIM], double *b, double *ucon) ;
	int k ;

	stress_calc(pr,tens,b,ucon) ;
	
	/* simple fluxes */
	f[RHO] = pr[RHO]*pr[UR] ;
	f[UU]  = tens[RR][TT] + f[RHO] ;
	f[UR]  = tens[RR][RR] ;
	f[UH]  = tens[RR][TH] ;
	f[UP]  = tens[RR][PH] ;
	f[BR]  = 0. ;
	f[BH]  = b[TH]*pr[UR] - b[RR]*pr[UH] ;
	f[BP]  = b[PH]*pr[UR] - b[RR]*pr[UP] ;

	PLOOP f[k] *= g ;
}

/* calculate fluxes in theta direction */
void primtohflux(double *pr, double *f)
{
	double b[NDIM],tens[NDIM][NDIM],ucon[NDIM] ;
	void stress_calc(double *pr, double tens[][NDIM], double *b, double *ucon) ;
	int k ;

	stress_calc(pr,tens,b,ucon) ;

	/* simple fluxes */
	f[RHO] = pr[RHO]*pr[UH] ;
	f[UU]  = tens[TH][TT] + f[RHO] ;
	f[UR]  = tens[TH][RR] ;
	f[UH]  = tens[TH][TH] ;
	f[UP]  = tens[TH][PH] ;
	f[BR]  = b[RR]*pr[UH] - b[TH]*pr[UR] ;
	f[BH]  = 0. ;
	f[BP]  = b[PH]*pr[UH] - b[TH]*pr[UP] ;

	PLOOP f[k] *= g ;

}

/* calculate "conserved" quantities */
void primtoU(double *pr, double *U)
{
	double b[NDIM],tens[NDIM][NDIM],ucon[NDIM] ;
	void stress_calc(double *pr, double tens[][NDIM], double *b, double *ucon) ;
	int k ;

	stress_calc(pr,tens,b,ucon) ;

	U[RHO] = pr[RHO]*ucon[TT] ;
	U[UU]  = tens[TT][TT] + U[RHO] ;
	U[UR]  = tens[TT][RR] ;
	U[UH]  = tens[TT][TH] ;
	U[UP]  = tens[TT][PH] ;
	U[BR]  = pr[BR] ;
	U[BH]  = pr[BH] ;
	U[BP]  = pr[BP] ;

	PLOOP U[k] *= g ;

}

/* calculate components of MHD stress tensor */
void stress_calc(double *pr, double tens[][NDIM], double *b, double *ucon)
{

	int j,k ;
	double r,u,P,w ;
	double b2 ;
	void ucon_calc(double *pr, double *ucon) ;
	void b_calc(double *pr, double *ucon, double *b) ;
	double ucov[NDIM] ;
	double bcov[NDIM] ;
	double delta(int j, int k) ;
	void lower(double *a, double *b) ;

        r = pr[RHO] ;
        u = pr[UU] ;
        P = (gam - 1.)*u ;
        w = P + r + u ;

        /* find contravariant components of four-velocity */
	ucon_calc(pr, ucon) ;
	lower(ucon,ucov) ;

        /* find contravariant components of four-vector b */
	b_calc(pr, ucon, b) ;
	lower(b,bcov) ;

        /* find magnetic pressure */
        b2 = 0. ;
        DLOOP b2 += b[j]*b[k]*Gcov[j][k] ;

        /* calculate MHD stress tensor */
	DLOOP tens[j][k] = (w + b2)*ucon[j]*ucov[k] + (P + b2/2.)*delta(j,k) - b[j]*bcov[k] ;

}

/* find contravariant time component of four-velocity */
void ucon_calc(double *pr, double *ucon)
{
	int j,k ;
	double AA,BB,CC ;
	double discr ;

        ucon[RR] = pr[UR] ;
        ucon[TH] = pr[UH] ;
        ucon[PH] = pr[UP] ;

	AA = Gcov[TT][TT] ;
	BB = 0. ;
        SLOOPA BB += 2.*Gcov[TT][j]*ucon[j] ;
	CC = 1. ;
        SLOOP CC += Gcov[j][k]*ucon[j]*ucon[k] ;

	discr = BB*BB - 4.*AA*CC ;
	if(discr < 0.) {
		fprintf(stderr,"failure: spacelike four-velocity %g\n",
			discr) ;
		fprintf(stderr,"%d %d %d\n",icurr,jcurr,pcurr) ;
		fprintf(stderr,"%g %g\n",rcurr,hcurr) ;
		fprintf(stderr,"%15.8g %15.8g %15.8g\n",
			ucon[RR],ucon[TH],ucon[PH]) ;
		exit(1) ;
	}
	else
        	ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;


	return ;
}

/* calculate magnetic field four-vector */
void b_calc(double *pr, double *ucon, double *b) 
{
	int j,k ;
	double *B ;

	B = &pr[BR-1] ;

        b[TT] = 0. ;
	for(j=1;j<NDIM;j++)
	for(k=0;k<NDIM;k++)
		b[TT] += B[j]*ucon[k]*Gcov[j][k] ;

        SLOOPA b[j] = (B[j] + b[TT]*ucon[j])/ucon[TT] ;

	return ;
}

/* add in source terms to equations of motion */
void source(double *ph, double *dU)
{
	double tens[NDIM][NDIM],b[NDIM],ucon[NDIM] ;
	void stress_calc(double *pr, double tens[][NDIM], double *b, double *ucon) ;
	int j,k ;

	stress_calc(ph,tens,b,ucon) ;

	/* contract with connection; this mixed form is different
	   in sign and in which index of the connection the stress tensor
	   is contracted with.  */
	PLOOP dU[k] = 0. ;
	DLOOP {
		dU[UU] += tens[j][k]*Conn[k][TT][j] ;
		dU[UR] += tens[j][k]*Conn[k][RR][j] ;
		dU[UH] += tens[j][k]*Conn[k][TH][j] ;
		dU[UP] += tens[j][k]*Conn[k][PH][j] ;
	}

	PLOOP dU[k] *= g ;

	/* done! */
}

/* returns b^2 (i.e., twice magnetic pressure) */
double b2_calc(double *pr)
{
	int j,k ;
	double ucon[NDIM],b[NDIM],b2 ;
	void ucon_calc(double *pr, double *ucon) ;
	void b_calc(double *pr, double *ucon, double *b) ;

        /* find contravariant components of four-velocity */
	ucon_calc(pr, ucon) ;

        /* find contravariant components of four-vector b */
	b_calc(pr, ucon, b) ;

        /* find magnetic pressure */
        b2 = 0. ;
        DLOOP b2 += b[j]*b[k]*Gcov[j][k] ;

	return( b2 ) ;

}

void lower(double *ucon, double *ucov)
{
        int j,k ;

        DLOOPA ucov[j] = 0. ;
        DLOOP ucov[j] += Gcov[j][k]*ucon[k] ;

        return ;
}


