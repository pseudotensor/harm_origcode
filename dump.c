
#include "decs.h"

void dump(FILE *fp)
{
	int i,j,k,l ;
	double divb,b2_calc(double *pr) ;
	void primtoU(double *pr, double *U) ;
	double vrchar(double *pr) ;
	double vhchar(double *pr) ;
	double tens_em[NDIM][NDIM],tens_matt[NDIM][NDIM],b[NDIM],ucon[NDIM] ;
        void sp_stress_calc(double *pr, double tens_matt[][NDIM], double tens_em[][NDIM], double *b, double *ucon) ;
	double U[NPR] ;

	fprintf(fp,"%d %d %10.5g %10.5g %10.5g\n",NR,NH,Rin,Rout,t) ;

	/*
	ZSLOOP(-2,NR+1,-2,NR+1) {
	*/
	ZSLOOP(0,NR-1,0,NR-1) {

		fprintf(fp,"%15.7g %15.7g",Rin + (i + 0.5)*dr,(j+0.5)*dh) ;
		PLOOP fprintf(fp,"%15.7g ",p[i][j][k]) ;

                /* flux-ct defn; corner-centered.  Use
		   only interior corners */
		if(i > 0 && j > 0 && i < NR && j < NH) {
			divb = fabs( 0.5*(
				p[i][j][BR]*gdet[i][j][CENT]
				+ p[i][j-1][BR]*gdet[i][j-1][CENT]
				- p[i-1][j][BR]*gdet[i-1][j][CENT]
				- p[i-1][j-1][BR]*gdet[i-1][j-1][CENT]
				)/dr +
				0.5*(
				p[i][j][BH]*gdet[i][j][CENT]
				+ p[i-1][j][BH]*gdet[i-1][j][CENT]
				- p[i][j-1][BH]*gdet[i][j-1][CENT]
				- p[i-1][j-1][BH]*gdet[i-1][j-1][CENT]
				)/dh) ;
		}
		else divb = 0. ;

		fprintf(fp,"%15.7g ",divb) ;

		fprintf(fp,"%15.7g ",b2_calc(p[i][j])) ;

		sp_stress_calc(p[i][j],tens_matt,tens_em,b,ucon) ;

		fprintf(fp,"%15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g %15.7g ",
			tens_matt[0][0], tens_em[0][0], tens_matt[0][3], tens_em[0][3],
			tens_matt[1][0], tens_em[1][0], tens_matt[1][3], tens_em[1][3],
			tens_matt[2][0], tens_em[2][0], tens_matt[2][3], tens_em[2][3]) ; 
		/*
		primtoU(p[i][j],U) ;
		PLOOP fprintf(fp,"%15.7g ", U[k]) ;
		*/

		/*
		GSET(i,j,CENT)
		fprintf(fp,"%15.7g %15.7g ",vrchar(p[i][j]),vhchar(p[i][j])) ;

		for(k=0;k<NDIM;k++)
		for(l=0;l<NDIM;l++) fprintf(fp,"%15.7g ",Tetr_con[k][l]) ;
		*/
                
		fprintf(fp,"\n") ;
	}
}

/* calculate components of MHD stress tensor */
void sp_stress_calc(double *pr, double tens_matt[][NDIM], double tens_em[][NDIM], double *b, double *ucon)
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
        DLOOP tens_em[j][k] = - b[j]*bcov[k] ;
        DLOOP tens_matt[j][k] = (w + b2)*ucon[j]*ucov[k] + (P + b2/2.)*delta(j,k) ;

}

