
#include "decs.h"

/* diagnostics subroutine */

void diag(call_code)
int call_code ;
{
	char dfnam[100],ifnam[100] ;
	int i,j ;
	FILE *dump_file,*image_file ;
	void dump(FILE *fp) ;
	void image(FILE *fp) ;
	double U[NPR],pp,e,rmed,divb,divbmax,e_fin,m_fin ;
	int imax,jmax ;
	void primtoU(double *p, double *U) ;

	static int dump_cnt,image_cnt ;
	static double e_init,m_init,tdump,timage ;
	static FILE *ener_file ;

	if(call_code==0) {
		ener_file = fopen("ener.out","w") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}
		tdump = t + SMALL ;
		timage = t + SMALL ;
		dump_cnt = 0 ;
		image_cnt = 0 ;
	}

	/* calculate conserved quantities */
	pp = 0. ;
	e = 0. ;
	rmed = 0. ;
	divbmax = 0. ;
	imax = 0 ; 
	jmax = 0 ;
	ZLOOP {
		g = gdet[i][j][CENT] ;
		Gcon = gcon[i][j][CENT] ;
		Gcov = gcov[i][j][CENT] ;
		primtoU(p[i][j],U) ;

		rmed += U[RHO]*g*dV ;
		pp += U[UP]*g*dV ;
		e += U[UU]*g*dV ;

		/* flux-ct defn */
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

		if(divb > divbmax && i > 0 && j > 0) {
			imax = i ;
			jmax = j ;
			divbmax = divb ;
		}
	}
	fprintf(stderr,"divbmax: %d %d %g\n",imax,jmax,divbmax) ;

	if(call_code == 0) {
		e_init = e ;
		m_init = rmed ;
	}
	else if (call_code == 2) {
		e_fin = e ;
		m_fin = rmed ;
		fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
			e_init,e_fin,(e_fin-e_init)/e_init) ;
		fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
			m_init,m_fin,(m_fin-m_init)/m_init) ;
	}
	else {
		fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
			t,rmed,pp,e,p[NR/2][NH/2][UU]*pow(p[NR/2][NH/2][RHO],-gam),
			p[NR/2][NH/2][UU]) ;
		fprintf(ener_file,"%15.8g %15.8g %15.8g ",mdot,edot,ldot) ;
		fprintf(ener_file,"\n") ;
		fflush(ener_file) ;
	}


	/* dump at regular intervals */
	if(t >= tdump || call_code == 2 || call_code == 0) {
		fprintf(stderr,"dumping.\n") ;
		/* make regular dump file */
		sprintf(dfnam,"dumps/dump%03d",dump_cnt) ;
		dump_file = fopen(dfnam,"w") ;

		if(dump_file==NULL) {
			fprintf(stderr,"error opening dump file\n") ;
			exit(2) ;
		}

		dump(dump_file) ;
		fclose(dump_file) ;

		dump_cnt++ ;
		tdump += DTd ;
	}
	
	/* image dump at regular intervals */
	if(t >= timage) {
		/* make regular dump file */
		sprintf(ifnam,"images/im%04d",image_cnt) ;
		image_file = fopen(ifnam,"w") ;

		if(image_file==NULL) {
			fprintf(stderr,"error opening image file\n") ;
			exit(2) ;
		}

		image(image_file) ;
		fclose(image_file) ;

		image_cnt++ ;
		timage += DTi ;
	}
}
