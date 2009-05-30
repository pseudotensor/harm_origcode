
/* 
 * calculate du/dp analytically.  
 * p is the full (NPR) element primitive variables
 * alpha is be a 5x5 matrix containing dU(1-5)/dp(1-5).
 * used only by nutoprim
 * 
 *
 * cfg 7-11-01
 *
 * cfg 8-23-01 heavy repairs
 *
 */

#include "decs.h"

void dudp_calc(double *pr, double **alpha)
{
	static double ucon[NDIM] ;
	static double b[NDIM] ;
	static double dutdui[NDIM] ;
	static double dbtdui[NDIM] ;
	static double db2dui[NDIM] ;
	static double dbiduj[NDIM][NDIM] ;
	static double duudud[NDIM][NDIM] ;
	static double tmp1[NDIM],tmp2[NDIM] ;
	static double tmp3[6][6] ;
	double w,b2 ;
	int j,k,l ;
	void ucon_calc(double *pr, double *ucon) ;
	void b_calc(double *pr, double *ucon, double *b) ;
	void dutdui_calc(double *ucon,double *dutdui) ;
	void dbtdui_calc(double *dutdui,double *pr,double *dbtdui) ;
	void dbiduj_calc(double *dbtdui,double *dutdui,double *ucon,double *b,double dbiduj[][NDIM]) ;
	void db2dui_calc(double dbiduj[][NDIM],double *b,double *db2dui) ;
	void duudud_calc(double *ucon, double duudud[][NDIM]) ;
	double contract(double *vcon, double *wcon) ;
	void lower(double *v1, double *v2) ;
	void raise(double *v1, double *v2) ;

	for(j=1;j<=5;j++)
	for(k=1;k<=5;k++) alpha[j][k] = 0. ;

	ucon_calc(pr,ucon) ;

	/*
	fprintf(stderr,"rcurr,hcurr: %g %g\n",rcurr,hcurr) ;
	fprintf(stderr,"%15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
			pr[0], pr[1], pr[2], pr[3], pr[4], pr[5], pr[6], pr[7]) ; 
	*/

	alpha[RHO+1][RHO+1] = ucon[TT] ;
	alpha[UU+1][RHO+1] = ucon[TT]*ucon[TT] ;
	alpha[UR+1][RHO+1] = ucon[TT]*ucon[RR] ;
	alpha[UH+1][RHO+1] = ucon[TT]*ucon[TH] ;
	alpha[UP+1][RHO+1] = ucon[TT]*ucon[PH] ;

	alpha[UU+1][UU+1] = gam*ucon[TT]*ucon[TT] + (gam - 1.)*Gcon[TT][TT] ;
	alpha[UR+1][UU+1] = gam*ucon[TT]*ucon[RR] + (gam - 1.)*Gcon[TT][RR] ;
	alpha[UH+1][UU+1] = gam*ucon[TT]*ucon[TH] + (gam - 1.)*Gcon[TT][TH] ;
	alpha[UP+1][UU+1] = gam*ucon[TT]*ucon[PH] + (gam - 1.)*Gcon[TT][PH] ;

	b_calc(pr,ucon,b) ;
	dutdui_calc(ucon,dutdui) ;
	dbtdui_calc(dutdui,pr,dbtdui) ;
	dbiduj_calc(dbtdui,dutdui,ucon,b,dbiduj) ;
	db2dui_calc(dbiduj,b,db2dui) ;
	b2 = contract(b,b) ;
	w = pr[RHO] + gam*pr[UU] + b2 ;

	alpha[RHO+1][UR+1] = pr[0]*dutdui[RR] ;
	alpha[UU+1][UR+1] = db2dui[RR]*ucon[TT]*ucon[TT] 
		+ 2.*w*dutdui[RR]*ucon[TT]
		+ 0.5*db2dui[RR]*Gcon[TT][TT] 
		- 2.*dbtdui[RR]*b[TT]  ;
	alpha[UR+1][UR+1] = db2dui[RR]*ucon[TT]*ucon[RR] 
		+ w*dutdui[RR]*ucon[RR]
		+ w*ucon[TT]
		+ 0.5*db2dui[RR]*Gcon[TT][RR] 
		- dbtdui[RR]*b[RR] 
		- b[TT]*dbiduj[RR][RR] ;
	alpha[UH+1][UR+1] = db2dui[RR]*ucon[TT]*ucon[TH] 
		+ w*dutdui[RR]*ucon[TH]
		+ 0.5*db2dui[RR]*Gcon[TT][TH] 
		- dbtdui[RR]*b[TH] 
		- b[TT]*dbiduj[TH][RR] ;
	alpha[UP+1][UR+1] = db2dui[RR]*ucon[TT]*ucon[PH] 
		+ w*dutdui[RR]*ucon[PH]
		+ 0.5*db2dui[RR]*Gcon[TT][PH] 
		- dbtdui[RR]*b[PH] 
		- b[TT]*dbiduj[PH][RR] ;

	alpha[RHO+1][UH+1] = pr[0]*dutdui[TH] ;
	alpha[UU+1][UH+1] = db2dui[TH]*ucon[TT]*ucon[TT] 
		+ 2.*w*dutdui[TH]*ucon[TT]
		+ 0.5*db2dui[TH]*Gcon[TT][TT] 
		- dbtdui[TH]*b[TT] 
		- b[TT]*dbiduj[TT][TH] ;
	alpha[UR+1][UH+1] = db2dui[TH]*ucon[TT]*ucon[RR] 
		+ w*dutdui[TH]*ucon[RR]
		+ 0.5*db2dui[TH]*Gcon[TT][RR] 
		- dbtdui[TH]*b[RR] 
		- b[TT]*dbiduj[RR][TH] ;
	alpha[UH+1][UH+1] = db2dui[TH]*ucon[TT]*ucon[TH] 
		+ w*dutdui[TH]*ucon[TH]
		+ 0.5*db2dui[TH]*Gcon[TT][TH] 
		+ w*ucon[TT]
		- dbtdui[TH]*b[TH] 
		- b[TT]*dbiduj[TH][TH] ;
	alpha[UP+1][UH+1] = db2dui[TH]*ucon[TT]*ucon[PH] 
		+ w*dutdui[TH]*ucon[PH]
		+ 0.5*db2dui[TH]*Gcon[TT][PH] 
		- dbtdui[TH]*b[PH] 
		- b[TT]*dbiduj[PH][TH] ;

	alpha[RHO+1][UP+1] = pr[0]*dutdui[PH] ;
	alpha[UU+1][UP+1] = db2dui[PH]*ucon[TT]*ucon[TT] 
		+ 2.*w*dutdui[PH]*ucon[TT]
		+ 0.5*db2dui[PH]*Gcon[TT][TT] 
		- dbtdui[PH]*b[TT] 
		- b[TT]*dbiduj[TT][PH] ;
	alpha[UR+1][UP+1] = db2dui[PH]*ucon[TT]*ucon[RR] 
		+ w*dutdui[PH]*ucon[RR]
		+ 0.5*db2dui[PH]*Gcon[TT][RR] 
		- dbtdui[PH]*b[RR] 
		- b[TT]*dbiduj[RR][PH] ;
	alpha[UH+1][UP+1] = db2dui[PH]*ucon[TT]*ucon[TH] 
		+ w*dutdui[PH]*ucon[TH]
		+ 0.5*db2dui[PH]*Gcon[TT][TH] 
		- dbtdui[PH]*b[TH] 
		- b[TT]*dbiduj[TH][PH] ;
	alpha[UP+1][UP+1] = db2dui[PH]*ucon[TT]*ucon[PH] 
		+ w*dutdui[PH]*ucon[PH]
		+ 0.5*db2dui[PH]*Gcon[TT][PH] 
		+ w*ucon[TT]
		- dbtdui[PH]*b[PH] 
		- b[TT]*dbiduj[PH][PH] ;

	/* mixed indices version of code, i.e. conserved variables
	   2-5 are T^t_j rather than T^{tj}.  Convert from original
	   results for unmixed indices */
	for(k=1;k<=5;k++) {
		for(j=2;j<=5;j++) tmp1[j-2] = alpha[j][k] ;
		lower(tmp1,tmp2) ;
		for(j=2;j<=5;j++) alpha[j][k] = tmp2[j-2] ;
	}

	/* this bit of legacy code can be uncommented if
	   the rest-mass flux is subtracted out from the energy 
	   flux, which may be numerically convenient (tests are unconclusive) */
	j = 2 ; for(k=1;k<=5;k++) alpha[j][k] += alpha[1][k] ;
	/*
	*/

	/* N.B.: all the conserved variables contain a factor of \sqrt{det(g_{\mu\nu})} */
	for(j=1;j<=5;j++)
	for(k=1;k<=5;k++) alpha[j][k] *= g ;

	return ;
}

void dutdui_calc(double *ucon,double *dutdui) 
{
	int j ;
	double ucov[NDIM] ;
	void lower(double *ucon, double *ucov) ;

	lower(ucon,ucov) ;
	SLOOPA dutdui[j] = -ucov[j]/ucov[0] ;

	return ;
}

void dbtdui_calc(double *dutdui, double *pr, double *dbtdui) 
{
	int j ;
	double B[NDIM],Bcov[NDIM] ;
	void lower(double *ucon, double *ucov) ;

	B[0] = 0. ;
	SLOOPA B[j] = pr[j+BR-1] ;

	lower(B,Bcov) ;
	SLOOPA dbtdui[j] = Bcov[j] + Bcov[0]*dutdui[j] ;

	return ;
}

void dbiduj_calc(double *dbtdui,double *dutdui,double *ucon, double *b, double dbiduj[][NDIM]) 
{
	int j,k ;

	DLOOP dbiduj[j][k] = 0. ;

	SLOOP dbiduj[j][k] = -b[j]*dutdui[k]/ucon[TT] 
			+ ucon[j]*dbtdui[k]/ucon[TT] ;

	SLOOPA dbiduj[j][j] += b[TT]/ucon[TT] ;

	SLOOPA dbiduj[TT][j] = dbtdui[j] ;

	return ;
}

void db2dui_calc(double dbiduj[][NDIM],double *b,double *db2dui) 
{
	int j,k ;
	double bcov[NDIM] ;
	void lower(double *ucon, double *ucov) ;

	lower(b,bcov) ;
	DLOOPA db2dui[j] = 0. ;
	DLOOP db2dui[j] += 2.*bcov[k]*dbiduj[k][j] ;

	return ;
}

double contract(double *vcon1, double *vcon2)
{
	int j,k ;
	double n ;

	n = 0. ;
	DLOOP n += Gcov[j][k]*vcon1[j]*vcon2[k] ;
	return(n) ;

}

double covtract(double *vcov1, double *vcov2)
{
	int j,k ;
	double n ;

	n = 0. ;
	DLOOP n += Gcon[j][k]*vcov1[j]*vcov2[k] ;
	return(n) ;

}

void duudud_calc(double *ucon, double duudud[][NDIM]) 
{
	int i,j ;
	void lower(double *vcon, double *vcov) ;
	double ucov[NDIM] ;

	lower(ucon,ucov) ;

	for(i=1;i<=3;i++)
	for(j=1;j<=3;j++) 
		duudud[i][j] = Gcon[i][j] + Gcon[i][TT]*(-ucon[j]/ucon[TT]) ;

}
