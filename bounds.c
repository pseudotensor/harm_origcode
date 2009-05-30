
#include "decs.h"

/* bound array containing entire set of primitive variables */

void bound_prim( double prim[][NH+4][NPR] )
{
	void bound_var( double var[][NH+4][NPR], int k, int r_in_bc, int r_out_bc,
		int th_in_bc, int th_out_bc) ;
	void appl_infl_switch_out( double var[][NH+4][NPR] ) ;

	bound_var(prim, RHO, OUTFLOW, OUTFLOW, SYMM,  SYMM) ;
	bound_var(prim, UU,  OUTFLOW, OUTFLOW, SYMM,  SYMM) ;
 
	bound_var(prim, UR,  OUTFLOW, OUTFLOW, SYMM,  SYMM) ;
	bound_var(prim, UH,  OUTFLOW, OUTFLOW, ASYMM, ASYMM) ;
	bound_var(prim, UP,  OUTFLOW, OUTFLOW, SYMM,  SYMM) ;

	bound_var(prim, BR,  OUTFLOW, OUTFLOW, SYMM,  SYMM) ;
	bound_var(prim, BH,  OUTFLOW, OUTFLOW, ASYMM, ASYMM) ;
	bound_var(prim, BP,  OUTFLOW, OUTFLOW, SYMM,  SYMM) ;

	/* apply switches to ghost zones */
	appl_infl_switch_out(prim) ;
}

#define UNK_BC	{fprintf(stderr,"unknown BC\n") ; exit(1) ;}

/* bound individual primitive variable */

void bound_var( double var[][NH+4][NPR], int k, 
	int r_in_bc, int r_out_bc,
	int th_in_bc, int th_out_bc
	)
{
	int i,j ;

	/* inner r boundary condition */
	if(r_in_bc == OUTFLOW) {
		for(j=0;j<NH;j++) {
			var[-1][j][k] = var[0][j][k] 
				* (gdet[0][j][CENT]/gdet[-1][j][CENT]) ;
			var[-2][j][k] = var[0][j][k] 
				* (gdet[0][j][CENT]/gdet[-2][j][CENT]) ;
		}
	}
	else UNK_BC ;

	/* outer r BC */
	if(r_out_bc == OUTFLOW) {
		for(j=0;j<NH;j++) {
			var[NR][j][k] = var[NR-1][j][k] ;
			var[NR+1][j][k] = var[NR-1][j][k] ;
		}
	}
	else UNK_BC ;

	/* north pole BC */
	if(th_in_bc == SYMM || th_in_bc == ASYMM) {
		for(i=-2;i<=NR+1;i++) {
			var[i][-1][k] = var[i][0][k] ;
			var[i][-2][k] = var[i][1][k] ;
		}
	}
	else UNK_BC ;

	/* south pole BC */
	if(th_out_bc == SYMM || th_out_bc == ASYMM) {
		for(i=-2;i<=NR+1;i++) {
			var[i][NH][k]   = var[i][NH-1][k] ;
			var[i][NH+1][k] = var[i][NH-2][k] ;
		}
	}
	else UNK_BC ;

}

/* inflow switch for boundary zones.  Note that it must
   be applied to ucon, *not* ucov */
void appl_infl_switch_out( double var[][NH+4][NPR] ) 
{
	int i,j ;
	void ucon_calc(double *pr, double *ucon) ;
	void b_calc(double *pr, double *ucon, double *b) ;
	double ucon[NDIM],b[NDIM] ;

	for(j=-2;j<NH+2;j++) 
	for(i=NR;i<NR+2;i++) {
		GSET(i,j,CENT) ;
		ucon_calc(var[i][j],ucon) ;
		if(ucon[RR] < 0.) {
			ucon[RR] = 0. ;
			var[i][j][UR] = ucon[RR] ;
			var[i][j][UH] = ucon[TH] ;
			var[i][j][UP] = ucon[PH] ;
		}
	}

	for(i=-2;i<NR+2;i++) {
		for(j=-2;j<0;j++) {
			GSET(i,j,CENT) ;
			ucon_calc(var[i][j],ucon) ;
			b_calc(var[i][j],ucon,b) ;
			ucon[TH] *= -1. ;
			b[TH]    *= -1. ;

			var[i][j][UR] = ucon[RR] ;
			var[i][j][UH] = ucon[TH] ;
			var[i][j][UP] = ucon[PH] ;

			var[i][j][BR] = b[RR]*ucon[TT] - ucon[RR]*b[TT] ;
			var[i][j][BH] = b[TH]*ucon[TT] - ucon[TH]*b[TT] ;
			var[i][j][BP] = b[PH]*ucon[TT] - ucon[PH]*b[TT] ;
		}
		for(j=NH;j<NH+2;j++) {
			GSET(i,j,CENT) ;
			ucon_calc(var[i][j],ucon) ;
			b_calc(var[i][j],ucon,b) ;
			ucon[TH] *= -1. ;
			b[TH]    *= -1 ;

			var[i][j][UR] = ucon[RR] ;
			var[i][j][UH] = ucon[TH] ;
			var[i][j][UP] = ucon[PH] ;

			var[i][j][BR] = b[RR]*ucon[TT] - ucon[RR]*b[TT] ;
			var[i][j][BH] = b[TH]*ucon[TT] - ucon[TH]*b[TT] ;
			var[i][j][BP] = b[PH]*ucon[TT] - ucon[PH]*b[TT] ;
		}
	}
}
