
#include <stdio.h>
#include <math.h>

/** MNEMONICS SECTION **/

/* boundary condition mnemonics */
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3

/* mnemonics for primitive vars; conserved vars */
#define RHO	0	
#define UU	1
#define UR	2
#define UH	3
#define UP	4
#define BR	5
#define BH	6
#define BP	7

/* mnemonics for dimensional indices */
#define TT	0	
#define RR	1
#define TH	2
#define PH	3

/* mnemonics for centering of grid functions */
#define RFACE	0	
#define HFACE	1
#define CORN	2
#define CENT	3

/* mnemonics for slope limiter */
#define MC	0
#define VANL	1
#define MINM	2

/** GLOBAL ARRAY SECTION **/

/* size of global arrays */
#define NR	64	/* number of zones */
#define NH	64	/* number of zones */
#define NPR	8	/* number of primitive variables */
#define NDIM	4	/* number of total dimensions.  Never changes */
#define NPG	4	/* number of positions on grid for grid functions */

extern double   a_p[NR+4][NH+4][NPR] ;	/* space for primitive vars */
extern double a_dqr[NR+4][NH+4][NPR] ;	/* slopes */
extern double a_dqh[NR+4][NH+4][NPR] ;	/* slopes */
extern double  a_Fr[NR+4][NH+4][NPR] ;	/* fluxes */
extern double  a_Fh[NR+4][NH+4][NPR] ;	/* fluxes */
extern double  a_ph[NR+4][NH+4][NPR] ;	/* half-step primitives */

/* grid functions */
extern double a_conn[NR+4][NH+4][NDIM][NDIM][NDIM] ;
extern double a_gcon[NR+4][NH+4][NPG][NDIM][NDIM] ;
extern double a_gcov[NR+4][NH+4][NPG][NDIM][NDIM] ;
extern double a_gdet[NR+4][NH+4][NPG] ;
extern double a_tetr_cov[NR+4][NH+4][NPG][NDIM][NDIM] ;
extern double a_tetr_con[NR+4][NH+4][NPG][NDIM][NDIM] ;

extern double (*   p)[NH+4][NPR] ;
extern double (* dqr)[NH+4][NPR] ;
extern double (* dqh)[NH+4][NPR] ;
extern double (*  Fr)[NH+4][NPR] ;
extern double (*  Fh)[NH+4][NPR] ;
extern double (*  ph)[NH+4][NPR] ;
extern double (* conn)[NH+4][NDIM][NDIM][NDIM] ;
extern double (* gcon)[NH+4][NPG][NDIM][NDIM] ;
extern double (* gcov)[NH+4][NPG][NDIM][NDIM] ;
extern double (* gdet)[NH+4][NPG] ;
extern double (* tetr_cov)[NH+4][NPG][NDIM][NDIM] ;
extern double (* tetr_con)[NH+4][NPG][NDIM][NDIM] ;

/* global variables that should be set to the
   value at the current location on the grid */
extern double (* Gcon)[NDIM] ;
extern double (* Gcov)[NDIM] ;
extern double (* Conn)[NDIM][NDIM] ;
extern double g ;
extern double (* Tetr_cov)[NDIM] ;
extern double (* Tetr_con)[NDIM] ;

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern double a ;
extern double gam ;

/* numerical parameters */
extern double Rin,Rout ;
extern double cour ;
extern double dr,dh,dV ;
extern double dt ;
extern double t,tf ;
extern double rcurr,hcurr ;

/* output parameters */
extern double DTd ;
extern double DTl ;
extern double DTi ;
extern int    DTr ;

/* global flags */
extern int fail ;
extern int lim ;

/* diagnostics */
extern double mdot ;
extern double edot ;
extern double ldot ;

/* numerical convenience */
#define SMALL	1.e-20 

/** MACROS **/

/* loop over all active zones */
#define ZLOOP for(i=0;i<NR;i++)for(j=0;j<NH;j++)

/* specialty loop */
extern int istart,istop,jstart,jstop ;
#define ZSLOOP(istart,istop,jstart,jstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

/* loop over Primitive variables */
#define PLOOP for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)

/* set global variables that indicate current local metric, etc. */
extern int icurr,jcurr,pcurr ;
extern int ihere,jhere,phere ;
#define GSET(ihere,jhere,phere) {Gcon = gcon[ihere][jhere][phere]; \
        Gcov = gcov[ihere][jhere][phere] ;\
        g = gdet[ihere][jhere][phere] ;\
        Conn = conn[ihere][jhere];\
	Tetr_cov = tetr_cov[ihere][jhere][phere];\
	Tetr_con = tetr_con[ihere][jhere][phere];\
        icurr = ihere ;\
	jcurr = jhere ;\
	pcurr = phere ;\
	rcurr = Rin + (ihere+0.5)*dr;\
	hcurr = 0. + (jhere+0.5)*dh;\
	}


extern double dminarg1,dminarg2; 
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

/* size of step in numerical derivative evaluations */
#define HSTEP	1.e-5

/** FIXUP PARAMETERS **/
#define RHOMIN	1.e-4
#define UUMIN	1.e-6

