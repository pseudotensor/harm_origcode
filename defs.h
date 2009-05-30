

double   a_p[NR+4][NH+4][NPR] ;  /* space for primitive vars */
double a_dqr[NR+4][NH+4][NPR] ;  /* slopes */
double a_dqh[NR+4][NH+4][NPR] ;  /* slopes */
double  a_Fr[NR+4][NH+4][NPR] ;  /* fluxes */
double  a_Fh[NR+4][NH+4][NPR] ;  /* fluxes */
double  a_ph[NR+4][NH+4][NPR] ;  /* half-step primitives */

/* grid functions */
double a_conn[NR+4][NH+4][NDIM][NDIM][NDIM] ;
double a_gcon[NR+4][NH+4][NPG][NDIM][NDIM] ;
double a_gcov[NR+4][NH+4][NPG][NDIM][NDIM] ;
double a_gdet[NR+4][NH+4][NPG] ;
double a_tetr_cov[NR+4][NH+4][NPG][NDIM][NDIM] ;
double a_tetr_con[NR+4][NH+4][NPG][NDIM][NDIM] ;

double (*   p)[NH+4][NPR] ;
double (* dqr)[NH+4][NPR] ;
double (* dqh)[NH+4][NPR] ;
double (*  Fr)[NH+4][NPR] ;
double (*  Fh)[NH+4][NPR] ;
double (*  ph)[NH+4][NPR] ;
double (* conn)[NH+4][NDIM][NDIM][NDIM] ;
double (* gcon)[NH+4][NPG][NDIM][NDIM] ;
double (* gcov)[NH+4][NPG][NDIM][NDIM] ;
double (* gdet)[NH+4][NPG] ;
double (* tetr_cov)[NH+4][NPG][NDIM][NDIM] ;
double (* tetr_con)[NH+4][NPG][NDIM][NDIM] ;

/* global variables that should be set to the
   value at the current location on the grid */
double (* Gcon)[NDIM] ;
double (* Gcov)[NDIM] ;
double (* Conn)[NDIM][NDIM] ;
double g ;
double (* Tetr_cov)[NDIM] ;
double (* Tetr_con)[NDIM] ;

/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
double a ;
double gam ;

/* numerical parameters */
double Rin,Rout ;
double cour ;
double dr,dh,dV ;
double dt ;
double t,tf ;
double rcurr,hcurr ;
int istart,istop,jstart,jstop ;
int icurr,jcurr,pcurr,ihere,jhere,phere ;
double dminarg1,dminarg2 ;

/* output parameters */
double DTd ;
double DTl ;
double DTi ;
int    DTr ;

/* global flags */
int    fail ;
int    lim ;

/* diagnostics */
double mdot = 0. ;
double edot = 0. ;
double ldot = 0. ;

