
#include "decs.h"

void set_arrays()
{
	int i,j,k ;

	p =   (double (*) [NH+4][NPR])(& (  a_p[2][2][0])) ;
	dqr = (double (*) [NH+4][NPR])(& (a_dqr[2][2][0])) ;
	dqh = (double (*) [NH+4][NPR])(& (a_dqh[2][2][0])) ;
	Fr =  (double (*) [NH+4][NPR])(& ( a_Fr[2][2][0])) ;
	Fh =  (double (*) [NH+4][NPR])(& ( a_Fh[2][2][0])) ;
	ph =  (double (*) [NH+4][NPR])(& ( a_ph[2][2][0])) ;

	/* everything must be initialized to zero */
	ZSLOOP(-2,NR+1,-2,NH+1) PLOOP {
			p[i][j][k]   = 0. ;
			ph[i][j][k]  = 0. ;
			dqr[i][j][k] = 0. ;
			dqh[i][j][k] = 0. ;
			Fr[i][j][k]  = 0. ;
			Fh[i][j][k]  = 0. ;
	}

	/* grid functions */
	conn = (double (*) [NH+4][NDIM][NDIM][NDIM])
		(& ( a_conn[2][2][0][0][0])) ;
	gcon = (double (*) [NH+4][NPG][NDIM][NDIM])
		(& ( a_gcon[2][2][0][0][0])) ;
	gcov = (double (*) [NH+4][NPG][NDIM][NDIM])
		(& ( a_gcov[2][2][0][0][0])) ;
	gdet = (double (*) [NH+4][NPG])
		(& ( a_gdet[2][2][0])) ;
	tetr_cov = (double (*) [NH+4][NPG][NDIM][NDIM])
		(& ( a_tetr_cov[2][2][0][0][0])) ;
	tetr_con = (double (*) [NH+4][NPG][NDIM][NDIM])
		(& ( a_tetr_con[2][2][0][0][0])) ;
	

}
