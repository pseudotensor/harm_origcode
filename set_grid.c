
#include "decs.h"

/* set up all grid functions */
void set_grid()
{
	int i,j ;
	void conn_func(double r, double th, double lconn[][NDIM][NDIM]) ;
	void gcov_func(double r, double th, double lgcov[][NDIM]) ;
	void gcon_func(double r, double th, double lgcon[][NDIM]) ;
	double gdet_func(double r, double th) ;
	void tetr_func(double tetr_cov[][NDIM],double tetr_con[][NDIM]) ;
	double r,th ;

	ZSLOOP(-2,NR+1,-2,NH+1) {
		
		/* zone-centered */
		r = Rin + (i + 0.5)*dr ;
		th = 0. + (j + 0.5)*dh ;
		gcov_func(r,th,gcov[i][j][CENT]) ;
		gcon_func(r,th,gcon[i][j][CENT]) ;
		gdet[i][j][CENT] = gdet_func(r,th) ;
		GSET(i,j,CENT) ;
		conn_func(r,th,conn[i][j]) ;
		tetr_func(tetr_cov[i][j][CENT],tetr_con[i][j][CENT]) ; 

		/* corner-centered */
		r = Rin + i*dr ;
		th = 0. + j*dh ;
		gcov_func(r,th,gcov[i][j][CORN]) ;
		gcon_func(r,th,gcon[i][j][CORN]) ;
		gdet[i][j][CORN] = gdet_func(r,th) ;
		GSET(i,j,CORN) ;
		tetr_func(tetr_cov[i][j][CORN],tetr_con[i][j][CORN]) ; 

		/* r-face-centered */
		r = Rin + i*dr ;
		th = 0. + (j + 0.5)*dh ;
		gcov_func(r,th,gcov[i][j][RFACE]) ;
		gcon_func(r,th,gcon[i][j][RFACE]) ;
		gdet[i][j][RFACE] = gdet_func(r,th) ;
		GSET(i,j,RFACE) ;
		tetr_func(tetr_cov[i][j][RFACE],tetr_con[i][j][RFACE]) ; 

		/* theta-face-centered */
		r = Rin + (i + 0.5)*dr ;
		th = 0. + j*dh ;
		gcov_func(r,th,gcov[i][j][HFACE]) ;
		gcon_func(r,th,gcon[i][j][HFACE]) ;
		gdet[i][j][HFACE] = gdet_func(r,th) ;
		GSET(i,j,HFACE) ;
		tetr_func(tetr_cov[i][j][HFACE],tetr_con[i][j][HFACE]) ; 

	}

	/* done! */
}
