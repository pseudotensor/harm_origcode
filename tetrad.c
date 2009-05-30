
#include "decs.h"

#if 1
/* find suitable orthonormal tetrad */

void tetr_func(double tetr_cov[][NDIM],double tetr_con[][NDIM]) 
{
	void tet_func(double metr[][NDIM], double tetr[][NDIM]) ;
	static double tmp[NDIM][NDIM] ;
	double mink(int i, int j) ;
	int j,k,l ;

	tet_func(Gcov,tetr_con) ;

	DLOOP tmp[j][k] = 0. ;
	DLOOP for(l=0;l<NDIM;l++) tmp[j][k] += tetr_con[j][l]*Gcov[l][k] ;
	DLOOP tetr_cov[j][k] = 0. ;
	DLOOP for(l=0;l<NDIM;l++) tetr_cov[j][k] += mink(j,l)*tmp[l][k] ;

#if 0
	/* tetr_ cov & con are inverse transposes of each other */
	DLOOP tmp[j+1][k+1] = tetr_cov[j][k] ;
	gaussj(tmp,NDIM,NULL,0) ;
	DLOOP tetr_con[j][k] = tmp[k+1][j+1] ;
#endif

	return ;
}

void tet_func(double metr[][NDIM], double tetr[][NDIM])
{
	char jobz,uplo ;
	int n,lda,lwork,info ;
	double a[NDIM][NDIM],w[NDIM],work[NDIM*NDIM*NDIM] ;
	int chk,liwork,iwork[5+5*NDIM] ;
	int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, 
			double *w, double *work, int *lwork, 
			int *iwork, int *liwork, int *info) ;
	int j,k ;

	jobz = 'V' ;
	uplo = 'U' ;
	n = NDIM ;
	lda = NDIM ;
	lwork = NDIM*NDIM*NDIM ;
	liwork = 5 + 5*NDIM ;

	DLOOP a[j][k] = metr[j][k] ;

	chk = dsyev_(
		&jobz, 		/* job: 'V' -> compute eigenvectors too */
		&uplo,		/* which part of a is stored, 'U' -> upper */
		&n,		/* order of matrix a */
		(double *)a,	/* matrix (row major order) */
		&lda,		/* leading dimension of a */
		w,		/* eigenvalues, ascending order */
		work,		/* workspace */
		&lwork,		/* size of workspace */
		iwork,		/* size of iwork */
		&liwork,	/* working array for optimal liwork */
		&info		/* successful? */
		) ;

	DLOOP tetr[j][k] = a[j][k]/sqrt(fabs(w[j])+SMALL) ;

	/*
	DLOOP fprintf(stderr,"chk3 %d %d %g %g\n",j,k,a[j][k],w[j]) ;
	*/

}
#endif
