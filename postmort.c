#include "decs.h"
#include "defs.h"

main()
{
	FILE *fp ;
	void init() ;
	double dum1,dum2,pr[NPR],U[NPR] ;
	int idum,i,j,k ; 
	void Utoprim(double *U, double *pr) ;

	init() ;

	fp = fopen("postmort.dat","r") ;

	fscanf(fp,"%d %d\n",&i,&j) ;
	for(k=0;k<NPR;k++) {
		fscanf(fp,"%d %lf %lf %lf %lf\n",
				&idum,&dum1,&pr[k],&U[k],&dum2) ;
	}

	fprintf(stderr,"\n") ;
	for(k=0;k<NPR;k++)
		fprintf(stderr,"%d %g %g\n",k,pr[k],U[k]) ;


	GSET(i,j,CENT) ;
	Utoprim(U,pr) ;

	fprintf(stderr,"\n") ;
	for(k=0;k<NPR;k++)
		fprintf(stderr,"%d %g %g\n",k,pr[k],U[k]) ;

	fclose(fp) ;
	exit(0) ;
}
