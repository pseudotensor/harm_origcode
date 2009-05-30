
/* restart functions; restart_init and restart_dump */

#include "decs.h"

void restart_write()
{
	FILE *fp ;
	int idum,i,j,k ;

	fp = fopen("dumps/rdump","w") ;
	if(fp == NULL) {
		fprintf(stderr,"can't open restart file\n") ;
		exit(2) ;
	}

	/* write out key global variables, in binary */
	idum = NR ;
	fwrite(&idum, sizeof(int),   1,fp) ;
	idum = NH ;
	fwrite(&idum, sizeof(int),   1, fp) ;
	fwrite(&t,    sizeof(double),1, fp) ;
	fwrite(&tf,   sizeof(double),1, fp) ;
	fwrite(&a,    sizeof(double),1, fp) ;
	fwrite(&gam,  sizeof(double),1, fp) ;
	fwrite(&Rin,  sizeof(double),1, fp) ;
	fwrite(&Rout, sizeof(double),1, fp) ;
	fwrite(&cour, sizeof(double),1, fp) ;
	fwrite(&dr,   sizeof(double),1, fp) ;
	fwrite(&dh,   sizeof(double),1, fp) ;
	fwrite(&dV,   sizeof(double),1, fp) ;
	fwrite(&dt,   sizeof(double),1, fp) ;
	fwrite(&DTd,  sizeof(double),1, fp) ;
	fwrite(&DTl,  sizeof(double),1, fp) ;
	fwrite(&DTi,  sizeof(double),1, fp) ;
	fwrite(&DTr,  sizeof(int),   1, fp) ;

	ZSLOOP(0,NR-1,0,NH-1) PLOOP fwrite(&p[i][j][k],    sizeof(double), NPR, fp) ;

	fclose(fp) ;

}

int restart_init()
{
	FILE *fp ;
	char ans[100] ;
	void set_arrays(), set_grid(), bound_prim(double pr[][NH+4][NPR]) ;
	void restart_read(FILE *fp) ;
	int strncmp(char *s1, char *s2, int n) ;

	/* set up global arrays */
	set_arrays() ;

	fp = fopen("dumps/rdump","r") ;
	if(fp == NULL) {
		fprintf(stderr,"no restart file\n") ;
		return(0) ;
	}
	else {
		fprintf(stderr,"restart file exists.  Use? [y|n]\n") ;
		fscanf(stdin,"%s",ans) ;
		if(strncmp(ans,"no",1) == 0) return(0) ;
	}

	restart_read(fp) ;
	fclose(fp) ;

	/* set metric functions */
	set_grid() ;

	/* bound */
	bound_prim(p) ;

	/* done! */
	return(1) ;

}

void restart_read(FILE *fp)
{
	int idum,i,j,k ;

	/* read in global variables, in binary */
	fread(&idum, sizeof(int),   1,fp) ;
	if(idum != NR) {
		fprintf(stderr,"error reading restart file; NR differs\n") ;
		exit(3) ;
	}
	fread(&idum, sizeof(int),   1, fp) ;
	if(idum != NH) {
		fprintf(stderr,"error reading restart file; NH differs\n") ;
		exit(4) ;
	}
	fread(&t,    sizeof(double),1, fp) ;
	fread(&tf,   sizeof(double),1, fp) ;
	fread(&a,    sizeof(double),1, fp) ;
	fread(&gam,  sizeof(double),1, fp) ;
	fread(&Rin,  sizeof(double),1, fp) ;
	fread(&Rout, sizeof(double),1, fp) ;
	fread(&cour, sizeof(double),1, fp) ;
	fread(&dr,   sizeof(double),1, fp) ;
	fread(&dh,   sizeof(double),1, fp) ;
	fread(&dV,   sizeof(double),1, fp) ;
	fread(&dt,   sizeof(double),1, fp) ;
	fread(&DTd,  sizeof(double),1, fp) ;
	fread(&DTl,  sizeof(double),1, fp) ;
	fread(&DTi,  sizeof(double),1, fp) ;
	fread(&DTr,  sizeof(int),   1, fp) ;

	ZSLOOP(0,NR-1,0,NH-1) PLOOP fread(&p[i][j][k],    sizeof(double), NPR, fp) ;

	return ;
}
