
#include "decs.h"

/* apply floors to density, internal energy */

void fixup(double (* pv)[NH+4][NPR]) 
{
	int i,j ;
	void pfixup(double *p) ;

	ZSLOOP(-2,NR+1,-2,NH+1) pfixup(pv[i][j]) ;

}

void pfixup(double *pr)
{
	/* floor on density (momentum *not* conserved) */
	if(pr[RHO] < RHOMIN) 
		pr[RHO] = RHOMIN ;

	/* floor on internal energy */
	if(pr[UU] < UUMIN)
		pr[UU] = UUMIN ;
}
