
#include "decs.h"
#include "defs.h"

int main(int argc,char *argv[])
{
	void init() ;
	void timestep() ;
	void step_ch() ;
	void diag(int call_code) ;
	void restart_write() ;
	int restart_init() ;

	double tnext ;
	int nstep ;

	/* perform initializations */
	/*
	if(!restart_init()) init() ;
	*/
	init() ;

	/* do initial diagnostics */
	diag(0) ;

	tnext = t ;
	nstep = 0 ;

	while(t < tf) {

		fprintf(stderr,"%10.5g %10.5g %8d\n",t,dt,nstep) ;

		/* step variables forward in time */
		step_ch() ;

		/* perform diagnostics */
		if(t >= tnext) {
			diag(1) ;
			tnext += DTl ;
		}

		/* restart dump */
		if(nstep%DTr == 0) restart_write() ;

		nstep++ ;
	}
	fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*NR*NH) ;

	/* do final diagnostics */
	diag(2) ;

	return(0) ;
}
