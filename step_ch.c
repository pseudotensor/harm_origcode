
#include "decs.h"

/* Toth/Rusanov/Yee  method */
#define FAC	0.5
#define SAFE	1.5

void step_ch()
{
	int i,j,k ;

	double ndt,dtr,dth,ndtr,ndth,cr,ch ;
	double tmp[NPR] ;
	static double tfail ;
	static double Uh[NPR],U[NPR],dU[NPR] ;
	static double pir[NPR],plr[NPR],prr[NPR],Flr[NPR],Frr[NPR],Ulr[NPR],Urr[NPR] ;
	static double pih[NPR],plh[NPR],prh[NPR],Flh[NPR],Frh[NPR],Ulh[NPR],Urh[NPR] ;

	static double Fr_ct[NR+1][NH+1][NPR] ;
	static double Fh_ct[NR+1][NH+1][NPR] ;

	void primtoU(double *pb, double *Ub) ;
	void Utoprim(double *Ua, double *pa) ;
	void primtorflux(double *pa, double *Fa) ;
	void primtohflux(double *pa, double *Fa) ;
	void source(double *pa, double *Ua) ;
	void bound_prim(double (* var)[NH+4][NPR]) ;
	void fixup(double (* var)[NH+4][NPR]) ;
	void pfixup(double *pr) ;

	double vrchar(double *pa) ;
	double vhchar(double *pa) ;
	double slope_lim(double y1, double y2, double y3) ;

	fprintf(stderr,"0") ;

	/** evaluate Woodward slopes of primitive variables **/
	ZSLOOP(-1,NR,-1,NH) PLOOP {
			dqr[i][j][k] = slope_lim(p[i-1][j][k],p[i][j][k],p[i+1][j][k]) ;
	}
	ZSLOOP(-1,NR,-1,NH) PLOOP {
			dqh[i][j][k] = slope_lim(p[i][j-1][k],p[i][j][k],p[i][j+1][k]) ;
	}

	fprintf(stderr,"1") ;
	/** perform Hancock predictor half-step **/
	ZLOOP {
		PLOOP {
			plr[k] = p[i][j][k] - 0.5*dqr[i][j][k] ;
			prr[k] = p[i][j][k] + 0.5*dqr[i][j][k] ;
			plh[k] = p[i][j][k] - 0.5*dqh[i][j][k] ;
			prh[k] = p[i][j][k] + 0.5*dqh[i][j][k] ;
		}

		GSET(i,j,RFACE) 
		primtorflux(plr,Flr) ;
		GSET(i+1,j,RFACE)
		primtorflux(prr,Frr) ;

		GSET(i,j,HFACE)
		primtohflux(plh,Flh) ;
		GSET(i,j+1,HFACE)
		primtohflux(prh,Frh) ;

		/* find source terms */
		GSET(i,j,CENT)
		source(p[i][j],dU) ;

		primtoU(p[i][j],Uh) ;
		/*
		fprintf(stderr,"\n Fluxes\n") ;
		PLOOP fprintf(stderr,"%d %d %d %g %g %g %g\n",
				i,j,k, Frr[k], Flr[k], Frh[k], Flh[k]) ;
		fprintf(stderr,"\n Sources\n") ;
		PLOOP fprintf(stderr,"%d %d %d %g %g\n",
				i,j,k, Uh[k], dU[k]) ;
		*/

		PLOOP {
			Uh[k] += 0.5*dt*( 
				- (Frr[k] - Flr[k])/dr 
				- (Frh[k] - Flh[k])/dh 
				+ dU[k] 
				) ;
		}

		PLOOP ph[i][j][k] = p[i][j][k] ;	/* needed for Utoprim */
		Utoprim(Uh,ph[i][j]) ;
	}
	bound_prim(ph) ;
	fixup(ph) ;

	fprintf(stderr,"2") ;
	/** evaluate fluxes at t + 0.5*dt **/
	/* radial fluxes */
	ZSLOOP(0,NR,-1,NH) {
		PLOOP pir[k] = 0.5*(ph[i-1][j][k] + 0.5*dqr[i-1][j][k] +
			   	      ph[i][j][k]   - 0.5*dqr[i][j][k]) ;
		pfixup(pir) ;

		GSET(i,j,RFACE)
		primtorflux(pir, Fr[i][j]) ;
	}
	/* theta fluxes */
	ZSLOOP(-1,NR,0,NH) {
		PLOOP pih[k] = 0.5*(ph[i][j-1][k] + 0.5*dqh[i][j-1][k] +
			   	      ph[i][j][k]   - 0.5*dqh[i][j][k]) ;
		pfixup(pih) ;

		GSET(i,j,HFACE)
		primtohflux(pih, Fh[i][j]) ;
	}

	/** evaluate diffusive (Lax-Friedrichs like) correction to flux, at halfstep **/
	fprintf(stderr,"3") ;
	ndtr = 1.e9 ;
	/* radial fluxes */
	ZSLOOP(0,NR,-1,NH) {
		PLOOP {
			plr[k] = ph[i-1][j][k] + 0.5*dqr[i-1][j][k] ;
			prr[k] = ph[i][j][k]   - 0.5*dqr[i][j][k] ;
			pir[k] = 0.5*(plr[k] + prr[k]) ;
		}
		pfixup(plr) ;
		pfixup(prr) ;
		pfixup(pir) ;

		GSET(i,j,RFACE) 
		primtoU(plr,Ulr) ;
		primtoU(prr,Urr) ;
		cr = vrchar(pir) ;

		/* evaluate restriction on timestep */
		dtr = cour*dr/cr ;
		if(dtr < ndtr) ndtr = dtr ;

		PLOOP Fr[i][j][k] += FAC*cr*(Ulr[k] - Urr[k]) ;
	}
	ndth = 1.e9 ;
	/* theta fluxes */
	ZSLOOP(-1,NR,0,NH) {
		PLOOP {
			plh[k] = ph[i][j-1][k] + 0.5*dqh[i][j-1][k] ;
			prh[k] = ph[i][j][k]   - 0.5*dqh[i][j][k] ;
			pih[k] = 0.5*(plh[k] + prh[k]) ;
		}
		pfixup(plh) ;
		pfixup(prh) ;
		pfixup(pih) ;

		GSET(i,j,HFACE) 
		primtoU(plh,Ulh) ;
		primtoU(prh,Urh) ;
		ch = vhchar(pih) ;

		/* evaluate restriction on timestep */
		dth = cour*dh/ch ;
		if(dth < ndth) ndth = dth ;

		PLOOP Fh[i][j][k] += FAC*ch*(Ulh[k] - Urh[k]) ;
	}

	/** here's where the next timestep is set **/
	ndt = DMIN(ndtr,ndth) ;

#if 1
	/** replace B-field fluxes with flux-ct symmetrized fluxes **/
	fprintf(stderr,"4") ;
	/* radial fluxes */
	ZSLOOP(0,NR,0,NH-1) {
		Fr_ct[i][j][BR] = 0. ;
		Fr_ct[i][j][BH] = 0.125*( 
			2.*Fr[i][j][BH]
			+ Fr[i][j+1][BH]
			+ Fr[i][j-1][BH]
			- Fh[i][j][BR]
			- Fh[i-1][j][BR]
			- Fh[i][j+1][BR]
			- Fh[i-1][j+1][BR] ) ;
	}
	/* theta fluxes */
	ZSLOOP(0,NR-1,0,NH) {
		Fh_ct[i][j][BR] = 0.125*( 
			2.*Fh[i][j][BR]
			+ Fh[i+1][j][BR]
			+ Fh[i-1][j][BR]
			- Fr[i][j][BH]
			- Fr[i][j-1][BH]
			- Fr[i+1][j][BH]
			- Fr[i+1][j-1][BH] ) ;
		Fh_ct[i][j][BH] = 0. ;
	}
	ZSLOOP(0,NR,0,NH-1) {
		Fr[i][j][BR] = Fr_ct[i][j][BR] ;
		Fr[i][j][BH] = Fr_ct[i][j][BH] ;
	}
	ZSLOOP(0,NR-1,0,NH) {
		Fh[i][j][BR] = Fh_ct[i][j][BR] ;
		Fh[i][j][BH] = Fh_ct[i][j][BH] ;
	}
#endif

	/* evaluate some diagnostics */
	mdot = edot = ldot = 0. ;
	for(j=0;j<NH;j++) {
		mdot += Fr[0][j][RHO]*2.*M_PI*dh ;
		edot -= Fr[0][j][UU] *2.*M_PI*dh ;
		ldot += Fr[0][j][UP] *2.*M_PI*dh ;
	}

	fprintf(stderr,"5") ;
	ZLOOP {
		/* find time-centered source terms */
		GSET(i,j,CENT)
		source(ph[i][j],dU) ;

		/* evolve conserved quantities */
		primtoU(p[i][j], U) ;
		PLOOP {  
			U[k] += dt*(
			- (Fr[i+1][j][k] - Fr[i][j][k])/dr 
			- (Fh[i][j+1][k] - Fh[i][j][k])/dh
			+  dU[k]
			);
		}

		/* recover new primitive variables; save the
		 * old primitive variables in ph in case the
		 * timestep fails */
		PLOOP {
			tmp[k] = p[i][j][k] ;
			p[i][j][k] = ph[i][j][k] ; 	
			ph[i][j][k] = tmp[k] ;
		}
		Utoprim(U, p[i][j]) ;
	}

	/* apply boundary conditions */
	fprintf(stderr,"6\n") ;
	bound_prim(p) ;
	fixup(p) ;

	/* if timestep fails, go back and do it with a smaller timestep */
	if(fail) {
		fail = 0 ;
		tfail = t ;
		lim = MINM ;	/* use more robust minmod limiter; this
				   should probably be removed */
		dt *= 0.25 ;	/* reduce timestep */
		ZLOOP PLOOP p[i][j][k] = ph[i][j][k] ;
	}
	else {
		/* increment time */
		t += dt ;

		/* set next timestep */
		if(ndt > SAFE*dt) ndt = SAFE*dt ;
		dt = ndt ;
		if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */
	}
	if(t - tfail > 1.) lim = MC ;	/* reset limiter if it's been long enough 
					   since a failure */

	/* done! */
}

