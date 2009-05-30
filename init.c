
/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"

void init()
{
	int i,j,k ;
	void set_arrays(),set_grid() ;
	void bound_prim( double prim[][NH+4][NPR] ) ;
	double r,th,R,sth,cth ;
	double b2_calc(double *p) ;
        void bltoks(double *pr, int i, int j) ;
	double ur,uh,up,u,rho ;
	double ranc(int seed) ;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;
	double rho2 ;

	/* for magnetic field */
	double A[NR+1][NH+1] ;
	double rho_av,rhomax,umax,beta,b2_ij,b2_max,norm,q,beta_act ;

	/* some physics parameters */
	gam = 5./3. ;

	/* disk parameters (use fishbone.m to select new solutions) */
	/*
	a = 0.5 ;
	l = 4.6 ;
	rin = 5.0 ;
	kappa = 1.e-3 ;
	beta = 40.*4.e2 ;
	*/
	/*
	*/

	
        a = 0. ;
        l = 5. ;
        rin = 5.7 ;
        kappa = 1.e-3 ;
        beta = 4.e2 ;
	

	/* numerical parameters */
	fail = 0 ;
	lim = MC ;	        /* used montonized central-difference limiter */
	cour = ranc(1) ; 	/* start with known random seed */
	cour = 0.4 ;
	dt = 0.0001 ;
	Rin = 1.7 ;
	Rout = 20. ;

	/* set up arrays,details */
	dr = (Rout - Rin)/NR ;	
	dh = M_PI/NH ;
	dV = dr*dh ;
	t = 0. ;
	fprintf(stderr,"Rmin = %10.5g\n",Rin-2*dr) ;
	set_arrays() ;
	set_grid() ;

	/* output parameters */
	tf = 2000.0 ;
	DTd = 16. ;	/* dumping frequency, in units of M */
	DTl = 4. ;	/* logfile frequency, in units of M */
	DTi = 4. ; 	/* image file frequ., in units of M */
	DTr = 32 ; 	/* restart file frequ., in timesteps */

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,NR,0,NH) A[i][j] = 0. ;
	ZSLOOP(0,NR-1,0,NH-1) {
		r = (i + 0.5)*dr + Rin ;
		th = (j + 0.5)*dh ;
		sth = sin(th) ;
		cth = cos(th) ;
		R = r*sth ;

		/* calculate lnh */
		DD = r*r - 2.*r + a*a ;
		AA = (r*r + a*a)*(r*r + a*a) - DD*a*a*sth*sth ;
		SS = r*r + a*a*cth*cth ;

		thin = M_PI/2. ;
		sthin = sin(thin) ;
		cthin = cos(thin) ;
		DDin = rin*rin - 2.*rin + a*a ;
		AAin = (rin*rin + a*a)*(rin*rin + a*a) 
			- DDin*a*a*sthin*sthin ;
		SSin = rin*rin + a*a*cthin*cthin ;

		if(r >= rin) {
			lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
				(AA*sth*AA*sth)))/(SS*DD/AA)) 
				- 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
				- 2.*a*r*l/AA 
				- (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
				(AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin)) 
				- 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
					(AAin*AAin*sthin*sthin)) 
				- 2.*a*rin*l/AAin ) ;
		}


		/* regions outside torus */
		if(lnh < 0. || r < rin) {
			rho = RHOMIN ;
                        u = UUMIN ;
                        ur = -1./sqrt(r) ;
                        uh = 0. ;
			rho2 = r*r + a*a*cth*cth ;
			up = a*ur/rho2 ;	/* => zero angular momentum */

			p[i][j][RHO] = rho ;
			p[i][j][UU] = u ;
			p[i][j][UR] = ur ;
			p[i][j][UH] = uh ;
			p[i][j][UP] = up ;
		}
		/* region inside magnetized torus; u^i is calculated in
		 * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
		 * so it needs to be transformed at the end */
		else { 
			hm1 = exp(lnh) - 1. ;
			rho = pow(hm1*(gam - 1.)/(kappa*gam),
						1./(gam - 1.)) ; 
			u = kappa*pow(rho,gam)/(gam - 1.) ;
			ur = 0. ;
			uh = 0. ;

			/* calculate u^phi */
			expm2chi = SS*SS*DD/(AA*AA*sth*sth) ;
			up1 = sqrt((-1. + sqrt(1. + 4.*l*l*expm2chi))/2.) ;
			up = 2.*a*r*sqrt(1. + up1*up1)/sqrt(AA*SS*DD) +
				sqrt(SS/AA)*up1/sth ;

			p[i][j][RHO] = rho*(1. + 0.e-2*(ranc(0)-0.5)) ;
			if(rho > rhomax) rhomax = rho ;
			p[i][j][UU] = u ;
			if(u > umax && r > rin) umax = u ;
			p[i][j][UR] = ur*(1. + 0.e-2*(ranc(0)-0.5)) ;
			p[i][j][UH] = uh*(1. + 0.e-2*(ranc(0)-0.5)) ;
			p[i][j][UP] = up ;

        		/* transform four-velocity to kerr-schild */
        		bltoks(p[i][j],i,j) ;
		}

		p[i][j][UR] += 1.e-3*(ranc(0)-0.5) ;
		p[i][j][UH] += 1.e-3*(ranc(0)-0.5) ;
		p[i][j][UP] += 1.e-3*(ranc(0)-0.5) ;

		p[i][j][BR] = 0. ;
		p[i][j][BH] = 0. ;
		p[i][j][BP] = 0. ;
	}
	fprintf(stderr,"rhomax: %g\n",rhomax) ;
	ZSLOOP(0,NR-1,0,NH-1) {
		p[i][j][RHO] /= rhomax ;
		p[i][j][UU]  /= rhomax ;
	}
	rhomax = 1. ;
	umax /= rhomax ;
	bound_prim(p) ;

	/* first find corner-centered vector potential */
	ZSLOOP(0,NR,0,NH) {

		r = (i + 0.5)*dr + Rin ;
		th = (j + 0.5)*dh ;
		sth = sin(th) ;

		R = r*sth ;

		GSET(i,j,CORN) ;
		A[i][j] = (R > rin) ? g*(R - rin) : 0. ;
	}

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	b2_max = 0. ;
	ZLOOP {
		GSET(i,j,CENT) 
		p[i][j][BR] =  (A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1])/(2.*dh*g) ;
		p[i][j][BH] = -(A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/(2.*dr*g) ;
		p[i][j][BP] = 0. ;

		b2_ij = b2_calc(p[i][j]) ;
		if(b2_ij > b2_max) b2_max = b2_ij ;
	}
	fprintf(stderr,"initial b2_max: %g\n",b2_max) ;

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.)*umax/(0.5*b2_max) ;
	fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
	norm = sqrt(beta_act/beta) ;
	b2_max = 0. ;
	ZLOOP {
		p[i][j][BR] *= norm ;
		p[i][j][BH] *= norm ;

		GSET(i,j,CENT) 
		b2_ij = b2_calc(p[i][j]) ;
		if(b2_ij > b2_max) b2_max = b2_ij ;
	}
	beta_act = (gam - 1.)*umax/(0.5*b2_max) ;
	fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;

	/* enforce boundary conditions */
	bound_prim(p) ;

}
