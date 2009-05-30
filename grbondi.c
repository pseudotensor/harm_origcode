
/* 
 * solves eigenvalue problem for spherically symmetric accretion 
 * in the Schwarzschild metric.
 *
 * NB: assumes gam = 5/3 
 *
 * cfg 6-9-01 
 *
 * 6-17-01: removed gam=5/3 restriction
 *  	    rationalized algorithm & functions
 *
 * 6-18-01: added bondi_trace to find solution
 *
 */

#include <stdio.h>
#include <math.h>

/* 
 * bondi_trace
 *
 * Given solution to the eigenvalue problem,
 * find solution with a given edot at a given radius.
 * Useful for tracing out a solution to the bondi problem
 * Use after calling bondi_solve.
 *
 * output: 
 * 	returns ur, the radial four-velocity at radius r
 *
 * input: 
 * 	K      entropy constant (p = K rho^(gam))
 * 	gam    adiabatic index
 * 	edotf  energy accretion rate (from bondi_solve)
 * 	r      current radius
 * 	rs     sonic radius.  Used for selecting correct branch
 * 		of solution
 *	urs    radial four-velocity of solution at sonic radius.
 *
 */

double bondi_trace(double K, double gam, double edotf, 
		double r, double rs, double urs)
{
	double tol,ur,dur ;
	double edot,edot_last ;
	double edot_calc(),dedur_calc() ;

	/* set convergence criterion */
	tol = 1.e-10 ;

	/* set initial guess that will automatically 
	 * select the correct branch of the solution */
	if(r > rs) {
		ur = 0.1*urs*sqrt(rs/r) ;
		dur = -0.01 ;
	}
	else {
		ur = 10.*urs*sqrt(rs/r) ;
		dur = 0.01 ;
	}

	/* step along in ur until root is found */
	edot = edot_calc(r,ur,gam,K) ;
	edot_last = edot ;
	while(fabs(dur) > tol) {
		edot_last = edot ;
		edot = edot_calc(r,ur,gam,K) ;

		/* just passed a zero; turn around and go back */
		if( (edot_last - edotf) * (edot - edotf) < 0. )  
			dur *= -0.5 ; 
		/* headed in the wrong direction; turn around */
		else if( fabs(edot_last - edotf) < fabs(edot - edotf) )
			dur *= -1. ;

		/* take step */
		ur += ur*dur ;
	}
	return(ur) ;
}


/*
 * bond_solve
 *
 *
 * assumes:
 * 	polytropic equation of state
 * 	mass accretion rate = 1. (code actually solves for an
 * 	 excretion flow, but solution is symmetric under change
 * 	 of sign of u^r).
 *
 *
 * output: 
 * 	Rs    sonic radius
 * 	Urs   radial four-velocity of fluid at sonic radius
 * 	Edot  energy accretion rate
 *
 * input: 
 * 	K    entropy constant (p = K rho^(gam))
 * 	gam  adiabatic index 
 *
 * method:
 * 	search for zeros of both d(edot)/dr *and* d(edot)/dur,
 * 	which one can prove corresponds to the sonic point.
 *
 */

void bondi_solve(double K, double gam, double *Rs, 
		double *Urs, double *Edot) 
{
	double tolb ;
	double edot_calc(),dedr_calc(),dedur_calc(),d2edur2_calc(),d2edr2_calc(),d2edrdur_calc() ;
	double dedr,dedr_last,dedur,d2edur2,d2edr2,d2edrdur ;
	double e_r,e_u,ep_r,ep_u,lambda ;
	double rs,urs,norm ;
	int nstep = 0 ;
	int dir,last_dir ;

	/* these are a few magic numbers that depend (insensitively) 
	   on the functions dedr and dedur */
	lambda = 0.02 ;
	tolb = 1.e-10 ;

	/* here's a generic guess */
	rs = 1.e5 ;
	urs = 0.001 ;

	/* move along locus dedur = 0 until a zero of dedr is crossed */
	last_dir = 1 ;
	dir = 1 ;
	dedr = dedr_calc(rs,urs,gam,K) ;
	dedr_last = dedr ;
	while(fabs(lambda) > tolb) {

		/* find the locus of zeros for dedur at current rs via Newton-Raphson */
		dedur = dedur_calc(rs,urs,gam,K) ;
		while(fabs(dedur) > tolb) {
			d2edur2 = d2edur2_calc(rs,urs,gam,K) ;
			urs = urs - dedur/d2edur2 ; 
			dedur = dedur_calc(rs,urs,gam,K) ;
		}

		/* this is the unit vector parallel to the gradient in dedur */
		d2edrdur = d2edrdur_calc(rs,urs,gam,K) ;
		d2edur2 = d2edur2_calc(rs,urs,gam,K) ;
		norm = sqrt(d2edur2*d2edur2 + d2edrdur*d2edrdur) ;
		e_r = d2edrdur/norm ;
		e_u = d2edur2/norm ;

		/* this is the unit vector perp to the gradient in dedur */
		ep_r = -e_u ;
		ep_u = e_r ;

		/* default direction is to step inward */
		if(ep_r > 0.) {
			ep_r *= -1 ;
			ep_u *= -1 ;
		}

		/* set step size */
		dedr_last = dedr ;
		dedr = dedr_calc(rs,urs,gam,K) ;
		if(dedr*dedr_last < 0)
			lambda *= -0.5 ;

		/* check that we're going the right way */
		if(nstep == 0) {
			d2edr2 = d2edr2_calc(rs,urs,gam,K) ;
			if(d2edr2*dedr > 0.)
				lambda *= -1. ;
		}

		/* take a step along locus dedur == 0 */
		rs += lambda*rs*ep_r ;
		urs += lambda*urs*ep_u ;

		/* check for nonsense values */
		if(rs < 2.*(1. + lambda)) {
			fprintf(stderr,"error: too close to event horizon\n") ;
			exit(1) ;
		}
		if(urs < 0.) {
			fprintf(stderr,"error: radial four-velocity changed signs\n") ;
			exit(1) ;
		}

		/*
		fprintf(stderr,"> %g %g %g %g %g\n",rs,urs,dedr,dedur,lambda) ;
		*/
		nstep++ ;
	}

	/* done! */
	*Rs = rs ;
	*Urs = urs ;
	*Edot = edot_calc(rs,urs,gam,K) ;

	fprintf(stderr,"nstep: %d\n",nstep) ;

	return ;
}

/*
 * The following functions were generated from Mathematica.
 * They have not been optimized to maximize numerical accuracy,
 * but they should be.
 *
 */

double edot_calc(double r, double ur, double g, double K)
{

	return(
		((pow(4*M_PI,g) - g*(pow(4*M_PI,g) + 4*K*M_PI*pow(1/(pow(r,2)*ur),-1 + g)))*
       		sqrt(1 - 2/r + pow(ur,2)))/((-1 + g)*pow(4*M_PI,g))
	) ;

}

double dedr_calc(double r,double ur,double g, double K) 
{

	return(
		(pow(16.*M_PI,g) + pow(2.,3 + 2.*g)*pow(g,2)*K*M_PI*
		pow(1/(pow(r,2)*ur),-1. + g)*(-2. + r + r*pow(ur,2)) - 
		g*(pow(16*M_PI,g) + 4.*K*M_PI*pow(1/(pow(r,2)*ur),-1. + g)*
		(-3.*pow(4.,g) + pow(2.,1 + 2*g)*r*(1. + pow(ur,2)))))/
		((-1. + g)*pow(16.*M_PI,g)*pow(r,2)*sqrt(1. - 2./r + pow(ur,2)))
	) ;

}

double dedur_calc(double r, double ur,double g, double K) 
{

	return(
		(pow(4.*M_PI,g)*ur + 4*pow(g,2)*K*M_PI*r*pow(1/(pow(r,2)*ur),g)
		 *(-2 + r + r*pow(ur,2)) - 
		g*(pow(4*M_PI,g)*ur + 4*K*M_PI*r*pow(1/(pow(r,2)*ur),g)*
		(-2 + r + 2*r*pow(ur,2))))/
		((-1 + g)*pow(4*M_PI,g)*sqrt(1 - 2/r + pow(ur,2)))
	) ;
}

double d2edr2_calc(double r, double ur, double g, double K) 
{

	return(
		((16*(-1 + g)*g*K*M_PI*pow(1/(pow(r,2)*ur),g)*ur)/r - 
		8*g*(1 - 3*g + 2*pow(g,2))*K*M_PI*pow(1/(pow(r,2)*ur),g)*ur*
		(1 - 2/r + pow(ur,2)) + 
		((-pow(4*M_PI,g) + g*(pow(4*M_PI,g) + 
		4*K*M_PI*pow(1/(pow(r,2)*ur),-1 + g)))*
		(-3 + 2*r*(1 + pow(ur,2))))/(pow(r,3)*(-2 + r + r*pow(ur,2))))/
		((-1 + g)*pow(4*M_PI,g)*sqrt(1 - 2/r + pow(ur,2)))
	) ;

}

double d2edur2_calc(double r, double ur, double g, double K)
{
	return(
		(8*(-1 + g)*g*K*M_PI*pow(1/(pow(r,2)*ur),-1 + g) - 
		((-2 + r)*(-pow(4*M_PI,g) + 
		g*(pow(4*M_PI,g) + 4*K*M_PI*pow(1/(pow(r,2)*ur),-1 + g))))/
		(-2 + r + r*pow(ur,2)) - 
		(4*(-1 + g)*pow(g,2)*K*M_PI*r*pow(1/(pow(r,2)*ur),g)*
		(-2 + r + r*pow(ur,2)))/ur)/
		((-1 + g)*pow(4*M_PI,g)*sqrt(1 - 2/r + pow(ur,2)))
	) ;

}

double d2edrdur_calc(double r, double ur, double g, double K) 
{
	return(
		(-(((1 + g*(-1 - K*pow(4*M_PI,1 - g)*pow(1/(pow(r,2)*ur),-1 + g)))*ur)/
		((-1 + g)*pow(r,2))) + 
		g*K*pow(4*M_PI,1 - g)*pow(1/(pow(r,2)*ur),g)*
		(1 - 2/r + pow(ur,2)) + 
		pow(2,3 - 2*g)*g*K*pow(M_PI,1 - g)*pow(1/(pow(r,2)*ur),g)*
		pow(ur,2)*(-2 + r + r*pow(ur,2)) - 
		(pow(2,3 - 2*g)*g*K*pow(M_PI,1 - g)*pow(1/(pow(r,2)*ur),g)*
		pow(-2 + r + r*pow(ur,2),2))/r - 
		(pow(2,3 - 2*g)*(-2 + g)*g*K*pow(M_PI,1 - g)*
		pow(1/(pow(r,2)*ur),g)*pow(-2 + r + r*pow(ur,2),2))/r)/
		pow(1 - 2/r + pow(ur,2),1.5)
	) ;
}

