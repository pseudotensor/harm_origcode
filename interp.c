
#include "decs.h"

double slope_lim(double y1,double y2,double y3) 
{
	double Dqm,Dqp,Dqc,s ;

	/* woodward, or monotonized central, slope limiter */
	if(lim == MC) {
		Dqm = 2.*(y2 - y1) ;
		Dqp = 2.*(y3 - y2) ;
		Dqc = 0.5*(y3 - y1) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else {
			if(fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return(Dqm) ;
			else if(fabs(Dqp) < fabs(Dqc))
				return(Dqp) ;
			else
				return(Dqc) ;
		}
	}
	/* van leer slope limiter */
	else if(lim == VANL) {
		Dqm = (y2 - y1) ;
		Dqp = (y3 - y2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else
			return(2.*s/(Dqm+Dqp)) ;
	}
	/* minmod slope limiter (crude but robust) */
	else if(lim == MINM) {
		Dqm = (y2 - y1) ;
		Dqp = (y3 - y2) ;
		s = Dqm*Dqp ;
		if(s <= 0.) return 0. ;
		else if(fabs(Dqm) < fabs(Dqp)) return Dqm ;
		else return Dqp ;
	}
}

