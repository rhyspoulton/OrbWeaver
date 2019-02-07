
#include "orbweaver.h"

double Integrand(double a,void *p){
	double H0 = 100*Cosmo.h;
	double HT = (3.08568025e+19/(H0 * 31556926))/1e+09;
	return (a*HT)/sqrt(Cosmo.omegaR + Cosmo.omegaM*a + Cosmo.omegaK*(a*a) + Cosmo.omegaL*(a*a*a*a));
}

// Function to find the age of the universe given a value of the scalefactor
double GetUniverseAge(double scalefactor){

	//Setup the integrand function
	gsl_function F;
	F.function = &Integrand;
	double epsabs=1e-4;
	double epsrel=1e-4;

	//Setup the return values
	double error, result;
  	size_t neval;

  	//Now do the intergration using a non-adaptive Gauss-Kronrod integration
	gsl_integration_qng(&F,0,scalefactor,epsabs,epsrel,&result,&error,&neval);

	return result;
}

double GetScaleFactor(double uniage){
	/*
	Function to find the scalefactor at a given uniage, the error for this is typically less than 0.1%:
	https://physics.stackexchange.com/questions/417609/explicit-equation-for-scale-factor-how-long-valid
	*/

	double H0 = 100*Cosmo.h;
	double HT = (3.08568025e+19/(H0 * 31556926))/1e+09;
	return pow( Cosmo.omegaM/Cosmo.omegaL, 1.0/3.0) * pow(sinh(( (3.0*sqrt(Cosmo.omegaL)) / (2.0*HT) ) * uniage) , 2.0/3.0);
}
