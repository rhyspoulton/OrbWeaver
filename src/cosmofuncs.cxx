
#include "orbweaver.h"

double Integrand(double a,void *p){
	double H0 = 100*Cosmo.h;
	double HT = (3.08568025e+19/(H0 * 31556926))/1e+09;
	return (a*HT)/sqrt(Cosmo.omegaR + Cosmo.omegaM*a + Cosmo.omegaK*(a*a) + Cosmo.omegaL*(a*a*a*a));
}

double *GenerateUniAges(vector<double> &scalefactors){

	//Find the number of points need to intergrate to
	int numpoints = scalefactors.size();

	//Setup the integrand function
	gsl_function F;
	F.function = &Integrand;
	double epsabs=1e-4;
	double epsrel=1e-4;

	//Setup the return values
	double *results;
	results = new double[numpoints];
	double error;
  	size_t neval;

  	//Now do the intergration using a non-adaptive Gauss-Kronrod integration
  	for(int i=0;i<numpoints;i++)
		gsl_integration_qng(&F,0,scalefactors[i],epsabs,epsrel,&results[i],&error,&neval);

	return results;
}