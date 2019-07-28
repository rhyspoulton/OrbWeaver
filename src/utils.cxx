#include "orbweaver.h"


vector<double> computeAngles(vector<double> prevpos, OrbitData orbitdata){

	/*
	
	Function to calculate the Euler angles of orbit

	*/


	double x[3],y[3],z[3], mag;
	vector<double> angles(3);

	//Get the x-axis which is along the semi-major axis
	x[0] = orbitdata.xrel - prevpos[0];
	x[1] = orbitdata.yrel - prevpos[1];
	x[2] = orbitdata.zrel - prevpos[2];

	//The z-axis which is the angular momentum
	z[0] = orbitdata.lxrel_ave;
	z[1] = orbitdata.lyrel_ave;
	z[2] = orbitdata.lzrel_ave;

	//Find the cross product of these to to find the y-axis
	y[0] = (x[1] * z[2]) - (x[2] * z[1]);
	y[1] = -((x[0] * z[2]) - (x[2] * z[0]));
	y[2] = (x[0] * z[1]) - (x[1] * z[0]);


	//Now we can find the x vector that is perpendicular to the y and z vectors
	x[0] = (y[1] * z[2]) - (y[2] * z[1]);
	x[1] = -((y[0] * z[2]) - (y[2] * z[0]));
	x[2] = (y[0] * z[1]) - (y[1] * z[0]);

	//Normalize the vectors
	mag = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	x[0] /= mag;
	x[1] /= mag;
	x[2] /= mag;
	mag = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	y[0] /= mag;
	y[1] /= mag;
	y[2] /= mag;
	mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
	z[0] /= mag;
	z[1] /= mag;
	z[2] /= mag;

	//Compute the Euler angles
	angles[0] = acos(-z[1]/sqrt(1-z[2]*z[2]));
	angles[1] = acos(z[2]);
	angles[2] = acos(y[2]/sqrt(1-z[2]*z[2]));

	return angles;
}

double Menc(double r, double Mvir, double c, double Rvir){
	
	//Calculated the scale radius
	double rs =Rvir/c;

	//Calulate the enclosed mass
	return (Mvir/(log(1+c)-c/(1+c)))*(log(1+r/rs)-r/rs/(1+r/rs));
}