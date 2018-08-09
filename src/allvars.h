
#ifndef ALLVARS_H
#define ALLVARS_H

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/timeb.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>


///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>

///if using OpenMP API
#ifdef USEOPENMP
#include <omp.h>
#endif


///if using HDF API
#ifdef USEHDF
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif

///if using ADIOS API
#ifdef USEADIOS
#include "adios.h"
#endif

//Define the amount of fields to read in from the hdf5 file 
#define NHDFFIELDS 14

//Comoving or physical flags
#define COMOVING 0
#define PHYSICAL 1

using namespace std;
using namespace Math;
using namespace NBody;

//Global varible to store the return code when the code terminates
extern int EXIT_CODE;

struct Options
{
	
	// Filenames
	char *fname;

	//The name of the orbweaver configuration
	char *configname;

	//The starting snapshot
	int isnap;

	//The ending snapshot
	int fsnap;

	//Number of snapshots
	int numsnaps;

	// The value which the halo ID snapvalue is offset by
	long long TEMPORALHALOIDVAL;

	//How verbose the code is
	int iverbose;

	//Default options
	Options()
	{
		fname=NULL;
		configname=NULL;
		numsnaps=fsnap-isnap+1;
		iverbose=0;
	}
};

#ifdef USEHDF




//Struct to load the halo into
struct HaloData{

	//The ID of the halo
	unsigned long long id;

	//The descendant of the halo
	unsigned long long descendant;

	//The progenitor of the halo
	unsigned long long progenitor;

	//The ID of the host of this halo
	long long hosthaloid;

	//The ID of the halo which it is currently orbiting
	long long orbitinghaloid; 

	//Virial mass of the halo
	double mvir;

	//Virial radius of the halo
	double rvir;

	//Position of the halo
	double x,y,z;

	//Velocity of the halo
	double vx,vy,vz;

	//Rmax of the halo
	double rmax;

	//Vmax of the halo
	double vmax;

	//Flag to mark if the halo has been process or no
	bool doneflag;
	
};



// Struct to store all the data for the halos
struct SnapData{
	Int_t numhalos;
	double scalefactor;
	double uniage; // Age of the universe at this snapshot
	vector<HaloData> Halo;

	SnapData(){
		numhalos=0;
		scalefactor=0.0;
		Halo.empty();
	};	
	~SnapData(){
		for(int i;i<Halo.size();i++)
			delete &Halo[i];
	};
};

struct OrbitData{

	//The number of orbits the halo is on
	float numorbits;

	//The scalefactor when pass 3 rvir for first time
	float timepass_R3;

	//The scalefactor when pass 2 rvir for first time
	float timepass_R2;

	//The scalefactor when pass 1 rvir for first time
	float timepass_R1;

	//The orbital period in difference in scalefactors
	float orbitperiod;

	//Closest approach for this halo to the halo it is orbiting
	float closestapproach;

	//The mass of the halo when passes into its hosts Rvir 
	float massatinfall;

	//Vmax of the halo when it passes into its hosts Rvir
	float vmaxinfall;

	//Mass loss rate
	float masslossrate;

};

struct UnitsData{
	bool distFlag;
	double length;
	double mass;
	double velocity;
};


struct CosmoData{
	double h;
	double boxsize;
	double omegaM;
	double omegaL;
	double omegaR;
	double omegaK;
};

extern UnitsData Units;
extern CosmoData Cosmo;

struct HDFCatalogNames{

	//Set the name for the cosmology group
	string cosmohdrname;

	//Set the name for the units group
	string unithdrname;

	//Set the base name for the snapshot
	string grpbasename;

	// Store the names of the attributes in the cosmology header
	vector<H5std_string> cosmoattrnames;

	// Store the names of the attributes in the cosmology header
	vector<H5std_string> unitattrnames;

	// Store the names of the attributes that needs to be read in for each snapshot
	vector<H5std_string> snapattrnames;



	// Store the names of the datasets that needs to be read in
	vector<H5std_string> datasetnames;
	// Store their datatype
	vector<PredType> datasettypes;




	HDFCatalogNames(){

		grpbasename = "Snap_%03d";
		cosmohdrname = "/Header/Cosmology";
		unithdrname = "/Header/Units";

		cosmoattrnames.push_back("BoxSize");
		cosmoattrnames.push_back("Hubble_param");

		unitattrnames.push_back("Comoving_or_Physical");
		unitattrnames.push_back("Length_unit_to_kpc");
		unitattrnames.push_back("Mass_unit_to_solarmass");
		unitattrnames.push_back("Velocity_unit_to_kms");

		snapattrnames.push_back("NHalos");
		snapattrnames.push_back("scalefactor");

		datasetnames.push_back("Xc");
		datasetnames.push_back("Yc");
		datasetnames.push_back("Zc");
		datasetnames.push_back("ID");
		datasetnames.push_back("Head");
		datasetnames.push_back("Tail");
		// datasetnames.push_back("hostHaloID");
		datasetnames.push_back("OrbitingHaloID");
		datasetnames.push_back("Mass_200crit");
		datasetnames.push_back("R_200crit");
		datasetnames.push_back("VXc");
		datasetnames.push_back("VYc");
		datasetnames.push_back("VZc");
		datasetnames.push_back("Rmax");
		datasetnames.push_back("Vmax");

		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::STD_U64LE);
		datasettypes.push_back(PredType::STD_U64LE);
		datasettypes.push_back(PredType::STD_U64LE);
		// datasettypes.push_back(PredType::STD_I64LE);
		datasettypes.push_back(PredType::STD_I64LE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);


	};

	
};


#endif //USEHDF5

#endif //ifndef ALLVARS_H






