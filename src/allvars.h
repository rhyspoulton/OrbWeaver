
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


//Define the output chunksize for the datasets
#define HDFOUTCHUNKSIZE 8192

//Define the amount of fields to read in from the hdf5 file 
#define NHDFFIELDS 19
#define NHDFFIELDSOUT 44

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

	//The original ID of the halo
	unsigned long long origid;

	//The descendant of the halo
	unsigned long long descendant;

	//The progenitor of the halo
	unsigned long long progenitor;

	//The ID of the host of this halo
	long long hosthaloid;

	//The ID of the halo which it is currently orbiting
	long long orbitinghaloid; 

	//Virial mass of the halo
	double mass;

	//Virial radius of the halo
	double rvir;

	//Position of the halo
	double x,y,z;

	//Velocity of the halo
	double vx,vy,vz;

	//Angular momentum of the halo
	double lx,ly,lz;

	//Rmax of the halo
	double rmax;

	//Vmax of the halo
	double vmax;

	//Concentration of the halo
	double cnfw;

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

	//The orbit number
	int orbitID;

	//The ID of the halo
	unsigned long long haloID;

	//The ID of its host halo
	unsigned long long hosthaloID;

	//Type of datapoint this is 4 = apocenter, 3 = crossing 3x host's rvir, 2 = crossing 2x host's rvir, 1 = crossing 1x host's rvir, 0 = pericenter
	float entrytype;

	//The number of orbits the halo has undergone
	float numorbits;

	//The orbital period
	float orbitperiod;

	//Closest approach for this halo to the halo it is orbiting
	float closestapproach;

	//The orbital eccentricity
	float orbitecc;

	//The orbital eccentricity from ratio of passage distances
	float orbiteccratio;

	//The peri-centric distance
	float rperi;

	//The apo-centric distance
	float rapo;

	//Mass loss rate
	float masslossrate;

	//The longitude of the ascending node with respect to the intial orbital plane
	float longascnode;

	//The inclination with respect to the intial orbital plane
	float inc;

	//The argument of pariapsis with respect to the intial orbital plane
	float argpariap;

	//How aligned the average orbital angular momentum is with the average host's angular momentum
	float hostalignment;

	//The scalefactor
	float scalefactor;

	//The age of the universe
	float uniage;

	//The x position
	float x;

	//The y position
	float y;

	//The z position
	float z;

	//The x velocity
	float vx;

	//The y velocity
	float vy;

	//The z velocity
	float vz;

	//The viral mass of the orbiting halo
	float mass;

	//The vmax of the orbiting halo
	float vmax;

	//The peak vmax up to this point
	float vmaxpeak;

	//The rmax of the orbiting halo
	float rmax;

	//The concentration of the orbiting halo
	float cnfw;

	//The tangential velocity of the halo with respect to its host
	float vtan;

	//The relative x position to the host
	float xrel;

	//The relative y position to the host
	float yrel;

	//The relative z position to the host
	float zrel;

	//The relative x velocity to the host
	float vxrel;

	//The relative y velocity to the host
	float vyrel;

	//The relative z velocity to the host
	float vzrel;

	//The average orbital angular momentum in x-direction since last passage
	float lxrel;

	//The average orbital angular momentum in y-direction since last passage
	float lyrel;

	//The average orbital angular momentum in z-direction since last passage
	float lzrel;

	//The viral radius of the host
	float rvirhost;

	//The viral mass of the host halo
	float masshost;

	//The vmax of the host halo
	float vmaxhost;

	//The rmax of the host halo
	float rmaxhost;

	//The concentration of the host halo
	float cnfwhost;
};

struct OrbitProps{
	//Flag to keep track if the halo is on a bound orbit
	bool orbitingflag;

	//Value to keep track of the time the halo started to orbit
	int prevpassagesnap;

	//Value to keep track of the time of the previous apo/peri-centric pasage
	double prevpassagetime;

	//Use to calculate the average mass loss rate
	double masslossrate;

	double mu;

	// Store the angular momentum vectors so the average can be calculated
	double lx;
	double ly;
	double lz;
	double hostlx;
	double hostly;
	double hostlz;
	double ltot;

	double E;

	OrbitProps(){
		orbitingflag = false;
		prevpassagesnap = 0;
		prevpassagetime = 0.0;
		mu=0.0;
		lx=0.0;
		ly=0.0;
		lz=0.0;
		ltot=0.0;
		E=0.0;

	};
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
	double G;
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
		cosmoattrnames.push_back("Gravity");
		cosmoattrnames.push_back("Omega_Lambda");
		cosmoattrnames.push_back("Omega_m");
		cosmoattrnames.push_back("Omega_r");
		cosmoattrnames.push_back("Omega_k");

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
		datasetnames.push_back("origID");
		datasetnames.push_back("Head");
		datasetnames.push_back("Tail");
		// datasetnames.push_back("hostHaloID");
		datasetnames.push_back("OrbitingHaloID");
		datasetnames.push_back("Mass_tot");
		datasetnames.push_back("R_200crit");
		datasetnames.push_back("VXc");
		datasetnames.push_back("VYc");
		datasetnames.push_back("VZc");
		datasetnames.push_back("Lx");
		datasetnames.push_back("Ly");
		datasetnames.push_back("Lz");
		datasetnames.push_back("Rmax");
		datasetnames.push_back("Vmax");
		datasetnames.push_back("cNFW");

		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::STD_U64LE);
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
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);
		datasettypes.push_back(PredType::NATIVE_DOUBLE);

	};
};


struct HDFOutputNames{

	// Store the names of the datasets
	vector<H5std_string> datasetnames;
	// Store their datatype
	vector<PredType> datasettypes;

	HDFOutputNames(){

		datasetnames.push_back("orbitID");
		datasettypes.push_back(PredType::STD_I32LE);
		datasetnames.push_back("haloID");
		datasettypes.push_back(PredType::STD_U64LE);
		datasetnames.push_back("hosthaloID");
		datasettypes.push_back(PredType::STD_U64LE);
		datasetnames.push_back("entrytype");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("numorbits");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitalperiod");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("closestapproach");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitecc");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbiteccratio");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rperi");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rapo");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("masslossrate");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LongAscNode");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Inclination");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("ArgPariap");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("HostAlignment");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("scalefactor");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("uniage");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("X");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Y");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Z");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VX");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VY");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VZ");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Mass");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmax");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmaxpeak");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rmax");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("cNFW");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vtan");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Xrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Yrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Zrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VXrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VYrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VZrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LXrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LYrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LZrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("R_200crit_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Mass_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmax_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rmax_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("cNFW_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
	};
};


#endif //USEHDF5

#endif //ifndef ALLVARS_H






