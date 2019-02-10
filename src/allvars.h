
#ifndef ALLVARS_H
#define ALLVARS_H

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/timeb.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multifit.h>


///if using HDF API
#include "H5Cpp.h"


#define _USE_MATH_DEFINES

//Define the output chunksize for the datasets
#define HDFOUTCHUNKSIZE 8192

//Define the amount of fields to read in from the hdf5 file 
#define NHDFFIELDSIN 21
#define NHDFFIELDSOUT 50

//Comoving or physical flags
#define COMOVING 0
#define PHYSICAL 1

using namespace std;
using namespace H5;

//Global varible to store the return code when the code terminates
extern int EXIT_CODE;

struct Options
{
	
	// Filenames
	char *fname;


	//Base name for the output
	char *outputbasename;

	//Number of snapshots
	int numsnaps;

	//The fraction of rvir host that a datapoint is created
	float fracrvircross;

	//The value which the halo ID snapvalue is offset by
	unsigned long long TEMPORALHALOIDVAL;

	//How verbose the code is
	int iverbose;

	//Default options
	Options()
	{
		fname=NULL;
		outputbasename=NULL;
		numsnaps=0;
		fracrvircross=0;
		TEMPORALHALOIDVAL=0;
		iverbose=0;
	}
};




//Struct to load the halo into
struct HaloData{

	//The ID of the halo
	unsigned long long id;

	//The original ID of the halo
	long long origid;

	//The descendant of the halo
	long long descendant;

	//The progenitor of the halo
	long long progenitor;

	//If te halo is a subhalo
	bool hostFlag;

	//The ID of the halo which it is currently orbiting
	long long orbitinghaloid; 

	//The number of particles in this halo
	unsigned long long npart;

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

//Struct to contain the interpolation routine data
struct SplineFuncs{

	//x
	gsl_spline *x;
	gsl_interp_accel *xacc;

	//y
	gsl_spline *y;
	gsl_interp_accel *yacc;

	//z
	gsl_spline *z;
	gsl_interp_accel *zacc;

	//vx
	gsl_spline *vx;
	gsl_interp_accel *vxacc;

	//vy
	gsl_spline *vy;
	gsl_interp_accel *vyacc;

	//vz
	gsl_spline *vz;
	gsl_interp_accel *vzacc;

	SplineFuncs(int nhalo){
		//Lets allocate the memory for the splines
		xacc = gsl_interp_accel_alloc();
		yacc = gsl_interp_accel_alloc();
		zacc = gsl_interp_accel_alloc();
		vxacc = gsl_interp_accel_alloc();
		vyacc = gsl_interp_accel_alloc();
		vzacc = gsl_interp_accel_alloc();

		x=gsl_spline_alloc (gsl_interp_cspline, nhalo);
		y=gsl_spline_alloc (gsl_interp_cspline, nhalo);
		z=gsl_spline_alloc (gsl_interp_cspline, nhalo);
		vx=gsl_spline_alloc (gsl_interp_cspline, nhalo);
		vy=gsl_spline_alloc (gsl_interp_cspline, nhalo);
		vz=gsl_spline_alloc (gsl_interp_cspline, nhalo);
	}

	~SplineFuncs(){
		gsl_interp_accel_free (xacc);
		gsl_interp_accel_free (yacc);
		gsl_interp_accel_free (zacc);
		gsl_interp_accel_free (vxacc);
		gsl_interp_accel_free (vyacc);
		gsl_interp_accel_free (vzacc);
		gsl_spline_free (x);
		gsl_spline_free (y);
		gsl_spline_free (z);
		gsl_spline_free (vx);
		gsl_spline_free (vy);
		gsl_spline_free (vz);
	}
};



// Struct to store all the data for the halos
struct SnapData{
	unsigned long long numhalos;
	double scalefactor;
	double uniage; // Age of the universe at this snapshot
	vector<HaloData> Halo;

	SnapData(){
		numhalos=0;
		scalefactor=0.0;
		Halo.empty();
	};	
	~SnapData(){
		Halo.clear();
		Halo.shrink_to_fit();
	};
};

struct OrbitData{

	//The orbit number
	unsigned long long orbitID;

	//The haloID in the orbit catalog
	unsigned long long orbithaloID;

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

	//Average orbital Energy
	float orbitalenergy;

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

	//The angle moved through since last orbit
	float phi;

	//How aligned the average orbital angular momentum is with the average host's angular momentum
	float hostalignment;

	//The time since crossing 1 rvir to the time the halo is lost
	float mergertimescale;

	//The scalefactor
	float scalefactor;

	//The age of the universe
	float uniage;

	//Flag to identify if this galaxy would of merged in a hydrodynamical simulation
	bool mergedflag;

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

	//The number of particles in the orbiting halo
	unsigned long long npart;

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

	//Boolean flag if this halo is a host halo (top of it spatial herachy)
	bool hostFlag;

	//The radial velocity of the halo with respect to its host
	float vrad;

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

	//The number of particles in the host halo
	unsigned long long nparthost;

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

	//Boolean flag if the host halo is a host halo (top of it spatial herachy)
	bool hostFlaghost;
};

struct OrbitProps{
	//Flag to keep track if the halo is on a bound orbit
	bool orbitingflag;

	//Store the number of orbits
	float numorbits;

	//Store the previous crossing point entry type
	float prevcrossingentrytype;

	//Store the previous passage entry type
	float prevpassageentrytype;

	//Store the snapshot which the previous crossing point happened
	int prevcrossingsnap;

	//Store the snapshot which the previous passage happened
	int prevpassagesnap;

	//Value to keep track of the time of the previous apo/peri-centric pasage
	double prevpassagetime;

	//Store the index of the previous passage
	int prevpassageindex;

	//Store the median orbital period, used to remove small oscilations
	double medianperiod;

	//Use to calculate the average mass loss rate
	double masslossrate;

	//The time the halo crossed 1Rvir
	double crossrvirtime;

	//The time merged if before the halo is lost
	double mergertime;

	//Store the reference angles
	double *refangles;

	//Store the radial vector for the previous position
	double prevpassagepos[3];

	//The closest approach; that the halo had to its host
	double closestapproach;

	// Store the reduced mass,angular momentum vectors for the orbiting halo 
	// and its host, total orbital angular momentum, orbital energy and 
	// gravitational velocity so the average can be calculated
	double mu;
	double lx;
	double ly;
	double lz;
	double hostlx;
	double hostly;
	double hostlz;
	double ltot;
	double E;
	double phi;

	OrbitProps(){
		orbitingflag = false;
        refangles = NULL;
		numorbits=0.0;
		prevcrossingentrytype=0.0;
		prevpassagesnap = 0;
		prevpassagetime = 0.0;
		prevpassageindex=0;
		crossrvirtime=0.0;
		prevpassageindex=0;
		mergertime=0.0;
		masslossrate=0.0;
		prevpassagepos[0]=0;
		prevpassagepos[1]=0;
		prevpassagepos[2]=0;
		closestapproach=numeric_limits<double>::max();
		mu=0.0;
		lx=0.0;
		ly=0.0;
		lz=0.0;
		ltot=0.0;
		E=0.0;
		phi=0;
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

	//Set the name for the header group
	string hdrname;

	//Set the name for the cosmology group
	string cosmohdrname;

	//Set the name for the units group
	string unithdrname;

	//Set the base name for the snapshot
	string grpbasename;

	// Store the names of the attributes in the header
	vector<H5std_string> hdrattrnames;

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
		hdrname = "/Header";
		cosmohdrname = "/Header/Cosmology";
		unithdrname = "/Header/Units";

		hdrattrnames.push_back("NSnaps");
		hdrattrnames.push_back("TEMPORALHALOIDVAL");

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

		datasetnames.push_back("ID");
		datasetnames.push_back("origID");
		datasetnames.push_back("Head");
		datasetnames.push_back("Tail");
		datasetnames.push_back("OrbitingHaloID");
		datasetnames.push_back("npart");
		datasetnames.push_back("hostFlag");
		datasetnames.push_back("Mass_200crit");
		datasetnames.push_back("R_200crit");
		datasetnames.push_back("Xc");
		datasetnames.push_back("Yc");
		datasetnames.push_back("Zc");
		datasetnames.push_back("VXc");
		datasetnames.push_back("VYc");
		datasetnames.push_back("VZc");
		datasetnames.push_back("Lx");
		datasetnames.push_back("Ly");
		datasetnames.push_back("Lz");
		datasetnames.push_back("Rmax");
		datasetnames.push_back("Vmax");
		datasetnames.push_back("cNFW");

		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasettypes.push_back(PredType::NATIVE_UINT8);
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
		datasettypes.push_back(PredType::STD_U64LE);
		datasetnames.push_back("orbithaloID");
		datasettypes.push_back(PredType::STD_I64LE);
		datasetnames.push_back("haloID");
		datasettypes.push_back(PredType::STD_I64LE);
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
		datasetnames.push_back("orbitecc__Wetzel2011");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbiteccratio");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitalEnergy");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rperi_Wetzel2011");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rapo_Wetzel2011");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("masslossrate");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LongAscNode");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Inclination");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("ArgPariap");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Phi");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("HostAlignment");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("MergerTimeScale");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("scalefactor");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("uniage");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("MergedFlag");
		datasettypes.push_back(PredType::NATIVE_UINT8);
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
		datasetnames.push_back("Npart");
		datasettypes.push_back(PredType::STD_I64LE);
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
		datasetnames.push_back("hostFlag");
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasetnames.push_back("Vrad");
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
		datasetnames.push_back("Npart_host");
		datasettypes.push_back(PredType::STD_I64LE);
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
		datasetnames.push_back("hostFlag_host");
		datasettypes.push_back(PredType::NATIVE_UINT8);
	};
};

#endif //ifndef ALLVARS_H






