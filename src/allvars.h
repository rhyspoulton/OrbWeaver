
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
#include <numeric>
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
#define NHDFFIELDSIN 26
#define NHDFFIELDSOUT 74

//Comoving or physical flags
#define PHYSICAL 0
#define COMOVING 1

//Parameters used to clean the orbits to remove apsis points that are "wobbles"
#define PHICLEANLIMIT (0.45 * M_PI)
#define ECCCLEANLIMIT 0.4
#define CLEANRATIOHOSTRADIUS 0.5 // The ratio of the host radius which cleaning is not done

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

	//The number of types of crossing entries
	int numtypeofcrossingentries;

	//The total number of types of different entries
	int totnumtypeofentries;

	//The fraction of rvir host that a datapoint is created
	float fracrvircross;

	//The file number currently working on
	int fileno;

	//The number of the hosts Rvir used to search for orbiting object,
	//this is also used as the upper limit for which crossing points are outputted.
	float numRvirSearch;

	//Flag to run the orbit cleaning
	bool iclean;

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
		numtypeofcrossingentries=0;
		totnumtypeofentries=0;
		fracrvircross=0;
		fileno=0;
		numRvirSearch=0.0;
		iclean=true;
		TEMPORALHALOIDVAL=0;
		iverbose=0;
	};
};




//Struct to load the halo into
struct HaloData{

	//The ID of the halo
	unsigned long long id;

	//The original ID of the halo
	unsigned long long origid;

	//The descendant of the halo
	long long descendant;

	//The progenitor of the halo
	long long progenitor;

	//The original root progenitor of the halo
	long long origrootprogenitor;

	//The original root descendant of the halo
	long long origrootdescendant;

	//If the halo is a subhalo
	bool fieldhalo;

	//If the orbit host merged in the simulation
	bool hostmerges;

	//The number ot substructure in the halo
	unsigned long long numsubstruct;

	//The ratio of mass in substructure for the halo
	double ratioofmassinsubsstruct;

	//The ID of the halo which it is currently orbiting
	long long orbitedhaloid;

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

	//Flag to say if this halo has been interpolated
	bool interpflag;

	//Flag to mark if the halo has been process or no
	bool doneflag;

	HaloData(){
		id=0;
		origid=0;
		descendant=0;
		progenitor=0;
		origrootprogenitor=0;
		origrootdescendant=0;
		fieldhalo=false;
		hostmerges=false;
		numsubstruct=0;
		ratioofmassinsubsstruct=0.0;
		orbitedhaloid=0;
		npart=0;
		mass=0.0;
		rvir=0.0;
		x=0.0;
		y=0.0;
		z=0.0;
		vx=0.0;
		vy=0.0;
		vz=0.0;
		lx=0.0;
		ly=0.0;
		lz=0.0;
		rmax=0.0;
		cnfw=0.0;
		interpflag=false;
		doneflag=false;
		interpflag=false;
	};
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
	};

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
	};
};



// Struct to store all the data for the halos
struct SnapData{
	unsigned long long numhalos;
	double scalefactor;
	double uniage; // Age of the universe at this snapshot
	double physboxsize; //The physical size of the box
	double Hz; // The Hubble parameter at this snapshot
	vector<HaloData> Halo;

	SnapData(){
		numhalos=0;
		scalefactor=0.0;
		uniage=0.0;
		Hz=0.0;
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

	//The orginal ID of the halo
	unsigned long long haloID;

	//The original Root Progenitor of the orbiting halo
	long long halorootprogenID;

	//The original Root Descendant of the orbiting halo
	long long halorootdescenID;

	//The original ID of its host halo
	long long orbitedhaloID;

	//The original Root Progenitor of the orbited halo from the halo catalog
	long long orbitedhaloorigrootprogenID;

	//The original Root Descendant of the orbited halo from the halo catalog
	long long orbitedhaloorigrootdescenID;

	//Type of datapoint this is 4 = apocenter, 3 = crossing 3x host's rvir, 2 = crossing 2x host's rvir, 1 = crossing 1x host's rvir, 0 = pericenter
	float entrytype;

	//Number of entries of the current entrytype
	int num_entrytype;

	//The number of orbits the halo has undergone
	float numorbits;

	//The total number of orbits that the satellite has undergone
	float totnumorbits;

	//The orbital period
	float orbitperiod;

	//Closest approach for this halo to the halo it is orbiting
	float closestapproach;

	//The scalefactor of the closest approach
	float closestapproachscalefactor;

	//The orbital eccentricity from ratio of passage distances
	float orbitecc_ratio;

	//The instananeous orbital energy
	float orbitalenergy_inst;

	//The average orbital Energy
	float orbitalenergy_ave;

	//The radius of a circular orbit with the same energy
	float rcirc;

	//The velocity of a circular orbit with the same energy
	float vcirc;

	//The angular momentum of a circular orbit with the same energy
	float jcirc;

	//The ratio of the orbital angular momentum to angular momentum of a circular orbit with the same energy
	float eta;

	//The peri-centric distance from orbital energy and angular momentum calculation
	float rperi_calc;

	//The apo-centric distance from orbital energy and angular momentum calculation
	float rapo_calc;

	//The orbital eccentricity from orbital energy and angular momentum calculation
	float orbitecc_calc;

	//Instantaneous mass loss rate
	float masslossrate_inst;

	//Average mass loss rate
	float masslossrate_ave;

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

	//The virial radius of the orbiting halo
	float rvir;

	//The minimum radius which half the mass is within
	float minrmax;

	//The minimum scale radius
	float minrscale;

	//The vmax of the orbiting halo
	float vmax;

	//The peak vmax up to this point
	float vmaxpeak;

	//The rmax of the orbiting halo
	float rmax;

	//The concentration of the orbiting halo
	float cnfw;

	//Boolean flag if this halo is a host halo (top of it spatial herachy)
	bool fieldhalo;

	//The number of substructure in this halo
	unsigned long long numsubstruct;

	//The ratio of mass in substructure for this halo
	float ratioofmassinsubsstruct;

	//The radial velocity of the halo with respect to its host
	float vrad;

	//The tangential velocity of the halo with respect to its host
	float vtan;

	//The radius of the orbit
	float r;

	//The relative x position to the host
	float xrel;

	//The relative y position to the host
	float yrel;

	//The relative z position to the host
	float zrel;

	//The total relative velocity
	float vrel;

	//The relative x velocity to the host
	float vxrel;

	//The relative y velocity to the host
	float vyrel;

	//The relative z velocity to the host
	float vzrel;

	//The total instantaneous orbital angular momentum
	float lrel_inst;

	//The instantaneous orbital angular momentum in x-direction
	float lxrel_inst;

	//The instantaneous orbital angular momentum in y-direction
	float lyrel_inst;

	//The instantaneous orbital angular momentum in z-direction
	float lzrel_inst;

	//The total average orbital angular momentum
	float lrel_ave;

	//The average orbital angular momentum in x-direction since last passage
	float lxrel_ave;

	//The average orbital angular momentum in y-direction since last passage
	float lyrel_ave;

	//The average orbital angular momentum in z-direction since last passage
	float lzrel_ave;

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
	bool fieldhalohost;

	//Boolean flag to mark if the host merged
	bool hostmerges;

	//The number of substructure in the host halo
	unsigned long long numsubstructhost;

	//The ratio of mass in substructure for the host ahlo
	float ratioofmassinsubsstructhost;

	OrbitData(){
		orbitID=0;
		orbithaloID=0;
		haloID=0;
		halorootprogenID=0;
		halorootdescenID=0;
		orbitedhaloID=0;
		orbitedhaloorigrootprogenID=0;
		orbitedhaloorigrootdescenID=0;
		entrytype=0.0;
		num_entrytype=0;
		numorbits=0;
		totnumorbits=0;
		orbitperiod=0.0;
		closestapproach=0.0;
		closestapproachscalefactor=0.0;
		orbitecc_calc=0.0;
		orbitecc_ratio=0.0;
		orbitalenergy_inst=0.0;
		orbitalenergy_ave=0.0;
		rcirc=0.0;
		vcirc=0.0;
		jcirc=0.0;
		eta=0.0;
		rperi_calc=0.0;
		rapo_calc=0.0;
		masslossrate_inst=0.0;
		masslossrate_ave=0.0;
		longascnode=0.0;
		inc=0.0;
		argpariap=0.0;
		phi=0.0;
		hostalignment=0.0;
		mergertimescale=0.0;
		scalefactor=0.0;
		uniage=0.0;
		mergedflag=false;
		x=0.0;
		y=0.0;
		z=0.0;
		vx=0.0;
		vy=0.0;
		vz=0.0;
		npart=0;
		mass=0.0;
		rvir=0.0;
		minrmax=0.0;
		minrscale=0.0;
		vmax=0.0;
		vmaxpeak=0.0;
		rmax=0.0;
		cnfw=0.0;
		fieldhalo=false;
		numsubstruct=0;
		ratioofmassinsubsstruct=0.0;
		vrad=0.0;
		vtan=0.0;
		r=0.0;
		xrel=0.0;
		yrel=0.0;
		zrel=0.0;
		vrel=0.0;
		vxrel=0.0;
		vyrel=0.0;
		vzrel=0.0;
		lrel_inst=0.0;
		lxrel_inst=0.0;
		lyrel_inst=0.0;
		lzrel_inst=0.0;
		lrel_ave=0.0;
		lxrel_ave=0.0;
		lyrel_ave=0.0;
		lzrel_ave=0.0;
		nparthost=0;
		rvirhost=0.0;
		masshost=0.0;
		vmaxhost=0.0;
		rmaxhost=0.0;
		cnfwhost=0.0;
		hostmerges=false;
		fieldhalohost=false;
		numsubstructhost=0;
		ratioofmassinsubsstructhost=0.0;
	};
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

	//Store the snapshot of the passage
	int passagesnap;

	//Store the snapshot which the previous passage happened
	int prevpassagesnap;

	//Value to keep track of the time of the previous apo/peri-centric pasage
	double prevpassagetime;

	//Store the index of this passage
	int passageindex;

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
	vector<double> refangles;

	//Store the radial vector for the previous position and the comoving radius
	vector<double> prevpassagepos;
	double prevpassagercomove;

	//The closest approach; that the halo had to its host
	double closestapproach;

	//The scalefactor that the closest approach happened at
	double closestapproachscalefactor;

	//The minimum rhalfmass and scale radius
	double minrmax;
	double minrscale;

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
	double E;
	double phi;

	OrbitProps(){
		orbitingflag = false;
		numorbits=0.0;
		prevcrossingentrytype=0.0;
		prevcrossingsnap = 0;
		passagesnap = 0;
		prevpassagesnap = 0;
		prevpassagetime = 0.0;
		prevpassageindex=0;
		crossrvirtime=0.0;
		passageindex=0;
		prevpassageindex=0;
		mergertime=0.0;
		masslossrate=0.0;
		refangles.resize(3);
		prevpassagepos.resize(3);
		prevpassagercomove=0.0;
		closestapproach=numeric_limits<double>::max();
		closestapproachscalefactor=0.0;
		minrmax=numeric_limits<double>::max();
		minrscale=numeric_limits<double>::max();
		mu=0.0;
		lx=0.0;
		ly=0.0;
		lz=0.0;
		hostlx=0.0;
		hostly=0.0;
		hostlz=0.0;
		E=0.0;
		phi=0;
	};
	~OrbitProps(){
		refangles.clear();
		refangles.shrink_to_fit();
		prevpassagepos.clear();
		prevpassagepos.shrink_to_fit();
	};
	OrbitProps &operator+=(OrbitProps &orbitprops){
		lx += orbitprops.lx;
		ly += orbitprops.ly;
		lz += orbitprops.lz;
		E += orbitprops.E;
		mu += orbitprops.mu;
		hostlx += orbitprops.hostlx;
		hostly += orbitprops.hostlx;
		hostlz += orbitprops.hostlx;
		masslossrate += orbitprops.masslossrate;
		phi += orbitprops.phi;
	}

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

		hdrattrnames.push_back("Fileno");
		hdrattrnames.push_back("NSnaps");
		hdrattrnames.push_back("TEMPORALHALOIDVAL");
		hdrattrnames.push_back("numRvirSearch");

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
		datasetnames.push_back("OrigID");
		datasetnames.push_back("Head");
		datasetnames.push_back("Tail");
		datasetnames.push_back("OrigRootProgenID");
		datasetnames.push_back("OrigRootDescenID");
		datasetnames.push_back("OrbitedHaloID");
		datasetnames.push_back("npart");
		datasetnames.push_back("FieldHalo");
		datasetnames.push_back("hostMerges");
		datasetnames.push_back("numSubStruct");
		datasetnames.push_back("RatioOfMassinSubsStruct");
		datasetnames.push_back("Mass");
		datasetnames.push_back("Radius");
		datasetnames.push_back("X");
		datasetnames.push_back("Y");
		datasetnames.push_back("Z");
		datasetnames.push_back("Vx");
		datasetnames.push_back("Vy");
		datasetnames.push_back("Vz");
		datasetnames.push_back("Lx");
		datasetnames.push_back("Ly");
		datasetnames.push_back("Lz");
		datasetnames.push_back("Rmax");
		datasetnames.push_back("Vmax");
		datasetnames.push_back("cNFW");

		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasettypes.push_back(PredType::NATIVE_ULLONG);
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
		datasettypes.push_back(PredType::NATIVE_DOUBLE);

	};
};


struct HDFOutputNames{

	// Store the names of the datasets
	vector<H5std_string> datasetnames;
	// Store their datatype
	vector<PredType> datasettypes;

	HDFOutputNames(){

		datasetnames.push_back("OrbitID");
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasetnames.push_back("HaloID_orbweaver");
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasetnames.push_back("HaloID_orig");
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasetnames.push_back("HaloRootProgen_orig");
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasetnames.push_back("HaloRootDescen_orig");
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasetnames.push_back("OrbitedHaloID_orig");
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasetnames.push_back("OrbitedHaloRootProgen_orig");
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasetnames.push_back("OrbitedHaloRootDescen_orig");
		datasettypes.push_back(PredType::NATIVE_LLONG);
		datasetnames.push_back("entrytype");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("num_entrytype");
		datasettypes.push_back(PredType::NATIVE_INT);
		datasetnames.push_back("numorbits");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("totnumorbits");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitalperiod");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("closestapproach");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("closestapproachscalefactor");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitecc_ratio");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitalenergy_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitalenergy_ave");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("R_circ");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("V_circ");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("L_circ");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Eta");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rperi_calc");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rapo_calc");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("orbitecc_calc");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("masslossrate_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("masslossrate_ave");
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
		datasetnames.push_back("npart");
		datasettypes.push_back(PredType::NATIVE_ULLONG);
		datasetnames.push_back("Mass");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Radius");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("min_Rmax");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("min_Rscale");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmax");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmaxpeak");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rmax");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("cNFW");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("FieldHalo");
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasetnames.push_back("numSubStruct");
		datasettypes.push_back(PredType::STD_I64LE);
		datasetnames.push_back("RatioOfMassinSubsStruct");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vrad");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vtan");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("R");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Xrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Yrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Zrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VXrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VYrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("VZrel");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Lrel_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LXrel_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LYrel_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LZrel_inst");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Lrel_ave");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LXrel_ave");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LYrel_ave");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("LZrel_ave");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("npart_host");
		datasettypes.push_back(PredType::STD_I64LE);
		datasetnames.push_back("Radius_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Mass_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Vmax_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("Rmax_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("cNFW_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
		datasetnames.push_back("FieldHalo_host");
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasetnames.push_back("hostMerges");
		datasettypes.push_back(PredType::NATIVE_UINT8);
		datasetnames.push_back("numSubStruct_host");
		datasettypes.push_back(PredType::STD_I64LE);
		datasetnames.push_back("RatioOfMassinSubsStruct_host");
		datasettypes.push_back(PredType::NATIVE_FLOAT);
	};
};

#endif //ifndef ALLVARS_H






