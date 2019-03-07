

#include "orbweaver.h"


///Checks if file exits by attempting to get the file attributes
///If success file obviously exists.
///If failure may mean that we don't have permission to access the folder which contains this file or doesn't exist.
///If we need to do that level of checking, lookup return values of stat which will give you more details on why stat failed.
bool FileExists(const char *fname) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  intStat = stat(fname,&stFileInfo);
  if(intStat == 0) return  true;
  else return false;
}

vector<HaloData> ReadSnapshotData(int snap, Group snapgroup, Options &opt, vector<SnapData> &snapdata, HDFCatalogNames hdfnames){
	
	int ichunk, chunksize=8192;
	Attribute snapattr;
	DataSpace attrdataspace;
	IntType inttype;
	FloatType floattype;
	DataSet *halosdataset;
	DataSpace *halosdataspace;
	DataSpace idataspace;
	int rank;
	hsize_t dim[1];
	vector<HaloData> Halo;

	hsize_t filespacecount[2],filespaceoffset[2];

	//buffers to load data
	int *intbuff=new int[chunksize];
	unsigned long long *ulongbuff=new unsigned long long[chunksize];
	long long *longbuff=new long long[chunksize];
	unsigned int *uintbuff=new unsigned int[chunksize];
	float *floatbuff=new float[chunksize*3];
	double *doublebuff=new double[chunksize*3];
	bool *boolbuff = new bool[chunksize];

	//Keep track of which field that is being read
	int ifield, count,offset=0;
	long long numread;

	//Setup the hdf5 objects
	halosdataset =  new DataSet[NHDFFIELDSIN];
	halosdataspace = new DataSpace[NHDFFIELDSIN];


	if(opt.iverbose) cout<<"Loading snap "<<snap<<endl;
	// Get the number of halo attribute for the snapshot
	snapattr = snapgroup.openAttribute(hdfnames.snapattrnames[0]);
	attrdataspace = snapattr.getSpace();
	inttype = snapattr.getIntType();
	if (inttype.getSize()==sizeof(int)) {
		snapattr.read(PredType::NATIVE_INT,&intbuff[0]);
		snapdata[snap].numhalos=intbuff[0];
	}
	if (inttype.getSize()==sizeof(long long)) {
		snapattr.read(PredType::NATIVE_LONG,&longbuff[0]);
		snapdata[snap].numhalos=longbuff[0];
	}
	// Get the scalefactor attribute for the snapshot
	snapattr = snapgroup.openAttribute(hdfnames.snapattrnames[1]);
	attrdataspace = snapattr.getSpace();
	floattype = snapattr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		snapattr.read(PredType::NATIVE_FLOAT,&floatbuff[0]);
		snapdata[snap].scalefactor=floatbuff[0];
	}
	if (floattype.getSize()==sizeof(double)) {
		snapattr.read(PredType::NATIVE_DOUBLE,&doublebuff[0]);
		snapdata[snap].scalefactor=doublebuff[0];
	}

	numread = snapdata[snap].numhalos;

	//Setup the local halodata
	Halo.resize(numread);


	//First lets open up the datasets
	for(int i=0;i<NHDFFIELDSIN;i++){
		halosdataset[i] = snapgroup.openDataSet(hdfnames.datasetnames[i]);
		halosdataspace[i] = halosdataset[i].getSpace();
	}


	if(snapdata[snap].numhalos<chunksize) ichunk = numread;
	else ichunk = chunksize;
	

	for(int n=0;n<snapdata[snap].numhalos;n+=ichunk){


		if((numread-n<chunksize) && (numread-n>0)) ichunk = numread-n;

		// ID
		ifield=0;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(ulongbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);


		for(int nn=0;nn<ichunk;nn++) Halo[count++].id = ulongbuff[nn];

		// original ID
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(ulongbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].origid = ulongbuff[nn];

		// Descendant
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(longbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].descendant = longbuff[nn];

		// Progenitor
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(longbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].progenitor = longbuff[nn];


		// orbitedhaloid
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(longbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].orbitedhaloid = longbuff[nn];


		// orbitedhaloorigrootprogen
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(longbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].orbitedhaloorigrootprogen = longbuff[nn];


		// Npart
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(ulongbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].npart = ulongbuff[nn];

		// fieldhalo
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(boolbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].fieldhalo = boolbuff[nn];

		// hostMerges
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(boolbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].hostmerges = boolbuff[nn];

		// numSubStruct
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(ulongbuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].numsubstruct = ulongbuff[nn];

		// RatioOfMassinSubsStruct
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].ratioofmassinsubsstruct = doublebuff[nn];

		// M_200crit
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].mass = doublebuff[nn] * Units.mass;

		// R_200crit
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].rvir = doublebuff[nn] * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].rvir = doublebuff[nn] * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Xc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].x = doublebuff[nn] * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].x = doublebuff[nn] * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Yc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].y = doublebuff[nn] * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].y = doublebuff[nn] * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Zc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].z = doublebuff[nn] * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].z = doublebuff[nn] * Cosmo.h/snapdata[snap].scalefactor * Units.length;
	

		// VXc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].vx = doublebuff[nn] * Units.velocity;

		// VYc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].vy = doublebuff[nn] * Units.velocity;

		// VZc
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].vz = doublebuff[nn] * Units.velocity;

		// Lx
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].lx = doublebuff[nn] * Units.mass * Units.velocity * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].lx = doublebuff[nn] * Units.mass * Units.velocity * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Ly
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].ly = doublebuff[nn] * Units.mass * Units.velocity * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].ly = doublebuff[nn] * Units.mass * Units.velocity * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Lz
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].lz = doublebuff[nn] * Units.mass * Units.velocity * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].lz = doublebuff[nn] * Units.mass * Units.velocity * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Rmax
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		if(Units.distFlag==COMOVING) for(int nn=0;nn<ichunk;nn++) Halo[count++].rmax = doublebuff[nn] * Units.length;
		else for(int nn=0;nn<ichunk;nn++) Halo[count++].rmax = doublebuff[nn] * Cosmo.h/snapdata[snap].scalefactor * Units.length;

		// Vmax
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].vmax = doublebuff[nn] * Units.velocity;

		// cNFW
		ifield++;
		count=offset;
		rank = 1;
		dim[0] = ichunk;
		idataspace = DataSpace(rank,dim);
		filespacecount[0]=ichunk;filespacecount[1]=1;
		filespaceoffset[0]=n;filespaceoffset[1]=0;
		halosdataspace[ifield].selectHyperslab(H5S_SELECT_SET,filespacecount,filespaceoffset);
		halosdataset[ifield].read(doublebuff,hdfnames.datasettypes[ifield],idataspace,halosdataspace[ifield]);

		for(int nn=0;nn<ichunk;nn++) Halo[count++].cnfw = doublebuff[nn];

		offset+=ichunk;	
		}


	delete[] doublebuff;
	delete[] longbuff;


				
	
	return Halo;
}


void ReadHeader(Options &opt, H5File *Fhdf,HDFCatalogNames hdfnames){

	float floatbuff;
	double doublebuff;
	bool boolbuff;
	unsigned long long ullongbuf;
	int intbuff;
	Group hdrgroup;
	Group cosmogroup;
	Group unitgroup;
	Attribute attr;
	DataSpace attrdataspace;
	FloatType floattype;

	//Open up the header group and extract the information
	hdrgroup = Fhdf->openGroup(hdfnames.hdrname);

	//The file number for this file
	attr = hdrgroup.openAttribute(hdfnames.hdrattrnames[0]);
	attrdataspace = attr.getSpace();
	attr.read(PredType::NATIVE_INT,&intbuff);
	opt.fileno=intbuff;

	//The number of snapshots in the simulation
	attr = hdrgroup.openAttribute(hdfnames.hdrattrnames[1]);
	attrdataspace = attr.getSpace();
	attr.read(PredType::NATIVE_INT,&intbuff);
	opt.numsnaps=intbuff;

	//The temporal haloID value to make halo IDs temporally unique
	attr = hdrgroup.openAttribute(hdfnames.hdrattrnames[2]);
	attrdataspace = attr.getSpace();
	attr.read(PredType::STD_U64LE,&ullongbuf);
	opt.TEMPORALHALOIDVAL=ullongbuf;

	//Open up the Units header group and extract the information
	unitgroup = Fhdf->openGroup(hdfnames.unithdrname);

	// Physical or comoving flag
	attr = unitgroup.openAttribute(hdfnames.unitattrnames[0]);
	attrdataspace = attr.getSpace();
	attr.read(PredType::NATIVE_HBOOL,&boolbuff);
	Units.distFlag=boolbuff;

	// Length in Mpc
	attr = unitgroup.openAttribute(hdfnames.unitattrnames[1]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Units.length=floatbuff/1000.0;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Units.length=doublebuff/1000.0;
	}
	// Mass in Msol
	attr = unitgroup.openAttribute(hdfnames.unitattrnames[2]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Units.mass=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Units.mass=doublebuff;
	}

	// Velocity in km/s
	attr = unitgroup.openAttribute(hdfnames.unitattrnames[3]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Units.velocity=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Units.velocity=doublebuff;
	}

	//Open up the cosmology header group and extract the information
	cosmogroup = Fhdf->openGroup(hdfnames.cosmohdrname);

	// boxsize
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[0]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.boxsize=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.boxsize=doublebuff;
	}

	// h
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[1]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.h=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.h=doublebuff;
	}

	// G
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[2]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.G=floatbuff /Units.mass;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.G=doublebuff /Units.mass;
	}

	// Omega Lambda
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[3]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.omegaL=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.omegaL=doublebuff;
	}

	// Omega matter
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[4]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.omegaM=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.omegaM=doublebuff;
	}

	// Omega radiation
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[5]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.omegaR=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.omegaR=doublebuff;
	}

	// Omega k
	attr = cosmogroup.openAttribute(hdfnames.cosmoattrnames[6]);
	attrdataspace = attr.getSpace();
	floattype = attr.getFloatType();
	if (floattype.getSize()==sizeof(float)) {
		attr.read(PredType::NATIVE_FLOAT,&floatbuff);
		Cosmo.omegaK=floatbuff;
	}
	if (floattype.getSize()==sizeof(double)) {
		attr.read(PredType::NATIVE_DOUBLE,&doublebuff);
		Cosmo.omegaK=doublebuff;
	}

	return;
}

void ReadData(Options &opt, vector<SnapData> &snapdata){

	H5File *Fhdf;
	Fhdf=new H5File;
	HDFCatalogNames hdfnames;
	char buff[100];
	const char* grpbasename;
	string grpname;
	Group snapgroup;
	int i;



	// Check if the file exists
	if(!FileExists(opt.fname)){

		printf("IOError: Cannot open file `%s'\n",opt.fname);
		// TerminateRun(2);
		EXIT_CODE = 2;
		throw exception();
	}

	try
	{
		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		Exception::dontPrint();

		//Open up the file
		Fhdf->openFile(opt.fname,H5F_ACC_RDONLY);

		//Read the hdf header info
		ReadHeader(opt,Fhdf,hdfnames);

		//Now the header information has been read, the number of snapshots is now known so snapdata can be resized
		snapdata.resize(opt.numsnaps);

		//Get the base name of the groups.
		grpbasename = hdfnames.grpbasename.c_str();

		for(int snap=0;snap<opt.numsnaps;snap++){

			// Open up the snapgroup
			sprintf(buff,grpbasename,snap);
			grpname = string(buff);
			snapgroup = Fhdf->openGroup(grpname);

			//Read the required data
			snapdata[snap].Halo = ReadSnapshotData(snap,snapgroup,opt,snapdata,hdfnames);

		}

		Fhdf->close();
	} 

	catch(GroupIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the H5File operations
	catch(FileIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the DataSet operations
	catch(DataSetIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the DataSpace operations
	catch(DataSpaceIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	return;

}

void WriteOrbitData(Options &opt, vector<OrbitData> &orbitdata){

	H5File file;
	H5std_string datasetname;
	DataSpace dataspace;
	DataSet dataset;
	Attribute attr;
	DataSpace attrspace;
	DSetCreatPropList hdfdatasetproplist;
	hsize_t dims[1], chunk_dims[1];
	int rank;
	HDFCatalogNames hdfnames;
	int itemp=0;
	HDFOutputNames hdfdatasetnames;
	char outfilename[1000];
	unsigned long long numentries = orbitdata.size();

	//Check the number of entries is >0
	if(numentries==0){
		cerr<<"Error: The number of entries to output is 0"<<endl;
		EXIT_CODE = 2;
		throw exception();
	}

	//Create the buffers to load the data into
	unsigned long long ullongbuf[numentries];
	long long llongbuf[numentries];
	int intbuff[numentries];
	float floatbuff[numentries];
	bool boolbuff[numentries];

	//Lets setup the info for the datasets
	dims[0] = numentries;
	rank=1;
	chunk_dims[0] = min((unsigned long long)HDFOUTCHUNKSIZE,(unsigned long long)numentries);

	int idataset = 0;

	try{
		// Turn off the auto-printing when failure occurs so that we can
		// handle the errors appropriately
		Exception::dontPrint();

		//Get the filename
		strcpy(outfilename,opt.outputbasename);
		strcat(outfilename,".orbweaver.orbitdata.hdf");
		cout<<"Outputting orbitdata to file "<<outfilename<<endl;

		//Open up the file
		file = H5File(outfilename,H5F_ACC_TRUNC);

		//Create the attributes for the file

		//The number of entries in the file
		attrspace=DataSpace(H5S_SCALAR);
		attr=file.createAttribute("Fileno", PredType::NATIVE_INT, attrspace);
		attr.write(PredType::NATIVE_INT,&opt.fileno);

		//The number of entries in the file
		attrspace=DataSpace(H5S_SCALAR);
		attr=file.createAttribute("Number_of_entries", PredType::STD_U64LE, attrspace);
		attr.write(PredType::STD_U64LE,&numentries);

		//The inputfile which was processed to produce this catalog
		attrspace=DataSpace(H5S_SCALAR);
		// Create new string datatype for attribute
		StrType strdatatype(PredType::C_S1, 1000);
		// Set up write buffer for attribute
		const H5std_string strwritebuf (opt.fname);
		attr = file.createAttribute("Input_preprocessed_catalog", strdatatype, attrspace);
		attr.write(strdatatype, strwritebuf);

		//The max orbitID in this catalog
		attrspace=DataSpace(H5S_SCALAR);
		attr=file.createAttribute("Max_orbitID", PredType::STD_U64LE, attrspace);
		attr.write(PredType::STD_U64LE,&orbitdata[numentries-1].orbitID);


		/* orbitID */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].orbitID;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbithaloID */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].orbithaloID;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* haloID */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].haloID;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbitedhaloID */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) llongbuf[j] = orbitdata[j].orbitedhaloID;
		dataset.write(llongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbitedhaloorigrootprogen */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) llongbuf[j] = orbitdata[j].orbitedhaloorigrootprogen;
		dataset.write(llongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* entrytype */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].entrytype;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* numentrytype */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) intbuff[j] = orbitdata[j].num_entrytype;
		dataset.write(intbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* numorbits */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].numorbits;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbitperiod */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].orbitperiod;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* closestapproach */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].closestapproach * orbitdata[j].scalefactor / Cosmo.h;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbitecc */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].orbitecc;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbiteccratio */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].orbiteccratio;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* orbitalenergy */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].orbitalenergy;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* rperi */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rperi * orbitdata[j].scalefactor / Cosmo.h;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* rapo */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rapo * orbitdata[j].scalefactor / Cosmo.h;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* masslossrate_inst */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].masslossrate_inst;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* masslossrate_ave */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].masslossrate_ave;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* longascnode */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].longascnode;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* inc */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].inc;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* argpariap */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].argpariap;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* Phi */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].phi;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* hostalignment */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].hostalignment;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* mergertimescale */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].mergertimescale;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* scalefactor */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].scalefactor;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* uniage */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].uniage;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* mergedflag */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) boolbuff[j] = orbitdata[j].mergedflag;
		dataset.write(boolbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* x */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].x * orbitdata[j].scalefactor / Cosmo.h;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* y */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].y * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* z */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].z * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vx */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vx;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vy */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vy;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vz */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vz;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* npart */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].npart;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* mass */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].mass;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* rvir */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rvir;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vmax */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vmax;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vmaxpeak */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vmaxpeak;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* rmax */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rmax * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* cnfw */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].cnfw;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* fieldhalo */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) boolbuff[j] = orbitdata[j].fieldhalo;
		dataset.write(boolbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;


		/* numSubStruct */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].numsubstruct;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;


		/* RatioOfMassinSubStruct */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].ratioofmassinsubsstruct;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vrad */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vrad;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vtan */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vtan;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* xrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].xrel * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* yrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].yrel * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* zrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].zrel * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vxrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vxrel;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vyrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vyrel;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vzrel */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vzrel;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lxrel_inst */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lxrel_inst * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lyrel_inst */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lyrel_inst * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lzrel_inst */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lzrel_inst * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lxrel_ave */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lxrel_ave * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lyrel_ave */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lyrel_ave * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* lzrel_ave */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].lzrel_ave * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* nparthost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].nparthost;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;


		/* rvirhost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rvirhost * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* masshost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].masshost;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* vmaxhost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].vmaxhost;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* rmaxhost */


		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].rmaxhost * orbitdata[j].scalefactor / Cosmo.h;;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* cnfwhost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].cnfwhost;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;


		/* fieldhalohost */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) boolbuff[j] = orbitdata[j].fieldhalohost;
		dataset.write(boolbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* hostMerges */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) boolbuff[j] = orbitdata[j].hostmerges;
		dataset.write(boolbuff,hdfdatasetnames.datasettypes[idataset]);
		idataset++;

		/* numSubStruct_host */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) ullongbuf[j] = orbitdata[j].numsubstructhost;
		dataset.write(ullongbuf,hdfdatasetnames.datasettypes[idataset]);
		idataset++;


		/* RatioOfMassinSubStruct_host */

		//Create the dataset
		dataspace = DataSpace(rank,dims);

		if(chunk_dims[0]>0){

			hdfdatasetproplist=DSetCreatPropList();
			// Modify dataset creation property to enable chunking
			hdfdatasetproplist.setChunk(rank, chunk_dims);
			// Set ZLIB (DEFLATE) Compression using level 6.
			hdfdatasetproplist.setDeflate(6);

			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace,hdfdatasetproplist);

		}
		else{
			dataset = file.createDataSet(hdfdatasetnames.datasetnames[idataset],hdfdatasetnames.datasettypes[idataset],dataspace);
		}

		//Write out the dataset
		for(unsigned int j=0; j<numentries;j++) floatbuff[j] = orbitdata[j].ratioofmassinsubsstructhost;
		dataset.write(floatbuff,hdfdatasetnames.datasettypes[idataset]);

		file.close();
	}
	catch(GroupIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the H5File operations
	catch(FileIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the DataSet operations
	catch(DataSetIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}

	// catch failure caused by the DataSpace operations
	catch(DataSpaceIException error)
	{
		error.printErrorStack();
		EXIT_CODE = 2;
		throw exception();
	}
}
