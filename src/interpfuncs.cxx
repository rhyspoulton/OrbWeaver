
#include "orbweaver.h"

double LogInterp(double prevdata, double nextdata, double f){
	return pow(nextdata,f) * pow(prevdata,1-f);
}
double LinInterp(double prevdata, double nextdata, double f){
	return prevdata + (nextdata - prevdata)*f;
}

//Use a cubic spline to interpolate the position and velocity of the halo
void InterpHaloPosVel(int nhalo, int ninterp, double *halouniages, double *interpuniages, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<SnapData> &snapdata, vector<HaloData> &interphalos){

	// Create a temporary dataset to store the known data points
	double tmpdata[nhalo];

	/*  Setup the interpolation routine  */

	//Intialize the acceleration and the spline routine
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nhalo);

	/*  x-pos  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
			tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].x;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].x = gsl_spline_eval(spline,interpuniages[i],acc);

	/*  y-pos  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].y;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].y = gsl_spline_eval(spline,interpuniages[i],acc);

	/*  z-pos  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].z;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].z = gsl_spline_eval(spline,interpuniages[i],acc);

	/*  x-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vx;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].vx = gsl_spline_eval(spline,interpuniages[i],acc);

	/*  y-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vy;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].vy = gsl_spline_eval(spline,interpuniages[i],acc);

	/*  z-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vz;

	//Intialize the data for the spline
	gsl_spline_init (spline, halouniages, tmpdata, nhalo);

	//Now do the interpolation by iterating through all the data
	for(int i = 0;i<ninterp;i++)
		interphalos[i].vz = gsl_spline_eval(spline,interpuniages[i],acc);

	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

}

HaloData InterpHaloProps(Options &opt, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<Int_t> &interpsnaps, vector<SnapData> &snapdata){

	//Set the number of halos and the number that need to interpolated
	Int_t nhalo = halosnaps.size(), ninterp = interpsnaps.size();

	//The halouniages has to be an array for the interpolation routine
	double halouniages[nhalo],f;
	double interpuniages[ninterp];
	int j=0;
	Int_t currentsnap,progensnap, progenindex, descsnap, descindex, orbitinghalosnap, orbitinghaloindex;

	//Lets extract the uniage for the interpolation routine for the snapshots that the halo is present and keep track
	for(int i = 0;i<nhalo;i++)
		halouniages[i] =  snapdata[halosnaps[i]].uniage;

	//The length of interpuniages is the amount of halos to be interpolated so can create the interpolated halos stuctures
	vector<HaloData> interphalos;
	interphalos.resize(ninterp);

	//Convert the snapshots which need to be interpolated into ages of the universe
	for(int i=0;i<ninterp;i++)
		interpuniages[i] = snapdata[interpsnaps[i]].uniage;

	InterpHaloPosVel(nhalo,ninterp,halouniages,interpuniages,halosnaps,haloindexes,snapdata,interphalos);

	for(int i = 0;i<nhalo-1;i++){

		//Note the progenitor and descendants for when interpolation needs to be done
		progensnap = halosnaps[i];
		progenindex = haloindexes[i];
		descsnap = halosnaps[i+1];
		descindex = haloindexes[i+1];

		//Lets move to descendant to see if it exist in the next snapshot
		currentsnap = progensnap+1;

		//Keep track of the halo that the progen is orbiting
		orbitinghalosnap = (Int_t)(snapdata[progensnap].Halo[progenindex].orbitinghaloid/opt.TEMPORALHALOIDVAL);
		orbitinghaloindex = (Int_t)(snapdata[progensnap].Halo[progenindex].orbitinghaloid%opt.TEMPORALHALOIDVAL-1);

		// Iterate until currentsnap==descsnap and extract
		while(currentsnap!=descsnap){

			//Find what will be the ID of this halo
			interphalos[j].id = currentsnap*opt.TEMPORALHALOIDVAL + snapdata[currentsnap].numhalos+1;

			//If at the first interpolated halo have the progenitor point to this halo
			//and this halo point back to the halo as a progenitor. Otherwise have it point
			//to the previous interpolated halo and have the previous halo point to this halo
			//as its descendant
			if(currentsnap==progensnap+1){
				snapdata[progensnap].Halo[progenindex].descendant = interphalos[j].id;
				interphalos[j].progenitor = snapdata[progensnap].Halo[progenindex].id;
			}
			else{
				interphalos[j-1].descendant = interphalos[j].id;
				interphalos[j].progenitor = interphalos[j-1].id;
			}

			/* Now do a logarithmic interpolation for all other properties */

			//Find the difference in the fraction of uniage to find where the intepolation should be done
			f = (snapdata[currentsnap].uniage - snapdata[halosnaps[i]].uniage)/(snapdata[halosnaps[i+1]].uniage - snapdata[halosnaps[i]].uniage);

			//Number of particles in the halo
			interphalos[j].npart = (unsigned long long)LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].npart,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].npart,f);

			//Virial mass of the halo
			interphalos[j].mass = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].mass,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].mass,f);

			//Virial radius of the halo
			interphalos[j].rvir = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].rvir,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].rvir,f);

			//Rmax of the halo
			interphalos[j].rmax = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].rmax,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].rmax,f);

			//Vmax of the halo
			interphalos[j].vmax = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].vmax,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].vmax,f);

			//cNFW of the halo
			interphalos[j].cnfw = LinInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].cnfw,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].cnfw,f);

			//Check if we will be at the snapshot of the orginal descendant, if so
			//lets have the interpolated halo point to it and the descendant point
			//to the interpolated halo
			if(currentsnap+1==descsnap){
				snapdata[descsnap].Halo[descindex].progenitor = interphalos[j].id;
				interphalos[j].descendant = snapdata[descsnap].Halo[descindex].id;
			}

			//Now need to set the halo it is orbiting to be the descendant of the host halo
			//that the previous halo was orbiting and then move on to descendant halo
			interphalos[j].orbitinghaloid = snapdata[orbitinghalosnap].Halo[orbitinghaloindex].descendant;
			orbitinghalosnap = (Int_t)(interphalos[j].orbitinghaloid/opt.TEMPORALHALOIDVAL);
			orbitinghaloindex = (Int_t)(interphalos[j].orbitinghaloid%opt.TEMPORALHALOIDVAL-1);

			if(currentsnap!=orbitinghalosnap)
				cout<<"Warning: this halo is set to orbit a host halo thats at a different snapshot, have the host halos be interpolated?"<<endl;

			//Iterate the snapshot
			currentsnap++;
			j++;
		}
	}
	//Now lets add this interpolated halo into the snapsdata at each of the interpolation
	//snapshots and also add to the number of halos at each snapshot
	for(int i=0;i<ninterp;i++){
		snapdata[interpsnaps[i]].Halo.push_back(interphalos[i]);

		//insert these halos into the halosnaps and indexes
		halosnaps.insert(halosnaps.begin()+interpsnaps[i]-halosnaps[0],interpsnaps[i]);
		haloindexes.insert(haloindexes.begin()+interpsnaps[i]-halosnaps[0],snapdata[interpsnaps[i]].numhalos);
		snapdata[interpsnaps[i]].numhalos++;
	}
}

void InterpPassageHaloProps(double interpuniage, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata){

	double f = (interpuniage - currentuniage)/(prevuniage - currentuniage);

	//Interpolate the orbiting halo's properties
	tmporbitdata.npart = (unsigned long long)LogInterp(prevorbitinghalo.npart,orbitinghalo.npart,f);
	tmporbitdata.mass = LogInterp(prevorbitinghalo.mass,orbitinghalo.mass,f);
	tmporbitdata.vmax = LogInterp(prevorbitinghalo.vmax,orbitinghalo.vmax,f);
	tmporbitdata.rmax = LogInterp(prevorbitinghalo.rmax,orbitinghalo.rmax,f);
	tmporbitdata.cnfw = LinInterp(prevorbitinghalo.cnfw,orbitinghalo.cnfw,f);

	//Interpolate the host halo properties
	tmporbitdata.nparthost = (unsigned long long)LogInterp(prevhosthalo.npart,hosthalo.npart,f);
	tmporbitdata.masshost = LogInterp(prevhosthalo.mass,hosthalo.mass,f);
	tmporbitdata.rvirhost = LogInterp(prevhosthalo.rvir,hosthalo.rvir,f);
	tmporbitdata.vmaxhost = LogInterp(prevhosthalo.vmax,hosthalo.vmax,f);
	tmporbitdata.rmaxhost = LogInterp(prevhosthalo.rmax,hosthalo.rmax,f);
	tmporbitdata.cnfwhost = LinInterp(prevhosthalo.cnfw,hosthalo.cnfw,f);
}

void InterpPassagePoints(vector<Int_t> halosnaps, vector<Int_t> haloindexes,vector<Int_t> hostindexes, vector<SnapData> &snapdata, vector<OrbitData> &branchorbitdata, OrbitProps &orbitprops){

	int nhalo = halosnaps.size();
	int ninterp = branchorbitdata.size();
	double halouniages[nhalo];
	double interpuniages[ninterp];
	double r;
	bool merged=false;

	//Lets extract the uniage for the interpolation routine for the snapshots that the halo is present and keep track
	for(int i = 0;i<nhalo;i++)
		halouniages[i] =  snapdata[halosnaps[i]].uniage;

	//Now extract all the universe age for all the points in branchorbitdata
	for(int i = 0; i<ninterp;i++)
		interpuniages[i] = branchorbitdata[i].uniage;


	vector<HaloData> interphalos;
	interphalos.resize(ninterp);
	vector<HaloData> interphosthalos;
	interphosthalos.resize(ninterp);

	InterpHaloPosVel(nhalo, ninterp, halouniages, interpuniages, halosnaps, haloindexes,snapdata,interphalos);
	InterpHaloPosVel(nhalo, ninterp, halouniages, interpuniages, halosnaps, hostindexes,snapdata,interphosthalos);

	for(int i = 0; i<ninterp;i++){


		//Do a periodicity correction
		if((interphalos[i].x - interphosthalos[i].x)>0.5*Cosmo.boxsize) interphalos[i].x-=Cosmo.boxsize;
		if((interphalos[i].y - interphosthalos[i].y)>0.5*Cosmo.boxsize) interphalos[i].y-=Cosmo.boxsize;
		if((interphalos[i].z - interphosthalos[i].z)>0.5*Cosmo.boxsize) interphalos[i].z-=Cosmo.boxsize;
		if((interphalos[i].x - interphosthalos[i].x)<-0.5*Cosmo.boxsize) interphalos[i].x+=Cosmo.boxsize;
		if((interphalos[i].y - interphosthalos[i].y)<-0.5*Cosmo.boxsize) interphalos[i].y+=Cosmo.boxsize;
		if((interphalos[i].z - interphosthalos[i].z)<-0.5*Cosmo.boxsize) interphalos[i].z+=Cosmo.boxsize;

		branchorbitdata[i].x = interphalos[i].x;
		branchorbitdata[i].y = interphalos[i].y;
		branchorbitdata[i].z = interphalos[i].z;
		branchorbitdata[i].vx = interphalos[i].vx;
		branchorbitdata[i].vy = interphalos[i].vy;
		branchorbitdata[i].vz = interphalos[i].vz;
		branchorbitdata[i].xrel = interphosthalos[i].x - interphalos[i].x;
		branchorbitdata[i].yrel = interphosthalos[i].y - interphalos[i].y;
		branchorbitdata[i].zrel = interphosthalos[i].z - interphalos[i].z;
		branchorbitdata[i].vxrel = interphosthalos[i].vx - interphalos[i].vx;
		branchorbitdata[i].vyrel = interphosthalos[i].vy - interphalos[i].vy;
		branchorbitdata[i].vzrel = interphosthalos[i].vz - interphalos[i].vz;

		//Lets see if after the interpolation of the halo, it has merged with its host 
		r = sqrt(branchorbitdata[i].xrel*branchorbitdata[i].xrel + branchorbitdata[i].yrel*branchorbitdata[i].yrel + branchorbitdata[i].zrel*branchorbitdata[i].zrel);
		if((branchorbitdata[i].entrytype!=0.0) & (branchorbitdata[i].entrytype!=-99)) cout<<branchorbitdata[i].entrytype<<" "<<r/branchorbitdata[i].rvirhost<<endl;
		if((branchorbitdata[i].entrytype!=0.0) & (branchorbitdata[i].entrytype!=-99) & (abs((r/branchorbitdata[i].rvirhost) - abs(branchorbitdata[i].entrytype))>5.0)){
			cout<<branchorbitdata[i].haloID<<" "<<branchorbitdata[i].entrytype<<" "<<r/branchorbitdata[i].rvirhost<<" "<<abs((r/branchorbitdata[i].rvirhost) - abs(branchorbitdata[i].entrytype))<<" "<<r<<" "<<branchorbitdata[i].npart<<" "<<branchorbitdata[i].nparthost<<endl;
		}

		if((r<0.1 * interphosthalos[i].rvir) & (interpuniages[i]<orbitprops.mergertime)){
			orbitprops.mergertime=interpuniages[i];
			merged=true;
		}

		if(merged){
			branchorbitdata[i].mergedflag=true;
		}
	}
}