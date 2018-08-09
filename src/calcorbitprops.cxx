

#include "orbweaver.h"

OrbitData CalcOrbitProps(HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &prevorbitdata, double scalefactor){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*Cosmo.boxsize){
		cout<<"he position was "<<orbitinghalo.x<<endl;
		orbitinghalo.x-=Cosmo.boxsize;
		cout<<"Corrected the position to "<<orbitinghalo.x<<endl;
	}

	if((orbitinghalo.y - hosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

	//Create a temporary orbitdata structure for the current halo
	OrbitData orbitdata = prevorbitdata;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r,vr;

	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = hosthalo.x - orbitinghalo.x;
	ry = hosthalo.y - orbitinghalo.y;
	rz = hosthalo.z - orbitinghalo.z;
	vrx = hosthalo.vx - orbitinghalo.vx;
	vry = hosthalo.vy - orbitinghalo.vy;
	vrz = hosthalo.vz - orbitinghalo.vz;
	r = sqrt(rx * rx + ry * ry + rz * rz);
	vr = (rx * vrx + ry * vry * rz * vrz) / r;

	//Lets check if we are at the base of this branch i.e. the progenitor ID is the same as the halo's
	if((orbitinghalo.progenitor==orbitinghalo.id) & (orbitinghalo.id!=0)){
		orbitdata.closestapproach = r;
		return orbitdata;
	}

	double prevrx,prevry,prevrz,prevvrx,prevvry,prevvrz,prevr,prevvr;

	//Lets find the same for the previous halo
	prevrx = prevhosthalo.x - prevorbitinghalo.x;
	prevry = prevhosthalo.y - prevorbitinghalo.y;
	prevrz = prevhosthalo.z - prevorbitinghalo.z;
	prevvrx = prevhosthalo.vx - prevorbitinghalo.vx;
	prevvry = prevhosthalo.vy - prevorbitinghalo.vy;
	prevvrz = prevhosthalo.vz - prevorbitinghalo.vz;
	prevr = sqrt(prevrx * prevrx + prevry * prevry + prevrz * prevrz);
	prevvr = (prevrx * prevvrx + prevry * prevvry * prevrz * prevvrz) / r;


	//Lets find if this halo has the closest approach so far
	if(r<prevorbitdata.closestapproach){
		// cout<<"Updating this halos closestapproach "<<r<<endl;
		prevorbitdata.closestapproach=0;
		orbitdata.closestapproach = r;
	}


	// Now lets see if the mulitplication of the two orbiting vector gives out a -ve number (has undergone 1/2 an orbit)
	if(vr*prevvr<0){
		orbitdata.numorbits = prevorbitdata.numorbits + 0.5;
	}


	// Find the mass loss rate in units of myr?

	return orbitdata;

}


//Use a cubic spline to interpolate the position and velocity of the halo
void InterpHaloPosVel(int nhalo, int ninterp, double *halouniages, vector<double> &interpuniages, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, SnapData *&snapdata, vector<HaloData> &interphalos){

	// Create a temporary dataset to store the known data points
	double tmpdata[nhalo];

	/*  Setup the interpolation routine  */

	//Intialize the acceleration and the spline routine
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, nhalo);

	//The length of interpuniages is the amount of halos to be interpolated
	ninterp = interpuniages.size();

	/*  x-pos  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		tmpdata[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].x;
	}

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


double LogInterp(double prevdata, double nextdata, double f){
	return pow(nextdata,f) * pow(prevdata,1-f);
}

HaloData InterpHaloProps(Options &opt, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<Int_t> &interpsnaps, SnapData *&snapdata){


	//Set the number of halos and the number that need to interpolated
	Int_t nhalo = halosnaps.size(), ninterp = interpsnaps.size();

	//The halouniages has to be an array for the interpolation routine
	double halouniages[nhalo],f;
	vector <double> interpuniages(ninterp);
	int j=0;
	Int_t currentsnap,progensnap, progenindex, descsnap, descindex, orbitinghalosnap, orbitinghaloindex;

	//Lets extract the uniage for the interpolation routine for the snapshots that the halo is present and keep track
	for(int i = 0;i<nhalo;i++)
		halouniages[i] =  snapdata[halosnaps[i]].uniage;

	//The length of interpuniages is the amount of halos to be interpolated so can create the interpolated halos stuctures
	vector<HaloData> interphalos;
	interphalos.reserve(ninterp);

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
			interphalos[j].id = currentsnap*opt.TEMPORALHALOIDVAL + snapdata[currentsnap].numhalos;

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

			//Virial mass of the halo
			interphalos[j].mvir = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].mvir,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].mvir,f);

			//Virial radius of the halo
			interphalos[j].rvir = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].rvir,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].rvir,f);

			//Rmax of the halo
			interphalos[j].rmax = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].rmax,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].rmax,f);

			//Vmax of the halo
			interphalos[j].vmax = LogInterp(snapdata[halosnaps[i]].Halo[haloindexes[i]].vmax,snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].vmax,f);


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
			snapdata[interpsnaps[i]].numhalos++;
	}

}


void ProcessHalo(Int_t snap, Int_t i, Options &opt, SnapData *&snapdata, vector<vector<OrbitData>> &orbitdata){

	unsigned long long descendantID = snapdata[snap].Halo[i].descendant;
	Int_t descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
	Int_t descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);
	unsigned long long haloID = snapdata[snap].Halo[i].id;
	Int_t halosnap = (Int_t)(haloID/opt.TEMPORALHALOIDVAL);
	Int_t haloindex = (Int_t)(haloID%opt.TEMPORALHALOIDVAL-1);


	//Keep track of this halo's snap and index for interpolation
	vector<Int_t> halosnaps;
	vector<Int_t> haloindexes;
	vector<Int_t> interpsnaps;

	//Store the index of the halo it is orbiting
	Int_t orbitinghaloindex;

	//Keep track of the previous halos halodata and orbitdata
	HaloData prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
	HaloData prevhosthalo;
	OrbitData prevorbitdata = {0};

	//Keep track of the snapshot
	Int_t currentsnap = snap;

	//Lets see if any interpolation needs to be done for this halo
	while(true){

		halosnaps.push_back(halosnap);
		haloindexes.push_back(haloindex);

		//Keep track of the snapshots that need to be interpolated and add to the orbitdata
		//for the halos that are to be interpolated
		while(currentsnap!=halosnap){
			interpsnaps.push_back(currentsnap);

			OrbitData tmporbitdata;
			orbitdata[currentsnap].push_back(tmporbitdata);

			currentsnap++;
		}

		//See if have reached the end of this branch
		if(descendantID==haloID) break;

		//Lets move to the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;


		//Extract its descendant
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);

		//Interate keeping track of the snapshots
		currentsnap++;
	}

	// If the interp snapshots contains snapshots then interpolation needs to be done
	if(interpsnaps.size()>0) InterpHaloProps(opt,halosnaps,haloindexes,interpsnaps,snapdata);

	//Reset the tree info to back at the base of the tree
	descendantID = snapdata[snap].Halo[i].descendant;
	descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
	descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);
	haloID = snapdata[snap].Halo[i].id;
	halosnap = (Int_t)(haloID/opt.TEMPORALHALOIDVAL);
	haloindex = (Int_t)(haloID%opt.TEMPORALHALOIDVAL-1);

	while(true){


		//Extract the halo it is orbiting at this snapshot
		orbitinghaloindex = (Int_t)(snapdata[halosnap].Halo[haloindex].orbitinghaloid%opt.TEMPORALHALOIDVAL-1);

		//Lets set this halos orbit data
		prevorbitdata = CalcOrbitProps(snapdata[halosnap].Halo[haloindex],snapdata[halosnap].Halo[orbitinghaloindex],prevorbitinghalo,prevhosthalo,prevorbitdata,snapdata[halosnap].scalefactor);
		prevhosthalo = snapdata[halosnap].Halo[orbitinghaloindex];
		orbitdata[halosnap][haloindex] = prevorbitdata;

		//Mark this halo as being done:
		snapdata[halosnap].Halo[haloindex].doneflag = true;

		//See if have reached the end of this branch
		if(descendantID==haloID) break;

		//Update the previous halo data
		prevorbitinghalo = snapdata[halosnap].Halo[haloindex];

		//lets move onto the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;

		//Extract its descendant
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);

	}


}

void ProcessOrbits(Options &opt, SnapData *&snapdata, vector<vector<OrbitData>> &orbitdata){

	int numsnaps =  opt.fsnap - opt.isnap+1;

	// Initilize the flag which marks the halo as being processed to false
	// and the vector to store the orbital data
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++){

		//Initilize a vector and reserve it (this intilizes the values to 0)
		vector<OrbitData> tmporbitdata;
		tmporbitdata.resize(snapdata[snap].numhalos);

		//Append the vector into the orbit data
		orbitdata.push_back(tmporbitdata);

		//Now go through all the halos marking their doneflag as false
		for(Int_t i=0;i<snapdata[snap].numhalos;i++){
			snapdata[snap].Halo[i].doneflag = false;

			//lets also convert to comoving distances if the distances are in physical units
			if(Units.distFlag==false){
				snapdata[snap].Halo[i].x*=Cosmo.h/snapdata[snap].scalefactor;
				snapdata[snap].Halo[i].y*=Cosmo.h/snapdata[snap].scalefactor;
				snapdata[snap].Halo[i].z*=Cosmo.h/snapdata[snap].scalefactor;
			}
		}
	}

	bool done = false;

	// Now lets start at the starting snapshot and walk up the tree
	// calculating the orbit relative to the halo which it was found
	// to be orbiting
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++){
		for(Int_t i=0;i<snapdata[snap].numhalos;i++){

			//Lets first check if this halo has been processed or is not orbiting a halo
			if((snapdata[snap].Halo[i].doneflag) | (snapdata[snap].Halo[i].orbitinghaloid==-1)) continue;

			ProcessHalo(snap,i,opt,snapdata,orbitdata);
			done = true;
			break;
			

		}
		if(opt.iverbose) cout<<"Done processing snap "<<snap<<endl;
		if(done) break;
	}
}