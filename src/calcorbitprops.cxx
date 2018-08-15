

#include "orbweaver.h"

void CalcOrbitProps(Int_t orbitID, int currentsnap, int prevsnap, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &orbitdata, OrbitData &tmporbitdata, SnapData *&snapdata, OrbitProps &orbitprops){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*Cosmo.boxsize){
		// cout<<"The position was "<<orbitinghalo.x<<endl;
		orbitinghalo.x-=Cosmo.boxsize;
		// cout<<"Corrected the position to "<<orbitinghalo.x<<endl;
	}

	if((orbitinghalo.y - hosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

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
	if(orbitinghalo.progenitor==orbitinghalo.id){
		tmporbitdata.closestapproach = r;
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

	//Define varibles for the calculations
	double mu, E, ecc;

	//Lets find if this halo has the closest approach so far
	if(r<tmporbitdata.closestapproach){
		tmporbitdata.closestapproach = r;
	}

	/* Now lets see if a new datapoint needs to be created if the halo has crossed through a interger number of rvir up to opt.numrvir */
	for(float i = 3.0;i>0.0;i-=0.5){

		// Less check to see if the previous halo was beyond the host Rvir and
		//if the current halo is within the host Rvir so it has infallen
		if(((prevr>i*prevhosthalo.rvir) & (r<i*hosthalo.rvir)) | ((prevr<i*prevhosthalo.rvir) & (r>i*hosthalo.rvir))){

			/* Store some properties of the orbit halo and its host at this point */

			//Store what orbitID number this is
			tmporbitdata.orbitID = orbitID;

			//Store how many rvir this entry is
			tmporbitdata.entrytype = i;

			//Store the scalefactor this happens at
			tmporbitdata.scalefactor = snapdata[currentsnap].scalefactor;

			//Set the orbit period as -1.0 here as only calculated at the passages
			tmporbitdata.orbitperiod = -1.0;

			//The orbting halo
			tmporbitdata.haloID = orbitinghalo.id;
			tmporbitdata.x = orbitinghalo.x;
			tmporbitdata.y = orbitinghalo.y;
			tmporbitdata.z = orbitinghalo.z;
			tmporbitdata.vx = orbitinghalo.vx;
			tmporbitdata.vy = orbitinghalo.vy;
			tmporbitdata.vz = orbitinghalo.vz;
			tmporbitdata.mass = orbitinghalo.mass;
			tmporbitdata.vmax = orbitinghalo.vmax;
			tmporbitdata.rmax = orbitinghalo.rmax;
			tmporbitdata.cnfw = orbitinghalo.cnfw;
			tmporbitdata.xrel = rx;
			tmporbitdata.yrel = ry;
			tmporbitdata.zrel = rz;
			tmporbitdata.vxrel = vrx;
			tmporbitdata.vyrel = vry;
			tmporbitdata.vzrel = vrz;

			//The host halo
			tmporbitdata.hosthaloID = hosthalo.id;
			tmporbitdata.rvirhost = hosthalo.rvir;
			tmporbitdata.masshost = hosthalo.mass;
			tmporbitdata.vmaxhost = hosthalo.vmax;
			tmporbitdata.rmaxhost = hosthalo.rmax;
			tmporbitdata.cnfwhost = hosthalo.cnfw;

			/* Calculate various properties to be outputted */

			// Find the change mass in units of Msun/Gyr
			tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

			//Any additional properties to be calculated here

			//Set the orbit period and eccentricty as -1.0 if the halo is not yet orbiting
			if(orbitprops.orbitingflag==false){
				tmporbitdata.orbitperiod = -1.0;
				tmporbitdata.orbitecc = -1.0;
				tmporbitdata.Lorbit = -1.0;
			}

			//Now append it into the orbitdata dataset
			orbitdata.push_back(tmporbitdata);

			//If has undergone a crossing then there is no need to check the other crossings
			return;
		}
	}

	/* Check if the halo has gone past pericenter or apocenter */

	//Check if has undergone its first peri-centric passage (orbitingflag==true) and has undergone
	//a change in its radial motion. Otherwise if the orbitingflag==false then check if the halo
	//has had a pericentric passage within the host halos virial radius which then switches on
	//the orbiting flag so the number of orbits is tracked
	if((orbitprops.orbitingflag) & (vr*prevvr<0)){
		// cout<<orbitinghalo.id<<" Gone 1/2 orbit "<<hosthalo.rvir<<endl;
		tmporbitdata.numorbits = tmporbitdata.numorbits + 0.5;

		//Check if lass passage was a peri-centric passage (passageflag==true)
		//so the next passgage will be apo-centric
		if(orbitprops.passageflag){

			/* Store some properties of the orbit halo and its host at this point */

			//Store what orbitID number this is
			tmporbitdata.orbitID = orbitID;

			//Mark this as a apo-centric passage
			tmporbitdata.entrytype = -1;

			//Store the scalefactor this happens at
			tmporbitdata.scalefactor = snapdata[currentsnap].scalefactor;

			//The orbting halo
			tmporbitdata.haloID = orbitinghalo.origid;
			tmporbitdata.x = orbitinghalo.x;
			tmporbitdata.y = orbitinghalo.y;
			tmporbitdata.z = orbitinghalo.z;
			tmporbitdata.vx = orbitinghalo.vx;
			tmporbitdata.vy = orbitinghalo.vy;
			tmporbitdata.vz = orbitinghalo.vz;
			tmporbitdata.mass = orbitinghalo.mass;
			tmporbitdata.vmax = orbitinghalo.vmax;
			tmporbitdata.rmax = orbitinghalo.rmax;
			tmporbitdata.cnfw = orbitinghalo.cnfw;
			tmporbitdata.xrel = rx;
			tmporbitdata.yrel = ry;
			tmporbitdata.zrel = rz;
			tmporbitdata.vxrel = vrx;
			tmporbitdata.vyrel = vry;
			tmporbitdata.vzrel = vrz;

			//The host halo
			tmporbitdata.hosthaloID = hosthalo.origid;
			tmporbitdata.rvirhost = hosthalo.rvir;
			tmporbitdata.masshost = hosthalo.mass;
			tmporbitdata.vmaxhost = hosthalo.vmax;
			tmporbitdata.rmaxhost = hosthalo.rmax;
			tmporbitdata.cnfwhost = hosthalo.cnfw;

			/* Calculate various properties to be outputted */

			// Find the change mass in units of Msun/Gyr
			tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

			//Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (snapdata[currentsnap].uniage - orbitprops.prevpassagetime);

			//The halos reduced mass
			mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);

			//The halos orbital energy
			E = 0.5 * mu * vr * vr - (Cosmo.G * orbitinghalo.mass * hosthalo.mass)/r;

			//The halos orbital angular momentum
			tmporbitdata.Lorbit = mu * r * vr;

			//The halos orbital eccentricity
			tmporbitdata.orbitecc = sqrt(1 + (2*E*tmporbitdata.Lorbit*tmporbitdata.Lorbit)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu));

			//Any additional properties to be calculated here

			//Now append it into the orbitdata dataset
			orbitdata.push_back(tmporbitdata);

			//Mark this passage as a apo-centric passage
			orbitprops.passageflag = false;

			//Update the previous passage time
			orbitprops.prevpassagetime = snapdata[currentsnap].uniage;

			//Keep track of the previous passage radial distance
			orbitprops.prevpassager = r;

			return;

		}
		else{

			/* Store some properties of the orbit halo and its host at this point */

			//Store what orbitID number this is
			tmporbitdata.orbitID = orbitID;

			//Mark this as a peri-centric passage
			tmporbitdata.entrytype = 0;

			//Store the scalefactor this happens at
			tmporbitdata.scalefactor = snapdata[currentsnap].scalefactor;

			//The orbting halo
			tmporbitdata.haloID = orbitinghalo.origid;
			tmporbitdata.x = orbitinghalo.x;
			tmporbitdata.y = orbitinghalo.y;
			tmporbitdata.z = orbitinghalo.z;
			tmporbitdata.vx = orbitinghalo.vx;
			tmporbitdata.vy = orbitinghalo.vy;
			tmporbitdata.vz = orbitinghalo.vz;
			tmporbitdata.mass = orbitinghalo.mass;
			tmporbitdata.vmax = orbitinghalo.vmax;
			tmporbitdata.rmax = orbitinghalo.rmax;
			tmporbitdata.cnfw = orbitinghalo.cnfw;
			tmporbitdata.xrel = rx;
			tmporbitdata.yrel = ry;
			tmporbitdata.zrel = rz;
			tmporbitdata.vxrel = vrx;
			tmporbitdata.vyrel = vry;
			tmporbitdata.vzrel = vrz;

			//The host halo
			tmporbitdata.hosthaloID = hosthalo.origid;
			tmporbitdata.rvirhost = hosthalo.rvir;
			tmporbitdata.masshost = hosthalo.mass;
			tmporbitdata.vmaxhost = hosthalo.vmax;
			tmporbitdata.rmaxhost = hosthalo.rmax;
			tmporbitdata.cnfwhost = hosthalo.cnfw;

			/* Calculate various properties to be outputted */

			//Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (snapdata[currentsnap].uniage - orbitprops.prevpassagetime);

			// Find the change mass in units of Msun/Gyr
			tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

			//The halos reduced mass
			mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);

			//The halos orbital energy
			E = 0.5 * mu * vr * vr - (Cosmo.G * orbitinghalo.mass * hosthalo.mass)/r;

			//The halos orbital angular momentum
			tmporbitdata.Lorbit = mu * r * vr;

			//The halos orbital eccentricity
			tmporbitdata.orbitecc = sqrt(1 + (2*E*tmporbitdata.Lorbit*tmporbitdata.Lorbit)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu));

			//Any additional properties to be calculated here

			//Now append it into the orbitdata dataset
			orbitdata.push_back(tmporbitdata);

			//Mark this passage as a peri-centric passage
			orbitprops.passageflag = true;

			//Update the previous passage time
			orbitprops.prevpassagetime = snapdata[currentsnap].uniage;

			//Keep track of the previous passage radial distance
			orbitprops.prevpassager = r;

			return;
		}
	}
	else if((orbitprops.orbitingflag==false) & (vr*prevvr<0) & (r<hosthalo.rvir)){
		// cout<<orbitinghalo.id<<" This halo has started to orbit "<<hosthalo.rvir<<endl;

		//Now lets mark this halo as on an orbit
		orbitprops.orbitingflag = true;

		/* Store some properties of the orbit halo and its host at this point */

		//Store what orbitID number this is
		tmporbitdata.orbitID = orbitID;

		//Mark this as a peri-centric passage
		tmporbitdata.entrytype = 0;

		//Store the scalefactor this happens at
		tmporbitdata.scalefactor = snapdata[currentsnap].scalefactor;

		//Mark this orbit as having a period of -1.0 as it cannot be calculated yet
		tmporbitdata.orbitperiod = -1.0;

		//The orbting halo
		tmporbitdata.haloID = orbitinghalo.origid;
		tmporbitdata.x = orbitinghalo.x;
		tmporbitdata.y = orbitinghalo.y;
		tmporbitdata.z = orbitinghalo.z;
		tmporbitdata.vx = orbitinghalo.vx;
		tmporbitdata.vy = orbitinghalo.vy;
		tmporbitdata.vz = orbitinghalo.vz;
		tmporbitdata.mass = orbitinghalo.mass;
		tmporbitdata.vmax = orbitinghalo.vmax;
		tmporbitdata.rmax = orbitinghalo.rmax;
		tmporbitdata.cnfw = orbitinghalo.cnfw;
		tmporbitdata.xrel = rx;
		tmporbitdata.yrel = ry;
		tmporbitdata.zrel = rz;
		tmporbitdata.vxrel = vrx;
		tmporbitdata.vyrel = vry;
		tmporbitdata.vzrel = vrz;

		//The host halo
		tmporbitdata.hosthaloID = hosthalo.origid;
		tmporbitdata.rvirhost = hosthalo.rvir;
		tmporbitdata.masshost = hosthalo.mass;
		tmporbitdata.vmaxhost = hosthalo.vmax;
		tmporbitdata.rmaxhost = hosthalo.rmax;
		tmporbitdata.cnfwhost = hosthalo.cnfw;

		/* Calculate various properties to be outputted */

		//Set the orbit period, angular momentum and eccentricty as -1.0 as these cannot be calculated yet
		tmporbitdata.orbitperiod = -1.0;
		tmporbitdata.orbitecc = -1.0;
		tmporbitdata.Lorbit = -1.0;

		// Find the change mass in units of Msun/Gyr
		tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

		//Any additional properties to be calculated here

		//Now append it into the orbitdata dataset
		orbitdata.push_back(tmporbitdata);

		//Mark this passage a pericentric passage
		orbitprops.passageflag = true;

		//Mark the time this passage happens
		orbitprops.prevpassagetime = snapdata[currentsnap].uniage;

		//Keep track of the previous passage radial distance
		orbitprops.prevpassager = r;

		return;
	}

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


double LogInterp(double prevdata, double nextdata, double f){
	return pow(nextdata,f) * pow(prevdata,1-f);
}
double LinInterp(double prevdata, double nextdata, double f){
	return prevdata + (nextdata - prevdata)*f;
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
		snapdata[interpsnaps[i]].numhalos++;
	}


}


void ProcessHalo(Int_t orbitID,Int_t snap, Int_t i, Options &opt, SnapData *&snapdata, vector<OrbitData> &orbitdata){

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
	OrbitData tmporbitdata = {0};
	Int_t prevsnap=halosnap-1;

	//Keep track of the properties of this orbit
	OrbitProps orbitprops;

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

	//If the interp snapshots contains snapshots then interpolation needs to be done
	//but only if the halo exists for more than 3 snapshots is it possible
	if((interpsnaps.size()>0) & (halosnaps.size()>3)) InterpHaloProps(opt,halosnaps,haloindexes,interpsnaps,snapdata);

	//Reset the tree info to back at the base of the tree
	descendantID = snapdata[snap].Halo[i].descendant;
	descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
	descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);
	haloID = snapdata[snap].Halo[i].id;
	halosnap = (Int_t)(haloID/opt.TEMPORALHALOIDVAL);
	haloindex = (Int_t)(haloID%opt.TEMPORALHALOIDVAL-1);

	// ofstream file;
	// file.open("../analysis/data/lininterp.dat");

	while(true){

		//Extract the halo it is orbiting at this snapshot
		orbitinghaloindex = (Int_t)(snapdata[halosnap].Halo[haloindex].orbitinghaloid%opt.TEMPORALHALOIDVAL-1);

		//Lets set this halos orbit data
		CalcOrbitProps(orbitID,halosnap,prevsnap,snapdata[halosnap].Halo[haloindex],snapdata[halosnap].Halo[orbitinghaloindex],prevorbitinghalo,prevhosthalo,orbitdata,tmporbitdata,snapdata,orbitprops);
		prevhosthalo = snapdata[halosnap].Halo[orbitinghaloindex];

		// if(find(interpsnaps.begin(), interpsnaps.end(), halosnap) != interpsnaps.end())
		// 	file<<-halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;
		// else
		// 	file<<halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;

		//Mark this halo as being done:
		snapdata[halosnap].Halo[haloindex].doneflag = true;

		//See if have reached the end of this branch
		if(descendantID==haloID) break;

		//Update the previous halo data and snapshot
		prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
		prevsnap = halosnap;

		//lets move onto the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;

		//Extract its descendant
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);

	}
	// file.close();
}

void ProcessOrbits(Options &opt, SnapData *&snapdata, vector<OrbitData> &orbitdata){

	int numsnaps =  opt.fsnap - opt.isnap+1;

	// Initilize the flag which marks the halo as being processed to false
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++)
		for(Int_t i=0;i<snapdata[snap].numhalos;i++)
			snapdata[snap].Halo[i].doneflag = false;


	Int_t orbitID = 0;

	// Now lets start at the starting snapshot and walk up the tree
	// calculating the orbit relative to the halo which it was found
	// to be orbiting
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++){
		for(Int_t i=0;i<snapdata[snap].numhalos;i++){

			//Lets first check if this halo has been processed or is not orbiting a halo
			if((snapdata[snap].Halo[i].doneflag) | (snapdata[snap].Halo[i].orbitinghaloid==-1)) continue;

			ProcessHalo(orbitID,snap,i,opt,snapdata,orbitdata);
			orbitID++;

		}
		if(opt.iverbose) cout<<"Done processing snap "<<snap<<endl;
	}
}