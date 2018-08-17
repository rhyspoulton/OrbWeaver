

#include "orbweaver.h"


void CalcOrbitProps(Int_t orbitID, int currentsnap, int prevsnap, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata, SnapData *&snapdata, OrbitProps &orbitprops){

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
	vr = (rx * vrx + ry * vry + rz * vrz) / r;

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
	prevvr = (prevrx * prevvrx + prevry * prevvry + prevrz * prevvrz) / r;

	//Define varibles for the calculations
	double mu, E, ecc, currentuniage;

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
			branchorbitdata.push_back(tmporbitdata);

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

			//The orbting halo
			tmporbitdata.haloID = orbitinghalo.origid;

			//The host halo
			tmporbitdata.hosthaloID = hosthalo.origid;

			//Store the scalefactor this happens at
			tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((vr/(vr - prevvr))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

			//From this scalefactor we can find the age of the universe
			currentuniage = GetUniverseAge(tmporbitdata.scalefactor);

			InterpPassageHaloProps(currentuniage,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata);

			/* Calculate various properties to be outputted */

			// Find the change mass in units of Msun/Gyr
			tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

			//Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (currentuniage - orbitprops.prevpassagetime);

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
			branchorbitdata.push_back(tmporbitdata);

			//Mark this passage as a apo-centric passage
			orbitprops.passageflag = false;

			//Update the previous passage time
			orbitprops.prevpassagetime = currentuniage;

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

			//The orbting halo
			tmporbitdata.haloID = orbitinghalo.origid;

			//The host halo
			tmporbitdata.hosthaloID = hosthalo.origid;

			//Store the scalefactor this happens at
			tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((vr/(vr - prevvr))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

			//From this scalefactor we can find the age of the universe
			currentuniage = GetUniverseAge(tmporbitdata.scalefactor);

			InterpPassageHaloProps(currentuniage,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata);

			/* Calculate various properties to be outputted */

			//Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (currentuniage - orbitprops.prevpassagetime);

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
			branchorbitdata.push_back(tmporbitdata);

			//Mark this passage as a peri-centric passage
			orbitprops.passageflag = true;

			//Update the previous passage time
			orbitprops.prevpassagetime = currentuniage;

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

		//Set this halo as done 1/2 orbit
		tmporbitdata.numorbits = 0.5;

		//Mark this as a peri-centric passage
		tmporbitdata.entrytype = 0;

		//The orbting halo
		tmporbitdata.haloID = orbitinghalo.origid;

		//The host halo
		tmporbitdata.hosthaloID = hosthalo.origid;

		//Store the scalefactor this happens at
		tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((vr/(vr - prevvr))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

		//From this scalefactor we can find the age of the universe
		currentuniage = GetUniverseAge(tmporbitdata.scalefactor);

		InterpPassageHaloProps(currentuniage,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata);

		/* Calculate various properties to be outputted */

		//Set the orbit period, angular momentum and eccentricty as -1.0 as these cannot be calculated yet
		tmporbitdata.orbitperiod = -1.0;
		tmporbitdata.orbitecc = -1.0;
		tmporbitdata.Lorbit = -1.0;

		// Find the change mass in units of Msun/Gyr
		tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/(snapdata[currentsnap].uniage - snapdata[prevsnap].uniage);

		//Any additional properties to be calculated here

		//Now append it into the orbitdata dataset
		branchorbitdata.push_back(tmporbitdata);

		//Mark this passage a pericentric passage
		orbitprops.passageflag = true;

		//Mark the time this passage happens
		orbitprops.prevpassagetime = snapdata[currentsnap].uniage;

		//Keep track of the previous passage radial distance
		orbitprops.prevpassager = r;

		return;
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
	vector<Int_t> hostindexes;
	vector<Int_t> interpsnaps;

	//Store the index of the halo it is orbiting
	Int_t orbitinghaloindex;

	//Keep track of the previous halos halodata and orbitdata
	HaloData prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
	HaloData prevhosthalo;
	vector<OrbitData> branchorbitdata;
	OrbitData tmporbitdata={0};
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
	// file.open("../analysis/data/circ.dat");

	while(true){

		//Extract the halo it is orbiting at this snapshot
		orbitinghaloindex = (Int_t)(snapdata[halosnap].Halo[haloindex].orbitinghaloid%opt.TEMPORALHALOIDVAL-1);
		hostindexes.push_back(orbitinghaloindex);

		//Lets set this halos orbit data
		CalcOrbitProps(orbitID,halosnap,prevsnap,snapdata[halosnap].Halo[haloindex],snapdata[halosnap].Halo[orbitinghaloindex],prevorbitinghalo,prevhosthalo,branchorbitdata,tmporbitdata,snapdata,orbitprops);
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

	//Now another interpolation can be done to find the actual positions of the passage points
	InterpPassagePoints(halosnaps,haloindexes,hostindexes,snapdata,branchorbitdata);



	//Now finished with this branches orbital calculations so it can be added
	//into the orbitdata vector that contains all halos, it only needs to be
	//moved rather than copied
	orbitdata.insert(orbitdata.end(),make_move_iterator(branchorbitdata.begin()),make_move_iterator(branchorbitdata.end()));
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