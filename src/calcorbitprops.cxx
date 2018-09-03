

#include "orbweaver.h"


void CalcOrbitProps(Int_t orbitID, int currentsnap, int prevsnap, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata, SnapData *&snapdata, OrbitProps &orbitprops){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*Cosmo.boxsize) orbitinghalo.x-=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r,vr,vrel;

	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = hosthalo.x - orbitinghalo.x;
	ry = hosthalo.y - orbitinghalo.y;
	rz = hosthalo.z - orbitinghalo.z;
	vrx = hosthalo.vx - orbitinghalo.vx;
	vry = hosthalo.vy - orbitinghalo.vy;
	vrz = hosthalo.vz - orbitinghalo.vz;
	r = sqrt(rx * rx + ry * ry + rz * rz);
	vr = (rx * vrx + ry * vry + rz * vrz) / r;
	vrel = sqrt(vrx*vrx + vry*vry + vrz*vrz);

	double mu, lx, ly, lz, e, count;
	// Find the average reduced mass, angular momentum and energy of this orbit
	// once the halo is found to be orbiting its host
	if(orbitprops.orbitingflag){

		//The halos reduced mass
		mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);
		orbitprops.mu += mu;

		//Angular momentum vectors
		lx = (ry * vrz) - (rz * vry);
		ly = -((rx * vrz) - (rz * vrx));
		lz = (rx * vry) - (ry * vrx);
		orbitprops.lx += lx;
		orbitprops.ly += ly;
		orbitprops.lz += lz;

		//The halos total orbital angular momentum
		orbitprops.ltot += mu * sqrt(lx*lx + ly*ly + lz*lz);

		//The halos orbital energy
		orbitprops.E += 0.5 * mu * vrel * vrel - (Cosmo.G * orbitinghalo.mass * hosthalo.mass)/r;
	}

	//Store the peak vmax in the orbiting halos history
	if(orbitinghalo.vmax>tmporbitdata.vmaxpeak)
		tmporbitdata.vmaxpeak = orbitinghalo.vmax;

	//Lets find if this halo has the closest approach so far or if we are at the base of this
	//branch i.e. the progenitor ID is the same as the halo's
	if(r<tmporbitdata.closestapproach)
		tmporbitdata.closestapproach = r;
	else if(orbitinghalo.progenitor==orbitinghalo.id)
		tmporbitdata.closestapproach = r;

	double prevrx,prevry,prevrz,prevvrx,prevvry,prevvrz,prevr,prevvr;

	//Lets find the same for the previous halo
	prevrx = prevhosthalo.x - prevorbitinghalo.x;
	prevry = prevhosthalo.y - prevorbitinghalo.y;
	prevrz = prevhosthalo.z - prevorbitinghalo.z;
	prevvrx = prevhosthalo.vx - prevorbitinghalo.vx;
	prevvry = prevhosthalo.vy - prevorbitinghalo.vy;
	prevvrz = prevhosthalo.vz - prevorbitinghalo.vz;
	prevr = sqrt(prevrx * prevrx + prevry * prevry + prevrz * prevrz);
	prevvr = (prevrx * prevvrx + prevry * prevvry + prevrz * prevvrz) / prevr;

	//Define varibles for the calculations
	double omega, deltat, ltot, E;


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

			//The difference in time since the previous snapshot
			deltat = snapdata[currentsnap].uniage - snapdata[prevsnap].uniage;

			// Find the change mass in units of Msun/Gyr
			tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/deltat;

			//Find the angular distance
			omega = acos((rx * prevrx + ry * prevry + rz * prevrz)/(r*prevr));

			//The tangential velocity of the orbiting halo with respect to its host
			tmporbitdata.vtan = r*(omega/deltat)*3.086e+19/3.15e+16;

			//Any additional properties to be calculated here

			//Set the orbit period and eccentricty as -1.0
			tmporbitdata.orbitperiod = -1.0;
			tmporbitdata.orbitecc = -1.0;
			tmporbitdata.lxrel = -1.0;
			tmporbitdata.lyrel = -1.0;
			tmporbitdata.lzrel = -1.0;

			//Now append it into the orbitdata dataset
			branchorbitdata.push_back(tmporbitdata);

			//If has undergone a crossing then there is no need to check the other crossings
			break;
		}
	}

	/* Check if the halo has gone past pericenter or apocenter */

	//Check if has undergone its first peri-centric passage (orbitingflag==true) and has undergone
	//a change in its radial motion. Otherwise if the orbitingflag==false then check if the halo
	//has had a pericentric passage within the host halos virial radius which then switches on
	//the orbiting flag so the number of orbits is tracked
	if((orbitprops.orbitingflag) & (vr*prevvr<0)){

		//Add 0.5 an orbit
		tmporbitdata.numorbits = tmporbitdata.numorbits + 0.5;

		/* Store some properties of the orbit halo and its host at this point */

		//Store what orbitID number this is
		tmporbitdata.orbitID = orbitID;

		//Set this as a passage point
		tmporbitdata.entrytype = 0.0;

		//The orbting halo
		tmporbitdata.haloID = orbitinghalo.origid;

		//The host halo
		tmporbitdata.hosthaloID = hosthalo.origid;

		//Store the scalefactor this happens at
		tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((vr/(vr - prevvr))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

		//From this scalefactor we can find the age of the universe
		tmporbitdata.uniage = GetUniverseAge(tmporbitdata.scalefactor);

		InterpPassageHaloProps(tmporbitdata.uniage,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata);

		/* Calculate various properties to be outputted */

		//Calculate the orbit period as 2x the previous passage
		tmporbitdata.orbitperiod = 2.0* (tmporbitdata.uniage - orbitprops.prevpassagetime);

		//The difference in time since the previous snapshot
		deltat = snapdata[currentsnap].uniage - snapdata[prevsnap].uniage;

		// Find the change mass in units of Msun/Gyr
		tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/deltat;

		//Find the average energy, total angular momentum and reduced mass
		E = orbitprops.E / (double)(currentsnap - orbitprops.prevpassagesnap);
		ltot = orbitprops.ltot / (double)(currentsnap - orbitprops.prevpassagesnap);
		mu = orbitprops.mu / (double)(currentsnap - orbitprops.prevpassagesnap);

		//The halos orbital eccentricity from the average properties
		tmporbitdata.orbitecc = sqrt(1.0 + ((2.0 * E * ltot*ltot)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu)));

		//Find the orbits apo and pericentric distances
		tmporbitdata.rperi = (ltot*ltot)/((1+tmporbitdata.orbitecc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

		//Find the orbits apo and pericentric distances
		tmporbitdata.rapo = (ltot*ltot)/((1-tmporbitdata.orbitecc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu); //1.0 / ((2.0/r) - (vr*vr)/(Cosmo.G*hosthalo.mass));

		//Find the angular distance
		omega = acos((rx * prevrx + ry * prevry + rz * prevrz)/(r*prevr));

		//The tangential velocity of the orbiting halo with respect to its host
		tmporbitdata.vtan = r*(omega/deltat)*3.086e+19/3.15e+16;

		//Lets output the average angular momentum since the last passage
		tmporbitdata.lxrel = orbitprops.lx / (double)(currentsnap - orbitprops.prevpassagesnap);
		tmporbitdata.lyrel = orbitprops.ly / (double)(currentsnap - orbitprops.prevpassagesnap);
		tmporbitdata.lzrel = orbitprops.lz / (double)(currentsnap - orbitprops.prevpassagesnap);

		//Reset the total angular momentum in the orbit props to zero
		orbitprops.lx = 0.0;
		orbitprops.ly = 0.0;
		orbitprops.lz = 0.0;
		orbitprops.E = 0.0;
		orbitprops.ltot = 0.0;
		orbitprops.mu = 0.0;

		//Any additional properties to be calculated here

		//Now append it into the orbitdata dataset
		branchorbitdata.push_back(tmporbitdata);

		//Update the previous passage time
		orbitprops.prevpassagetime = tmporbitdata.uniage;

		//Mark the snapshot that this passage happens
		orbitprops.prevpassagesnap = currentsnap;

		return;

	}
	else if((orbitprops.orbitingflag==false) & (vr*prevvr<0) & (r<3*hosthalo.rvir)){

		//Now lets mark this halo as on an orbit
		orbitprops.orbitingflag = true;

		/* Store some properties of the orbit halo and its host at this point */

		//Store what orbitID number this is
		tmporbitdata.orbitID = orbitID;

		//Set this as a passage point
		tmporbitdata.entrytype = 0.0;

		//Set this halo as done 1/2 orbit
		tmporbitdata.numorbits = 0.5;

		//The orbting halo
		tmporbitdata.haloID = orbitinghalo.origid;

		//The host halo
		tmporbitdata.hosthaloID = hosthalo.origid;

		//Store the scalefactor this happens at
		tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((prevvr/(prevvr - vr))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

		//From this scalefactor we can find the age of the universe
		tmporbitdata.uniage = GetUniverseAge(tmporbitdata.scalefactor);

		InterpPassageHaloProps(tmporbitdata.uniage,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata);

		/* Calculate various properties to be outputted */

		//Set the orbit period, angular momentum and eccentricty as -1.0 as these cannot be calculated yet
		tmporbitdata.orbitperiod = -1.0;
		tmporbitdata.orbitecc = -1.0;
		tmporbitdata.lxrel = -1.0;
		tmporbitdata.lyrel = -1.0;
		tmporbitdata.lzrel = -1.0;

		//The difference in time since the previous snapshot
		deltat = snapdata[currentsnap].uniage - snapdata[prevsnap].uniage;

		// Find the change mass in units of Msun/Gyr
		tmporbitdata.masslossrate = (orbitinghalo.mass - prevorbitinghalo.mass)/deltat;

		//Find the angular distance
		omega = acos((rx * prevrx + ry * prevry + rz * prevrz)/(r*prevr));

		//The tangential velocity of the orbiting halo with respect to its host
		tmporbitdata.vtan = r*(omega/deltat)*3.086e+19/3.15e+16;

		//Any additional properties to be calculated here

		//Now append it into the orbitdata dataset
		branchorbitdata.push_back(tmporbitdata);

		//Mark the time this passage happens
		orbitprops.prevpassagetime = tmporbitdata.uniage;

		//Mark the snapshot that this passage happens
		orbitprops.prevpassagesnap = currentsnap;
		return;
	}

}

double *computeAngles(double prevpos[3], OrbitData orbitdata){

	double x[3],y[3],z[3], mag;
	double* angles = new double[3];

	//Get the x-axis which is along the semi-major axis
	x[0] = orbitdata.xrel - prevpos[0];
	x[1] = orbitdata.yrel - prevpos[1];
	x[2] = orbitdata.zrel - prevpos[2];

	//The z-axis which is the angular momentum
	z[0] = orbitdata.lxrel;
	z[1] = orbitdata.lyrel;
	z[2] = orbitdata.lzrel;

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

void SetPassageType(vector<OrbitData> &branchorbitdata){

	double r,prevr=0,prevpos[3], *currangles, *refangles;
	int firstpassageindex,preventrytype=-2;
	vector<int> deleteindx;
	bool apofirst = false;

	//Now that the positions of the passages are set can now set the entry
	//type for these points and the orbit eccentricities at these points
	for(int i = 0; i<branchorbitdata.size();i++){

		//Only interested in passage points
		if(branchorbitdata[i].entrytype==0.0){

			//Find the radial distance to its host
			r = sqrt((branchorbitdata[i].xrel * branchorbitdata[i].xrel) + (branchorbitdata[i].yrel * branchorbitdata[i].yrel) + (branchorbitdata[i].zrel * branchorbitdata[i].zrel));

			//Can only do this if the previous radial distance has been set
			if(prevr>0){

				//Now we have the the semi-major axis vector
				if((refangles[0] + refangles[1] + refangles[2])==0.0)
					refangles = computeAngles(prevpos,branchorbitdata[i]);

				//Compute the current angles of the orbital plane
				currangles = computeAngles(prevpos,branchorbitdata[i]);

				//Compute the difference in the angles
				branchorbitdata[i].longascnode = currangles[0] - refangles[0];
				branchorbitdata[i].inc = currangles[1] - refangles[1];
				branchorbitdata[i].argpariap = currangles[2] - refangles[2];


				// If the current distance is greater that the previous radial distance then
				//this is a apo-centric passage, otherwise it is a pericentric passage
				if(r>prevr){

					//If the previous entry type was a apocentric, then delete this one and
					//wait until the halo goes past pericenter
					if((preventrytype==-1) | ((r-prevr)<0.05)){
						deleteindx.push_back(i);
						continue;
					}
					branchorbitdata[i].entrytype = -1;
					branchorbitdata[i].orbiteccratio = (r-prevr)/(r+prevr);

				}
				else{

					//If the previous entry type was a pericentric, then delete this one and
					//wait until the halo goes past apocenter
					if((preventrytype==0) | ((prevr-r)<0.05)){
						deleteindx.push_back(i);
						continue;
					}
					branchorbitdata[i].orbiteccratio = (prevr-r)/(prevr+r);

					//If the 2nd passage was a pericentric passage then set the initial passage as
					//a apocentric passage
					if(branchorbitdata[i].numorbits==1.0){
						cout<<firstpassageindex<<endl;
						branchorbitdata[firstpassageindex].entrytype = -1;
					}
				}
			}
			prevr = r;
			preventrytype = branchorbitdata[i].entrytype;
			prevpos[0] = branchorbitdata[i].xrel;
			prevpos[1] = branchorbitdata[i].yrel;
			prevpos[2] = branchorbitdata[i].zrel;

			//Store the index of the first passage and the relative position vector
			if(branchorbitdata[i].numorbits==0.5)
				firstpassageindex = i;

		}
	}

	for(int i =deleteindx.size()-1 ;i>=0; i--)
		branchorbitdata.erase(branchorbitdata.begin()+deleteindx[i]);
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
	HaloData prevorbitinghalo = {0};
	HaloData prevhosthalo = {0};
	vector<OrbitData> branchorbitdata;
	OrbitData tmporbitdata={0};
	Int_t prevsnap=halosnap;

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
		if(prevsnap!=halosnap)
			CalcOrbitProps(orbitID,halosnap,prevsnap,snapdata[halosnap].Halo[haloindex],snapdata[halosnap].Halo[orbitinghaloindex],prevorbitinghalo,prevhosthalo,branchorbitdata,tmporbitdata,snapdata,orbitprops);

		// if(find(interpsnaps.begin(), interpsnaps.end(), halosnap) != interpsnaps.end())
		// 	file<<-halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;
		// else
		// 	file<<halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;

		//Mark this halo as being done:
		snapdata[halosnap].Halo[haloindex].doneflag = true;

		//See if have reached the end of this branch
		if(descendantID==haloID) break;

		//Update the previous halo and its host
		prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
		prevhosthalo = snapdata[halosnap].Halo[orbitinghaloindex];
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
	vector<double> interpuniages;
	for(int i = 0; i<branchorbitdata.size();i++){
		//Only extract the scalfactor if this is a passage entry
		if(branchorbitdata[i].entrytype<=0){
			//Find the age of the universe
			interpuniages.push_back(GetUniverseAge(branchorbitdata[i].scalefactor));
		}
	}

	int nhalo = halosnaps.size();
	int ninterp = interpuniages.size();

	//Lets see if there is anything to interpolate but only if the halo exists
	//for more than 3 snapshots is it possible to interpolate
	if((ninterp>0) & (nhalo>3)) InterpPassagePoints(nhalo,ninterp,interpuniages,halosnaps,haloindexes,hostindexes,snapdata,branchorbitdata);

	//Now have set the distances for the passages then the entry types and ratio
	//eccentricities can be calculated
	SetPassageType(branchorbitdata);


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
	// Int_t snap = 55;
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++){
	// Int_t i = 991;
		for(Int_t i=0;i<snapdata[snap].numhalos;i++){

			//Lets first check if this halo has been processed or is not orbiting a halo
			if((snapdata[snap].Halo[i].doneflag) | (snapdata[snap].Halo[i].orbitinghaloid==-1)) continue;

			ProcessHalo(orbitID,snap,i,opt,snapdata,orbitdata);
			orbitID++;

		}
		if(opt.iverbose) cout<<"Done processing snap "<<snap<<endl;
	}
}