

#include "orbweaver.h"

void CalcOrbitProps(Options &opt,
	unsigned long long orbitID,
	int currentsnap, int prevsnap,
	HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo,
	vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata,
	vector<SnapData> &snapdata,
	int* num_entrytypes, OrbitProps &orbitprops, OrbitProps &prevorbitprops,
	SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*Cosmo.boxsize) orbitinghalo.x-=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

	if((prevorbitinghalo.x - prevhosthalo.x)>0.5*Cosmo.boxsize) orbitinghalo.x-=Cosmo.boxsize;
	if((prevorbitinghalo.y - prevhosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((prevorbitinghalo.z - prevhosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((prevorbitinghalo.x - prevhosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((prevorbitinghalo.y - prevhosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((prevorbitinghalo.z - prevhosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r,vrad,vrel;

	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = hosthalo.x - orbitinghalo.x;
	ry = hosthalo.y - orbitinghalo.y;
	rz = hosthalo.z - orbitinghalo.z;
	vrx = hosthalo.vx - orbitinghalo.vx;
	vry = hosthalo.vy - orbitinghalo.vy;
	vrz = hosthalo.vz - orbitinghalo.vz;
	r = sqrt(rx * rx + ry * ry + rz * rz);
	vrad = (rx * vrx + ry * vry + rz * vrz) / r;
	vrel = sqrt(vrx*vrx + vry*vry + vrz*vrz);

	//Store the peak vmax in the orbiting halos history
	if(orbitinghalo.vmax>tmporbitdata.vmaxpeak)
		tmporbitdata.vmaxpeak = orbitinghalo.vmax;

	//Store the minimum minrmax and rscale so far
	if(orbitinghalo.rmax<orbitprops.minrmax)
		orbitprops.minrmax = orbitinghalo.rmax;
	if((orbitinghalo.rvir/orbitinghalo.cnfw)<orbitprops.minrscale)
		orbitprops.minrscale = orbitinghalo.rvir/orbitinghalo.cnfw;

	//Now done the interpolation can check if this is the closest approach so far
	if(r<orbitprops.closestapproach)
		orbitprops.closestapproach=r;
	tmporbitdata.closestapproach=orbitprops.closestapproach;

	//Store what orbitID number this is
	tmporbitdata.orbitID = orbitID;

	//The orbting halo
	tmporbitdata.haloID = orbitinghalo.origid;

	//The original Root Progenitor for the orbiting halo
	tmporbitdata.halorootprogenID = orbitinghalo.origrootprogenitor;

	//The host halo
	tmporbitdata.orbitedhaloID = hosthalo.origid;

	//The original Root Progenitor of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootprogenID = hosthalo.origrootprogenitor;

	//Store the haloID this halo if it is not interpolated
	if(orbitinghalo.interpflag)
		tmporbitdata.orbithaloID = 0;
	else
		tmporbitdata.orbithaloID = orbitinghalo.id;

	double prevrx,prevry,prevrz,prevvrx,prevvry,prevvrz,prevr,prevvrad;

	//Lets find the same for the previous halo
	prevrx = prevhosthalo.x - prevorbitinghalo.x;
	prevry = prevhosthalo.y - prevorbitinghalo.y;
	prevrz = prevhosthalo.z - prevorbitinghalo.z;
	prevvrx = prevhosthalo.vx - prevorbitinghalo.vx;
	prevvry = prevhosthalo.vy - prevorbitinghalo.vy;
	prevvrz = prevhosthalo.vz - prevorbitinghalo.vz;
	prevr = sqrt(prevrx * prevrx + prevry * prevry + prevrz * prevrz);
	prevvrad = (prevrx * prevvrx + prevry * prevvry + prevrz * prevvrz) / prevr;


	double mu, lx, ly, lz, e, count, deltat;
	//The difference in time since the previous snapshot
	deltat = snapdata[currentsnap].uniage - snapdata[prevsnap].uniage;

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
		orbitprops.lx += mu * lx;
		orbitprops.ly += mu * ly;
		orbitprops.lz += mu * lz;

		//The halos total orbital angular momentum
		orbitprops.ltot += sqrt(lx*lx + ly*ly + lz*lz);

		//The halos orbital energy
		orbitprops.E += 0.5 * mu * vrel * vrel - (Cosmo.G * orbitinghalo.mass * hosthalo.mass)/r;

		//Keep track of the host's angular momentum vectors
		orbitprops.hostlx += hosthalo.lx;
		orbitprops.hostly += hosthalo.ly;
		orbitprops.hostlz += hosthalo.lz;

		// Find the change mass in units of Msun/Gyr
		orbitprops.masslossrate += (prevorbitinghalo.mass - orbitinghalo.mass)/deltat;

		//Find the total angle that the orbit has moved through since last orbit
		orbitprops.phi+=acos((rx*prevrx + ry*prevry + rz*prevrz)/(r*prevr));
	}

	//Define varibles for the calculations
	double ltot, f, vcirc, vcomp, vradx, vrady, vradz, vtanx, vtany, vtanz, prevpassager, semiMajor, keplarperi_wetzel2011od, *currangles;
	int prevpassageindex;


	/* Now lets see if a new datapoint needs to be created if the halo has crossed through a interger number of rvir up to opt.numRvirSearch */
	float numrvircrossing=0;
	int entrytypeindex=0;
	int j = 0;
	for(float i = opt.numRvirSearch;i>0.0;i-=opt.fracrvircross){

		// Less check to see if the previous halo was beyond the host Rvir and
		//if the current halo is within the host Rvir so it has infallen or if it
		//is the other way around so has outfallen
		if((prevr>i*prevhosthalo.rvir) & (r<i*hosthalo.rvir)){
			numrvircrossing = i;

			entrytypeindex = j + 2;
			num_entrytypes[entrytypeindex]++;

			//Prioritize first infall when crossing rvir since this will set the merger timescale
			if(((round(i *10)/10)==1.0) & (orbitprops.crossrvirtime==0.0)) 	break;
		}
		else if((prevr<i*prevhosthalo.rvir) & (r>i*hosthalo.rvir)){
			numrvircrossing = -i;

			entrytypeindex = opt.numtypeofcrossingentries + j + 2;
			num_entrytypes[entrytypeindex]++;
		}
		j++;
	}

	//Keep track if this halo and its host is top of spatial herachy
	tmporbitdata.fieldhalo = orbitinghalo.fieldhalo;
	tmporbitdata.fieldhalohost = hosthalo.fieldhalo;

	//Set a flag if the orbit host merges with another halo
	tmporbitdata.hostmerges = orbitinghalo.hostmerges;

	//Store the information about it and its host substructure
	tmporbitdata.numsubstruct = orbitinghalo.numsubstruct;
	tmporbitdata.ratioofmassinsubsstruct = orbitinghalo.ratioofmassinsubsstruct;
	tmporbitdata.numsubstructhost = hosthalo.numsubstruct;
	tmporbitdata.ratioofmassinsubsstructhost = hosthalo.ratioofmassinsubsstruct;

	//Compute the instaneous mass loss rate

	// Find the change mass in units of Msun/Gyr
	tmporbitdata.masslossrate_inst = (prevorbitinghalo.mass - orbitinghalo.mass )/deltat;

	//Check if a crossing point has happened and it is not the same as the previous crossing point. In addition check if the same crossing point happend 
	if((numrvircrossing!=0) & (numrvircrossing!=orbitprops.prevcrossingentrytype) & ((abs(numrvircrossing)!=abs(orbitprops.prevcrossingentrytype)) | (currentsnap>orbitprops.prevcrossingsnap+1))){

		/* Store some properties of the orbit halo and its host at this point */

		// cout<<"Before "<<numrvircrossing<<" "<<r/hosthalo.rvir<<endl;

		//Store how many rvir this entry is
		tmporbitdata.entrytype = numrvircrossing;

		//The number of this type of entry
		tmporbitdata.num_entrytype = num_entrytypes[entrytypeindex];

		//The current number of orbits
		tmporbitdata.numorbits = orbitprops.numorbits;

		tmporbitdata.uniage = InterpCrossingHaloProps(numrvircrossing,snapdata[currentsnap].uniage,snapdata[prevsnap].uniage,orbitinghalo,hosthalo,prevorbitinghalo,prevhosthalo,tmporbitdata,snapdata,splinefuncs,hostsplinefuncs);

		tmporbitdata.scalefactor = GetScaleFactor(tmporbitdata.uniage);

		//Update all the quantities for calculations
		rx = tmporbitdata.xrel;
		ry = tmporbitdata.yrel;
		rz = tmporbitdata.zrel;
		r = sqrt(rx*rx + ry*ry + rz*rz);
		vrx = tmporbitdata.vxrel;
		vry = tmporbitdata.vyrel;
		vrz = tmporbitdata.vzrel;

		//Set the orbit period as -1.0 here as only calculated at the passages
		tmporbitdata.orbitperiod = -1.0;

		/* Calculate various properties to be outputted */

		//Find the components of the radial vector
		vcomp = (rx*vrx + ry*vry + rz*vrz)/(r*r);
		vradx = vcomp * rx;
		vrady = vcomp * ry;
		vradz = vcomp * rz;
		tmporbitdata.vrad = vcomp;

		//Then can use these components to find the components of the tangential velocity
		vtanx = vrx - vradx;
		vtany = vry - vrady;
		vtanz = vrz - vradz;
		tmporbitdata.vtan = sqrt(vtanx*vtanx + vtany*vtany + vtanz*vtanz);

		//The halos reduced mass
		mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);

		//Angular momentum vectors
		lx = (ry * vrz) - (rz * vry);
		ly = -((rx * vrz) - (rz * vrx));
		lz = (rx * vry) - (ry * vrx);
		tmporbitdata.lxrel_inst += mu * lx;
		tmporbitdata.lyrel_inst += mu * ly;
		tmporbitdata.lzrel_inst += mu * lz;


		//Find the energy of the orbit
		tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * orbitinghalo.mass*hosthalo.mass)/r;

		//Find the radius of a circular orbit with the same energy
		tmporbitdata.rcirc = abs((Cosmo.G * orbitinghalo.mass * Menc(r,hosthalo.mass,hosthalo.cnfw,hosthalo.rvir))/(2.0 * tmporbitdata.orbitalenergy_inst));

		//Find the circular velocity of a circular orbit with the same energy
		vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

		//Can use these to find the orbital angular momentum of a circular orbit
		tmporbitdata.jcirc = mu * tmporbitdata.rcirc * vcirc;

		//Any additional properties to be calculated here

		//Set the orbit period and eccentricty as -1.0
		tmporbitdata.orbitperiod = -1.0;
		tmporbitdata.orbitecc_wetzel2011 = -1.0;
		tmporbitdata.lxrel_ave = -1.0;
		tmporbitdata.lyrel_ave = -1.0;
		tmporbitdata.lzrel_ave = -1.0;
		tmporbitdata.hostalignment = 0.0;

		//Store the previous crossing point entry type
		orbitprops.prevcrossingentrytype=numrvircrossing;

		//Store the snapshot which this crossing point happened
		orbitprops.prevcrossingsnap = currentsnap;

		//Now done the interpolation can check if this is the closest approach so far
		r = sqrt(tmporbitdata.xrel*tmporbitdata.xrel + tmporbitdata.yrel*tmporbitdata.yrel + tmporbitdata.zrel*tmporbitdata.zrel);
		if(r<orbitprops.closestapproach)
			orbitprops.closestapproach=r;
		tmporbitdata.closestapproach=orbitprops.closestapproach;

		//Now append it into the orbitdata dataset
		branchorbitdata.push_back(tmporbitdata);

		//Set the time when this halo first crosses 1 rvir of its host halo
		if(((round(numrvircrossing*10.0)/10.0)==1.0) & (orbitprops.crossrvirtime==0.0))
			orbitprops.crossrvirtime = tmporbitdata.uniage;

	}


	/* Check if the halo has gone past pericenter or apocenter */

	//Check if has undergone its first peri-centric passage (orbitingflag==true) and has undergone
	//a change in its radial motion. Otherwise if the orbitingflag==false then check if the halo
	//has had a pericentric passage within the host halos virial radius which then switches on
	//the orbiting flag so the number of orbits is tracked
	if((vrad*prevvrad<0) & (r<3.0*hosthalo.rvir)){

		//Add 0.5 an orbit
		orbitprops.numorbits = orbitprops.numorbits + 0.5;
		tmporbitdata.numorbits = orbitprops.numorbits;

		/* Store some properties of the orbit halo and its host at this point */

		//Set this as a passage point dependant if it is outgoing or infalling
		if(vrad>0){ //Peri-centric
			tmporbitdata.entrytype = 99;
			num_entrytypes[0]++;
			tmporbitdata.num_entrytype = num_entrytypes[0];
		}
		else{    //Apocentric
			tmporbitdata.entrytype = -99;
			num_entrytypes[1]++;
			tmporbitdata.num_entrytype = num_entrytypes[1];
		}

		//If the previous passage is the same entry type then remove it this can happen if the halo went outside of 3 Rvir host and underwent a passgae
		if((orbitprops.prevpassageentrytype==tmporbitdata.entrytype) & (orbitprops.orbitingflag)){
			branchorbitdata.erase(branchorbitdata.begin()+orbitprops.prevpassageindex);
			orbitprops.numorbits= orbitprops.numorbits - 0.5;
			tmporbitdata.numorbits=orbitprops.numorbits;

			if(orbitprops.numorbits<=0.5){
				orbitprops.orbitingflag=false;
			}
			else{

				prevorbitprops.lx += orbitprops.lx;
				prevorbitprops.ly += orbitprops.ly;
				prevorbitprops.lz += orbitprops.lz;
				prevorbitprops.E += orbitprops.E;
				prevorbitprops.ltot += orbitprops.ltot;
				prevorbitprops.mu += orbitprops.mu;
				prevorbitprops.hostlx += orbitprops.hostlx;
				prevorbitprops.hostly += orbitprops.hostlx;
				prevorbitprops.hostlz += orbitprops.hostlx;
				prevorbitprops.masslossrate += orbitprops.masslossrate;
				prevorbitprops.phi += orbitprops.phi;

				 orbitprops = prevorbitprops;
			}
		}

		//Store the scalefactor this happens at
		tmporbitdata.scalefactor = exp(log(snapdata[currentsnap].scalefactor) -abs((vrad/(vrad - prevvrad))) * (log(snapdata[currentsnap].scalefactor/snapdata[prevsnap].scalefactor)));

		//From this scalefactor we can find the age of the universe
		tmporbitdata.uniage = GetUniverseAge(tmporbitdata.scalefactor);

		InterpSingleHaloProps(tmporbitdata.uniage, snapdata[currentsnap].uniage,snapdata[prevsnap].uniage, orbitinghalo, hosthalo, prevorbitinghalo, prevhosthalo, tmporbitdata, snapdata, splinefuncs, hostsplinefuncs);

		//Update all the quantities for calculations
		rx = tmporbitdata.xrel;
		ry = tmporbitdata.yrel;
		rz = tmporbitdata.zrel;
		r = sqrt(rx*rx + ry*ry + rz*rz);
		vrx = tmporbitdata.vxrel;
		vry = tmporbitdata.vyrel;
		vrz = tmporbitdata.vzrel;

		//The difference in time since the previous snapshot
		deltat = tmporbitdata.uniage - snapdata[prevsnap].uniage;

		//Find the components of the radial vector
		vcomp = (rx*vrx + ry*vry + rz*vrz)/(r*r);
		vradx = vcomp * rx;
		vrady = vcomp * ry;
		vradz = vcomp * rz;
		tmporbitdata.vrad = sqrt(vradx*vradx + vrady*vrady + vradz*vradz);

		//Then can use these components to find the components of the tangential velocity
		vtanx = vrx - vradx;
		vtany = vry - vrady;
		vtanz = vrz - vradz;
		tmporbitdata.vtan = sqrt(vtanx*vtanx + vtany*vtany + vtanz*vtanz);

		//The halos reduced mass
		mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);

		//Angular momentum vectors
		lx = (ry * vrz) - (rz * vry);
		ly = -((rx * vrz) - (rz * vrx));
		lz = (rx * vry) - (ry * vrx);
		tmporbitdata.lxrel_inst += mu * lx;
		tmporbitdata.lyrel_inst += mu * ly;
		tmporbitdata.lzrel_inst += mu * lz;

		//Find the energy of the orbit
		tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * orbitinghalo.mass*hosthalo.mass)/r;

		//Find the radius of a circular orbit with the same energy
		tmporbitdata.rcirc = abs((Cosmo.G * orbitinghalo.mass * Menc(r,hosthalo.mass,hosthalo.cnfw,hosthalo.rvir))/(2.0 * tmporbitdata.orbitalenergy_inst));

		//Find the circular velocity of a circular orbit with the same energy
		vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

		//Can use these to find the orbital angular momentum of a circular orbit
		tmporbitdata.jcirc = mu * tmporbitdata.rcirc * vcirc;


		/* Calculate various properties to be outputted if the halo is marked as orbiting */

		if(orbitprops.orbitingflag){

			prevpassager = sqrt(orbitprops.prevpassagepos[0]*orbitprops.prevpassagepos[0] + orbitprops.prevpassagepos[1]*orbitprops.prevpassagepos[1] + orbitprops.prevpassagepos[2]*orbitprops.prevpassagepos[2]);

			if(vrad>0)
				tmporbitdata.orbiteccratio = (prevpassager-r)/(prevpassager+r);
			else
				tmporbitdata.orbiteccratio = (r-prevpassager)/(r+prevpassager);

			// omega = acos((rx * orbitprops.prevpassagepos[0] + ry * orbitprops.prevpassagepos[1] + rz * orbitprops.prevpassagepos[2])/(r * sqrt(orbitprops.prevpassagepos[0]*orbitprops.prevpassagepos[0] + orbitprops.prevpassagepos[1]*orbitprops.prevpassagepos[1] + orbitprops.prevpassagepos[2]*orbitprops.prevpassagepos[2])));

			// //Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (tmporbitdata.uniage - orbitprops.prevpassagetime);

			// semiMajor = (prevpassager+r)/2.0; //sqrt((rx-orbitprops.prevpassagepos[0])*(rx-orbitprops.prevpassagepos[0]) + (ry-orbitprops.prevpassagepos[1])*(ry-orbitprops.prevpassagepos[1]) + (rz-orbitprops.prevpassagepos[2])*(rz-orbitprops.prevpassagepos[2]));

			// keplarperi_wetzel2011od = 2*3.142 * sqrt((semiMajor*semiMajor*semiMajor)/(Cosmo.G * (hosthalo.mass + orbitinghalo.mass))) *(3.086e+19/3.15e+16);

			//Remove any passages 
			if((currentsnap - orbitprops.prevpassagesnap) < 3){

				//Add the previous passage to the remove indexes
				branchorbitdata.erase(branchorbitdata.begin()+orbitprops.prevpassageindex);


				// If the object has only undergone one orbit then lets remove the passages and reset so the halo is 
				// no longer set to be orbiting, otherwise update the previous passage quantities
				if(tmporbitdata.numorbits==1.0){
					orbitprops.orbitingflag=false;

					//Reset the total angular momentum in the orbit props to zero
					orbitprops.lx = 0.0;
					orbitprops.ly = 0.0;
					orbitprops.lz = 0.0;
					orbitprops.E = 0.0;
					orbitprops.ltot = 0.0;
					orbitprops.mu = 0.0;
					orbitprops.hostlx = 0.0;
					orbitprops.hostly = 0.0;
					orbitprops.hostlz = 0.0;
					orbitprops.masslossrate = 0.0;
					orbitprops.phi=0.0;
					orbitprops.numorbits=0.0;
				}
				else{

					prevorbitprops.lx += orbitprops.lx;
					prevorbitprops.ly += orbitprops.ly;
					prevorbitprops.lz += orbitprops.lz;
					prevorbitprops.E += orbitprops.E;
					prevorbitprops.ltot += orbitprops.ltot;
					prevorbitprops.mu += orbitprops.mu;
					prevorbitprops.hostlx += orbitprops.hostlx;
					prevorbitprops.hostly += orbitprops.hostlx;
					prevorbitprops.hostlz += orbitprops.hostlx;
					prevorbitprops.masslossrate += orbitprops.masslossrate;
					prevorbitprops.phi += orbitprops.phi;

					orbitprops = prevorbitprops;

					//Remove 0.5 from the total number of orbits so far
					orbitprops.numorbits-=0.5;
				}

				return;
			}

			//Find the average energy, total angular momentum and reduced mass
			ltot = orbitprops.ltot / (double)(currentsnap - orbitprops.prevpassagesnap);
			mu = orbitprops.mu / (double)(currentsnap - orbitprops.prevpassagesnap);

			tmporbitdata.orbitalenergy_ave = orbitprops.E / (double)(currentsnap - orbitprops.prevpassagesnap);;

			//The halos orbital eccentricity from the average properties
			tmporbitdata.orbitecc_wetzel2011 = sqrt(1.0 + ((2.0 * tmporbitdata.orbitalenergy_ave * ltot*ltot)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu)));

			//Find the orbits apo and pericentric distances
			tmporbitdata.rperi_wetzel2011 = (ltot*ltot)/((1+tmporbitdata.orbitecc_wetzel2011) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

			//Find the orbits apo and pericentric distances
			tmporbitdata.rapo_wetzel2011 = (ltot*ltot)/((1-tmporbitdata.orbitecc_wetzel2011) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu); //1.0 / ((2.0/r) - (vr*vr)/(Cosmo.G*hosthalo.mass));

			//The average mass loss rate
			tmporbitdata.masslossrate_ave = orbitprops.masslossrate / (double)(currentsnap - orbitprops.prevpassagesnap);

			//Lets output the average angular momentum since the last passage
			tmporbitdata.lxrel_ave = orbitprops.lx / (double)(currentsnap - orbitprops.prevpassagesnap);
			tmporbitdata.lyrel_ave = orbitprops.ly / (double)(currentsnap - orbitprops.prevpassagesnap);
			tmporbitdata.lzrel_ave = orbitprops.lz / (double)(currentsnap - orbitprops.prevpassagesnap);

			//Find the average angular momentum of the host halo
			orbitprops.hostlx /= (double)(currentsnap - orbitprops.prevpassagesnap);
			orbitprops.hostly /= (double)(currentsnap - orbitprops.prevpassagesnap);
			orbitprops.hostlz /= (double)(currentsnap - orbitprops.prevpassagesnap);

			//Now we have the average angular momentum,the alignment with the host angular momentum can be computed
			tmporbitdata.hostalignment = acos(((orbitprops.hostlx*tmporbitdata.lxrel_ave) + (orbitprops.hostly*tmporbitdata.lyrel_ave)	+ (orbitprops.hostlz*tmporbitdata.lzrel_ave))/
				(sqrt(orbitprops.hostlx*orbitprops.hostlx + orbitprops.hostly*orbitprops.hostly + orbitprops.hostlz*orbitprops.hostlz) * sqrt(tmporbitdata.lxrel_ave*tmporbitdata.lxrel_ave + tmporbitdata.lyrel_ave*tmporbitdata.lyrel_ave + tmporbitdata.lzrel_ave*tmporbitdata.lzrel_ave)));

			/* Compute all the angles */


			//If the reference angles hasn't been set then lets set it
			if(orbitprops.refangles == NULL){
				orbitprops.refangles = computeAngles(orbitprops.prevpassagepos,tmporbitdata);
			}
			else{
				//Compute the current angles
				currangles = computeAngles(orbitprops.prevpassagepos,tmporbitdata);

				tmporbitdata.longascnode = currangles[0] - orbitprops.refangles[0];
				tmporbitdata.inc = currangles[1] - orbitprops.refangles[1];
				tmporbitdata.argpariap = currangles[2] - orbitprops.refangles[2];
			}

			//The total angle moved through since last orbit
			tmporbitdata.phi = orbitprops.phi;

			prevorbitprops = orbitprops;

			//Reset the total angular momentum in the orbit props to zero
			orbitprops.lx = 0.0;
			orbitprops.ly = 0.0;
			orbitprops.lz = 0.0;
			orbitprops.E = 0.0;
			orbitprops.ltot = 0.0;
			orbitprops.mu = 0.0;
			orbitprops.hostlx = 0.0;
			orbitprops.hostly = 0.0;
			orbitprops.hostlz = 0.0;
			orbitprops.masslossrate = 0.0;
			orbitprops.phi = 0.0;
		}
		//Any additional properties to be calculated here

		//Store the index of the previous passage
		orbitprops.prevpassageindex = branchorbitdata.size();

		//Update the previous passage time
		orbitprops.prevpassagetime = tmporbitdata.uniage;

		//Mark the snapshot that this passage happens
		orbitprops.prevpassagesnap = currentsnap;

		//Store the previous entrytype;
		orbitprops.prevpassageentrytype = tmporbitdata.entrytype;

		//Store the radial vector for this passage
		orbitprops.prevpassagepos[0] = rx;
		orbitprops.prevpassagepos[1] = ry;
		orbitprops.prevpassagepos[2] = rz;

		//Now done the interpolation can check if this is the closest approach so far
		r = sqrt(tmporbitdata.xrel*tmporbitdata.xrel + tmporbitdata.yrel*tmporbitdata.yrel + tmporbitdata.zrel*tmporbitdata.zrel);
		if(r<orbitprops.closestapproach)
			orbitprops.closestapproach=r;
		tmporbitdata.closestapproach=orbitprops.closestapproach;

		//If not set to be orbiting then update the flag to say this object is orbiting
		if(orbitprops.orbitingflag==false)
			orbitprops.orbitingflag=true;

		//Now append it into the orbitdata dataset
		branchorbitdata.push_back(tmporbitdata);
	}

}

void AddFinalEntry(Options &opt,
	unsigned long long orbitID,
	int currentsnap, int prevsnap,
	HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo,
	vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata,
	vector<SnapData> &snapdata,
	OrbitProps &orbitprops,
	bool mergedflag){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*Cosmo.boxsize) orbitinghalo.x-=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)>0.5*Cosmo.boxsize) orbitinghalo.y-=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*Cosmo.boxsize) orbitinghalo.z-=Cosmo.boxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*Cosmo.boxsize) orbitinghalo.x+=Cosmo.boxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*Cosmo.boxsize) orbitinghalo.y+=Cosmo.boxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*Cosmo.boxsize) orbitinghalo.z+=Cosmo.boxsize;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r, vrel, vcirc, lx, ly, lz, vcomp, vradx, vrady, vradz, vtanx, vtany, vtanz, mu, deltat;

	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = hosthalo.x - orbitinghalo.x;
	ry = hosthalo.y - orbitinghalo.y;
	rz = hosthalo.z - orbitinghalo.z;
	vrx = hosthalo.vx - orbitinghalo.vx;
	vry = hosthalo.vy - orbitinghalo.vy;
	vrz = hosthalo.vz - orbitinghalo.vz;
	r = sqrt(rx * rx + ry * ry + rz * rz);
	vrel = sqrt(vrx*vrx + vry*vry + vrz*vrz);


	//Store the peak vmax in the orbiting halos history
	if(orbitinghalo.vmax>tmporbitdata.vmaxpeak)
		tmporbitdata.vmaxpeak = orbitinghalo.vmax;

	//Now done the interpolation can check if this is the closest approach so far
	if(r<orbitprops.closestapproach)
		tmporbitdata.closestapproach=r;

	//Store what orbitID number this is
	tmporbitdata.orbitID = orbitID;

	//Store the haloID this halo if it is not interpolated
	if(orbitinghalo.interpflag)
		tmporbitdata.orbithaloID = 0;
	else
		tmporbitdata.orbithaloID = orbitinghalo.id;


	//Keep track if this halo and its host is top of spatial herachy
	tmporbitdata.fieldhalo = orbitinghalo.fieldhalo;
	tmporbitdata.fieldhalohost = hosthalo.fieldhalo;

	//Set a flag if the orbit host merges with another halo
	tmporbitdata.hostmerges = orbitinghalo.hostmerges;

	//Store the information about it and its host substructure
	tmporbitdata.numsubstruct = orbitinghalo.numsubstruct;
	tmporbitdata.ratioofmassinsubsstruct = orbitinghalo.ratioofmassinsubsstruct;
	tmporbitdata.numsubstructhost = hosthalo.numsubstruct;
	tmporbitdata.ratioofmassinsubsstructhost = hosthalo.ratioofmassinsubsstruct;

	//Compute the instaneous mass loss rate

	//The difference in time since the previous snapshot
	deltat = snapdata[currentsnap].uniage - snapdata[prevsnap].uniage;

	// Find the change mass in units of Msun/Gyr
	tmporbitdata.masslossrate_inst = (prevorbitinghalo.mass - orbitinghalo.mass )/deltat;


	//Set this as a merger entry
	tmporbitdata.entrytype = 0.0;

	//The current number of orbits
	tmporbitdata.numorbits = orbitprops.numorbits;

	//Set the orbit period as -1.0 here as only calculated at the passages
	tmporbitdata.orbitperiod = -1.0;

	//The orbting halo
	tmporbitdata.haloID = orbitinghalo.origid;

	//The original Root Progenitor for the orbiting halo
	tmporbitdata.halorootprogenID = orbitinghalo.origrootprogenitor;

	//The host halo
	tmporbitdata.orbitedhaloID = hosthalo.origid;

	//The original Root Progenitor of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootprogenID = hosthalo.origrootprogenitor;

	//The age of the universe at this point
	tmporbitdata.uniage = snapdata[currentsnap].uniage;

	//The scalefactor of the universe
	tmporbitdata.scalefactor = snapdata[currentsnap].scalefactor;

	if((orbitprops.crossrvirtime>0) & (mergedflag))
		tmporbitdata.mergertimescale = snapdata[currentsnap].uniage - orbitprops.crossrvirtime;

	/* Calculate various properties to be outputted */

	//Find the components of the radial vector
	vcomp = (rx*vrx + ry*vry + rz*vrz)/(r*r);
	vradx = vcomp * rx;
	vrady = vcomp * ry;
	vradz = vcomp * rz;
	tmporbitdata.vrad = vcomp;

	//Then can use these components to find the components of the tangential velocity
	vtanx = vrx - vradx;
	vtany = vry - vrady;
	vtanz = vrz - vradz;
	tmporbitdata.vtan = sqrt(vtanx*vtanx + vtany*vtany + vtanz*vtanz);

	//The halos reduced mass
	mu = (orbitinghalo.mass * hosthalo.mass) / (orbitinghalo.mass + hosthalo.mass);

	//Angular momentum vectors
	lx = (ry * vrz) - (rz * vry);
	ly = -((rx * vrz) - (rz * vrx));
	lz = (rx * vry) - (ry * vrx);
	tmporbitdata.lxrel_inst += mu * lx;
	tmporbitdata.lyrel_inst += mu * ly;
	tmporbitdata.lzrel_inst += mu * lz;

	//Find the energy of the orbit
	tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * orbitinghalo.mass*hosthalo.mass)/r;

	//Find the radius of a circular orbit with the same energy
	tmporbitdata.rcirc = abs((Cosmo.G * orbitinghalo.mass * Menc(r,hosthalo.mass,hosthalo.cnfw,hosthalo.rvir))/(2.0 * tmporbitdata.orbitalenergy_inst));

	//Find the circular velocity of a circular orbit with the same energy
	vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

	//Can use these to find the orbital angular momentum of a circular orbit
	tmporbitdata.jcirc = mu * tmporbitdata.rcirc * vcirc;

	//Any additional properties to be calculated here

	//Set the orbit period and eccentricty as -1.0
	tmporbitdata.orbitperiod = -1.0;
	tmporbitdata.orbitecc_wetzel2011 = -1.0;
	tmporbitdata.lxrel_ave = -1.0;
	tmporbitdata.lyrel_ave = -1.0;
	tmporbitdata.lzrel_ave = -1.0;
	tmporbitdata.hostalignment = 0.0;


	//The halo properties
	tmporbitdata.npart = orbitinghalo.npart;
	tmporbitdata.mass = orbitinghalo.mass;
	tmporbitdata.rvir = orbitinghalo.rvir;
	tmporbitdata.vmax = orbitinghalo.vmax;
	tmporbitdata.rmax = orbitinghalo.rmax;
	tmporbitdata.cnfw = orbitinghalo.cnfw;

	//Extract the host halo properties
	tmporbitdata.nparthost = hosthalo.npart;
	tmporbitdata.masshost = hosthalo.mass;
	tmporbitdata.rvirhost = hosthalo.rvir;
	tmporbitdata.vmaxhost = hosthalo.vmax;
	tmporbitdata.rmaxhost = hosthalo.rmax;
	tmporbitdata.cnfwhost = hosthalo.cnfw;

	//Extract the posistions and velocities
	tmporbitdata.x = orbitinghalo.x;
	tmporbitdata.y = orbitinghalo.y;
	tmporbitdata.z = orbitinghalo.z;
	tmporbitdata.vx = orbitinghalo.vx;
	tmporbitdata.vy = orbitinghalo.vy;
	tmporbitdata.vz = orbitinghalo.vz;
	tmporbitdata.xrel = rx;
	tmporbitdata.yrel = ry;
	tmporbitdata.zrel = rz;
	tmporbitdata.vxrel = vrx;
	tmporbitdata.vyrel = vry;
	tmporbitdata.vzrel = vrz;

	//Now append it into the orbitdata dataset
	branchorbitdata.push_back(tmporbitdata);
}


void ProcessHalo(Options &opt, unsigned long long orbitID, int snap, unsigned long long index, 
	vector<SnapData> &snapdata, vector<OrbitData> &orbitdata,
	HaloData prevorbitinghalo, HaloData prevhosthalo,
	vector<OrbitData> branchorbitdata, OrbitData tmporbitdata,
	OrbitProps orbitprops, OrbitProps prevorbitprops,
	int *num_entrytypes){

	unsigned long long haloID = snapdata[snap].Halo[index].id;
	int halosnap = (int)(haloID/opt.TEMPORALHALOIDVAL);
	unsigned long long haloindex = (unsigned long long)(haloID%opt.TEMPORALHALOIDVAL-1);
	unsigned long long descendantID = snapdata[snap].Halo[index].descendant;
	int descendantsnap = (int)(descendantID/opt.TEMPORALHALOIDVAL);
	unsigned long long descendantindex = (unsigned long long)(descendantID%opt.TEMPORALHALOIDVAL-1);
	unsigned long long descendantProgenID = snapdata[descendantsnap].Halo[descendantindex].progenitor;


	//Keep track of this halo's snap and index for interpolation
	vector<int> halosnaps;
	vector<unsigned long long> haloindexes;
	vector<unsigned long long> hostindexes;
	vector<int> interpsnaps;

	//Store the index of the halo it is orbiting
	unsigned long long orbitinghaloindex;
	int prevsnap=halosnap;

	//Keep track of the snapshot
	int currentsnap = snap;

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

		//See if have reached the end of this branch or have merged with its host
		if((descendantID==haloID) | (descendantProgenID!=haloID)) break;

		//Lets move to the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;


		//Extract its descendant and its progenitor
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (int)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (unsigned long long)(descendantID%opt.TEMPORALHALOIDVAL-1);
		descendantProgenID = snapdata[descendantsnap].Halo[descendantindex].progenitor;

		//Interate keeping track of the snapshots
		currentsnap++;
	}

	//If a halo exist less than 10 snapshots then lets remove it from the catalogue
	if(((halosnaps.size()+interpsnaps.size())<10) | (halosnaps.size()<5)) return;

	//Lets setup the interpolation functions for the Position and Velocity of the orbiting halo
	SplineFuncs splinefuncs(halosnaps.size());
	SetupPosVelInterpFunctions(halosnaps,haloindexes,snapdata,splinefuncs);


	//If the interp snapshots contains snapshots then interpolation needs to be done
	if(interpsnaps.size()>0) InterpHaloProps(opt,halosnaps,haloindexes,interpsnaps,snapdata,splinefuncs);

	for(int i=0;i<halosnaps.size();i++){
		hostindexes.push_back((unsigned long long)(snapdata[halosnaps[i]].Halo[haloindexes[i]].orbitedhaloid%opt.TEMPORALHALOIDVAL-1));
	}

	//Now the orbiting halo has been interpolated the interpolation functions can be setup for the host halo
	SplineFuncs hostsplinefuncs(halosnaps.size());
	SetupPosVelInterpFunctionsHost(halosnaps,hostindexes,haloindexes,snapdata,hostsplinefuncs);


	//Reset the tree info to back at the base of the tree
	haloID = snapdata[snap].Halo[index].id;
	halosnap = (int)(haloID/opt.TEMPORALHALOIDVAL);
	haloindex = (unsigned long long)(haloID%opt.TEMPORALHALOIDVAL-1);
	descendantID = snapdata[snap].Halo[index].descendant;
	descendantsnap = (int)(descendantID/opt.TEMPORALHALOIDVAL);
	descendantindex = (unsigned long long)(descendantID%opt.TEMPORALHALOIDVAL-1);
	descendantProgenID = snapdata[descendantsnap].Halo[descendantindex].progenitor;


	//Flag to keep track if the halo has merged
	bool mergedflag = false;


	// ofstream file;
	// file.open("../analysis/data/circ.dat");

	// ofstream file2;
	// file2.open("../analysis/data/data2.dat");

	while(true){

		//Reset all the data to zero
		tmporbitdata={};

		//Extract the halo it is orbiting at this snapshot
		orbitinghaloindex = (unsigned long long)(snapdata[halosnap].Halo[haloindex].orbitedhaloid%opt.TEMPORALHALOIDVAL-1);


		//Lets set this halos orbit data
		if(halosnap!=prevsnap)
			CalcOrbitProps(opt,
				orbitID,
				halosnap, prevsnap,
				snapdata[halosnap].Halo[haloindex], snapdata[halosnap].Halo[orbitinghaloindex], prevorbitinghalo,prevhosthalo,
				branchorbitdata,tmporbitdata,
				snapdata,
				num_entrytypes, orbitprops, prevorbitprops,
				splinefuncs, hostsplinefuncs);

		// if(find(interpsnaps.begin(), interpsnaps.end(), halosnap) != interpsnaps.end())
		// 	file<<-halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;
		// else
		// 	file<<halosnap<<" "<<snapdata[halosnap].Halo[haloindex].x - snapdata[halosnap].Halo[orbitinghaloindex].x<<" "<<snapdata[halosnap].Halo[haloindex].y - snapdata[halosnap].Halo[orbitinghaloindex].y<<" "<<snapdata[halosnap].Halo[haloindex].z - snapdata[halosnap].Halo[orbitinghaloindex].z<<endl;

		//Mark this halo as being done:
		snapdata[halosnap].Halo[haloindex].doneflag = true;

		// cout<<haloID<<" "<<descendantID<<" "<<descendantProgenID<<" "<<snapdata[halosnap].Halo[orbitinghaloindex].id<<" "<<snapdata[halosnap].Halo[orbitinghaloindex].descendant<<endl;
		//See if have reached the end of this branch or has merged with its host
		if(descendantID==haloID)
			break;
		else if(descendantProgenID!=haloID){
			mergedflag=true;
			break;
		}

		//Update the previous halo and its host
		prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
		prevhosthalo = snapdata[halosnap].Halo[orbitinghaloindex];
		prevsnap = halosnap;

		//lets move onto the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;

		//Extract its descendant and its progenitor
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (int)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (unsigned long long)(descendantID%opt.TEMPORALHALOIDVAL-1);
		descendantProgenID = snapdata[descendantsnap].Halo[descendantindex].progenitor;

	}

	//If there are no data to be outputted the can continue onto the next halo
	if(branchorbitdata.size()==0) return;
	// file2.close();
	// file.close();

	//Output a final entry point if the halo has undergone a crossing or passage point
	AddFinalEntry(opt,
				orbitID,
				halosnap, prevsnap,
				snapdata[halosnap].Halo[haloindex], snapdata[halosnap].Halo[orbitinghaloindex], prevorbitinghalo,
				branchorbitdata,tmporbitdata,
				snapdata,
				orbitprops,
				mergedflag);

	//Set the total number of orbits for all entries
	for(int i = 0; i<branchorbitdata.size(); i++){
		branchorbitdata[i].totnumorbits = orbitprops.numorbits;
		branchorbitdata[i].minrscale = orbitprops.minrscale;
		branchorbitdata[i].minrmax = orbitprops.minrmax;
	}

	double simtime = snapdata.back().uniage - snapdata.front().uniage;

	//If merged then set the MergedFlag == tur
	if(mergedflag) for(int i = 0; i<branchorbitdata.size(); i++) branchorbitdata[i].mergedflag=true;

	// CleanOrbits(branchorbitdata,simtime);

	//Now finished with this branches orbital calculations so it can be added
	//into the orbitdata vector that contains all halos, it only needs to be
	//moved rather than copied
	orbitdata.insert(orbitdata.end(),make_move_iterator(branchorbitdata.begin()),make_move_iterator(branchorbitdata.end()));
}

void ProcessOrbits(Options &opt, vector<SnapData> &snapdata, vector<OrbitData> &orbitdata){

	unsigned long long orbitID = 0;
	int *num_entrytypes = new int[opt.totnumtypeofentries];

	//Keep track of the previous halos halodata and orbitdata
	HaloData prevorbitinghalo;
	HaloData prevhosthalo;
	vector<OrbitData> branchorbitdata;
	OrbitData tmporbitdata;

	//Keep track of the properties of each orbit
	OrbitProps orbitprops;
	OrbitProps prevorbitprops;

	// Now lets start at the starting snapshot and walk up the tree
	// calculating the orbit relative to the halo which it was found
	// to be orbiting
	// int snap = 55;
	// int snap = 116;
	for(int snap=0;snap<opt.numsnaps;snap++){
	// unsigned long long i = 990;
	// unsigned long long i = 1177;
		for(unsigned long long i=0;i<snapdata[snap].numhalos;i++){

			// Lets first check if this halo has been processed or is not orbiting a halo
			if((snapdata[snap].Halo[i].doneflag) | (snapdata[snap].Halo[i].orbitedhaloid==-1)) continue;

			//Reset all data to zero
			for(int j=0; j<opt.totnumtypeofentries; j++) num_entrytypes[j]=0;
			prevorbitinghalo = {};
			prevhosthalo = {};
			branchorbitdata.clear();
			tmporbitdata = {};
			orbitprops = {};
			prevorbitprops = {};

			ProcessHalo(opt,orbitID,snap,i,
				snapdata,orbitdata,
				prevorbitinghalo, prevhosthalo,
				branchorbitdata, tmporbitdata,
				orbitprops,prevorbitprops,
				num_entrytypes);

			orbitID++;

		}
		if(opt.iverbose) cout<<"Done processing snap "<<snap<<endl;
	}
}
