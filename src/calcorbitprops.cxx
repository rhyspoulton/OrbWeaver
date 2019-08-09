

#include "orbweaver.h"

void CalcOrbitProps(Options &opt,
	unsigned long long orbitID,
	int currentsnap, int prevsnap,
	HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo,
	vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata,
	vector<SnapData> &snapdata,
	vector<int> num_entrytypes, OrbitProps &orbitprops, vector<OrbitProps> &passagesorbitprops,
	SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs){

	//First correct for periodicity compared to the host halo
	if((orbitinghalo.x - hosthalo.x)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.x-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.y - hosthalo.y)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.y-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.z-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.x+=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.y+=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.z+=snapdata[currentsnap].physboxsize;

	if((prevorbitinghalo.x - prevhosthalo.x)>0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.x-=snapdata[currentsnap].physboxsize;
	if((prevorbitinghalo.y - prevhosthalo.y)>0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.y-=snapdata[currentsnap].physboxsize;
	if((prevorbitinghalo.z - prevhosthalo.z)>0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.z-=snapdata[currentsnap].physboxsize;
	if((prevorbitinghalo.x - prevhosthalo.x)<-0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.x+=snapdata[currentsnap].physboxsize;
	if((prevorbitinghalo.y - prevhosthalo.y)<-0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.y+=snapdata[currentsnap].physboxsize;
	if((prevorbitinghalo.z - prevhosthalo.z)<-0.5*snapdata[currentsnap].physboxsize) prevorbitinghalo.z+=snapdata[currentsnap].physboxsize;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r,rcomove,vrad,vrel;

	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = orbitinghalo.x - hosthalo.x;
	ry = orbitinghalo.y - hosthalo.y;
	rz = orbitinghalo.z - hosthalo.z;
	r = sqrt(rx * rx + ry * ry + rz * rz);

	// Calculate the velocity relative to the host including the hubble flow
	vrx = orbitinghalo.vx - hosthalo.vx + r * snapdata[currentsnap].Hz;
	vry = orbitinghalo.vy - hosthalo.vy + r * snapdata[currentsnap].Hz;
	vrz = orbitinghalo.vz - hosthalo.vz + r * snapdata[currentsnap].Hz;
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

	//Can see if this is the closest approach, this needs to be compared in comoving
	rcomove = r * Cosmo.h / snapdata[currentsnap].scalefactor;
	if(rcomove<orbitprops.closestapproach)
		orbitprops.closestapproach = rcomove;
		orbitprops.closestapproachscalefactor = snapdata[currentsnap].scalefactor;

	//Put into the output data and convert back to physical
	tmporbitdata.closestapproach = orbitprops.closestapproach * orbitprops.closestapproachscalefactor/ Cosmo.h;
	tmporbitdata.closestapproachscalefactor = orbitprops.closestapproachscalefactor;

	//Store what orbitID number this is
	tmporbitdata.orbitID = orbitID;

	//The orbting halo
	tmporbitdata.haloID = orbitinghalo.origid;

	//The original Root Progenitor for the orbiting halo
	tmporbitdata.halorootprogenID = orbitinghalo.origrootprogenitor;

	//The original Root Descendant for the orbiting halo
	tmporbitdata.halorootdescenID = orbitinghalo.origrootdescendant;

	//The host halo
	tmporbitdata.orbitedhaloID = hosthalo.origid;

	//The original Root Progenitor of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootprogenID = hosthalo.origrootprogenitor;

	//The original Root Descendant of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootdescenID = hosthalo.origrootdescendant;

	//Store the haloID this halo if it is not interpolated
	if(orbitinghalo.interpflag)
		tmporbitdata.orbithaloID = 0;
	else
		tmporbitdata.orbithaloID = orbitinghalo.id;

	double prevrx,prevry,prevrz,prevvrx,prevvry,prevvrz,prevr,prevvrad;

	//Lets find the previous orbiting halo host distance
	prevrx = prevorbitinghalo.x - prevhosthalo.x;
	prevry = prevorbitinghalo.y - prevhosthalo.y;
	prevrz = prevorbitinghalo.z - prevhosthalo.z;
	prevr = sqrt(prevrx * prevrx + prevry * prevry + prevrz * prevrz);

	// Calculate the velocity relative to the host including the hubble flow
	prevvrx = prevorbitinghalo.vx - prevhosthalo.vx + prevr * snapdata[prevsnap].Hz;
	prevvry = prevorbitinghalo.vy - prevhosthalo.vy + prevr * snapdata[prevsnap].Hz;
	prevvrz = prevorbitinghalo.vz - prevhosthalo.vz + prevr * snapdata[prevsnap].Hz;
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
	double ltot, f, vcirc, jcirc, vcomp, vradx, vrady, vradz, vtanx, vtany, vtanz, prevpassager, hostlx, hostly, hostlz;
	vector<double> currangles;
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
		vrx = tmporbitdata.vxrel + r * GetH(tmporbitdata.scalefactor);
		vry = tmporbitdata.vyrel + r * GetH(tmporbitdata.scalefactor);
		vrz = tmporbitdata.vzrel + r * GetH(tmporbitdata.scalefactor);
		vrel = sqrt(vrx*vrx + vry*vry + vrz*vrz);

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
		mu = (tmporbitdata.mass * tmporbitdata.masshost) / (tmporbitdata.mass + tmporbitdata.masshost);

		//Angular momentum vectors
		lx = (ry * vrz) - (rz * vry);
		ly = -((rx * vrz) - (rz * vrx));
		lz = (rx * vry) - (ry * vrx);
		ltot = mu* sqrt(lx*lx + ly*ly + lz*lz);
		tmporbitdata.lxrel_inst += mu * lx;
		tmporbitdata.lyrel_inst += mu * ly;
		tmporbitdata.lzrel_inst += mu * lz;

		//Find the energy of the orbit
		tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * tmporbitdata.mass*tmporbitdata.masshost)/r;

		// only compute Rcirc and Jcirc if on a bound orbit (E<0)
		if(tmporbitdata.orbitalenergy_inst<0){

			//Find the radius of a circular orbit with the same energy
			tmporbitdata.rcirc = abs((Cosmo.G * tmporbitdata.mass * Menc(r,tmporbitdata.masshost,tmporbitdata.cnfwhost,tmporbitdata.rvirhost))/(2.0 * tmporbitdata.orbitalenergy_inst));


			//Find the circular velocity of a circular orbit with the same energy
			vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

			//Can use these to find the orbital angular momentum of a circular orbit
			jcirc = mu * tmporbitdata.rcirc * vcirc;

			//Compute the circularity for the orbit (eta)
			tmporbitdata.eta = ltot/jcirc;
		}

		//Any additional properties to be calculated here

		//The halos orbital eccentricity calculated from energy and angular momentum
		tmporbitdata.orbitecc_calc = sqrt(1.0 + ((2.0 * tmporbitdata.orbitalenergy_inst * ltot*ltot)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu)));

		//Find the orbits apo and pericentric distances
		tmporbitdata.rperi_calc = (ltot*ltot)/((1+tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

		//Find the orbits apo and pericentric distances
		tmporbitdata.rapo_calc = (ltot*ltot)/((1-tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

		//Set the orbit period and eccentricty as -1.0
		tmporbitdata.orbitperiod = -1.0;
		tmporbitdata.lxrel_ave = -1.0;
		tmporbitdata.lyrel_ave = -1.0;
		tmporbitdata.lzrel_ave = -1.0;
		tmporbitdata.hostalignment = 0.0;

		//Store the previous crossing point entry type
		orbitprops.prevcrossingentrytype=numrvircrossing;

		//Store the snapshot which this crossing point happened
		orbitprops.prevcrossingsnap = currentsnap;

		//Now done the interpolation can check if this is the closest approach so far which needs to be done in comoving
		rcomove = r * Cosmo.h / tmporbitdata.scalefactor;
		if(rcomove<orbitprops.closestapproach)
			orbitprops.closestapproach = rcomove;
			orbitprops.closestapproachscalefactor = tmporbitdata.scalefactor;

		//Put into the output data and convert back to physical
		tmporbitdata.closestapproach = orbitprops.closestapproach * orbitprops.closestapproachscalefactor/ Cosmo.h;
		tmporbitdata.closestapproachscalefactor = orbitprops.closestapproachscalefactor;

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


			//Remove 1 from the num entrytype
			tmporbitdata.num_entrytype--;
			if(vrad>0){ //Peri-centric
				num_entrytypes[0]--;
			}
			else{    //Apocentric
				num_entrytypes[1]--;
			}

			if(orbitprops.numorbits<=0.5){
				orbitprops.orbitingflag=false;
			}
			else{

				//Add the orbit average qantities to the previous passage orbit props
				passagesorbitprops.back().lx += orbitprops.lx;
				passagesorbitprops.back().ly += orbitprops.ly;
				passagesorbitprops.back().lz += orbitprops.lz;
				passagesorbitprops.back().E += orbitprops.E;
				passagesorbitprops.back().mu += orbitprops.mu;
				passagesorbitprops.back().hostlx += orbitprops.hostlx;
				passagesorbitprops.back().hostly += orbitprops.hostlx;
				passagesorbitprops.back().hostlz += orbitprops.hostlx;
				passagesorbitprops.back().masslossrate += orbitprops.masslossrate;
				passagesorbitprops.back().phi += orbitprops.phi;

				//Keep the closest approach, the scalefactor it happened at and the minmum Rscale and Rmax
				passagesorbitprops.back().closestapproach = orbitprops.closestapproach;
				passagesorbitprops.back().closestapproachscalefactor = orbitprops.closestapproachscalefactor;
				passagesorbitprops.back().minrmax = orbitprops.minrmax;
				passagesorbitprops.back().minrscale = orbitprops.minrscale;

				orbitprops = passagesorbitprops.back();

				//Delete the last entry in the passagesorbitprops
				passagesorbitprops.pop_back();
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
		vrx = tmporbitdata.vxrel + r * GetH(tmporbitdata.scalefactor);
		vry = tmporbitdata.vyrel + r * GetH(tmporbitdata.scalefactor);
		vrz = tmporbitdata.vzrel + r * GetH(tmporbitdata.scalefactor);

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
		mu = (tmporbitdata.mass * tmporbitdata.masshost) / (tmporbitdata.mass + tmporbitdata.masshost);

		//Angular momentum vectors
		lx = (ry * vrz) - (rz * vry);
		ly = -((rx * vrz) - (rz * vrx));
		lz = (rx * vry) - (ry * vrx);
		ltot = mu* sqrt(lx*lx + ly*ly + lz*lz);
		tmporbitdata.lxrel_inst += mu * lx;
		tmporbitdata.lyrel_inst += mu * ly;
		tmporbitdata.lzrel_inst += mu * lz;

		//Find the energy of the orbit
		tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * tmporbitdata.mass*tmporbitdata.masshost)/r;

		// only compute Rcirc and Jcirc if on a bound orbit (E<0)
		if(tmporbitdata.orbitalenergy_inst<0){

			//Find the radius of a circular orbit with the same energy
			tmporbitdata.rcirc = abs((Cosmo.G * tmporbitdata.mass * Menc(r,tmporbitdata.masshost,tmporbitdata.cnfwhost,tmporbitdata.rvirhost))/(2.0 * tmporbitdata.orbitalenergy_inst));


			//Find the circular velocity of a circular orbit with the same energy
			vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

			//Can use these to find the orbital angular momentum of a circular orbit
			jcirc = mu * tmporbitdata.rcirc * vcirc;

			//Compute the circularity for the orbit (eta)
			tmporbitdata.eta = ltot/jcirc;
		}

		//The halos orbital eccentricity calculated from energy and angular momentum
		tmporbitdata.orbitecc_calc = sqrt(1.0 + ((2.0 * tmporbitdata.orbitalenergy_inst * ltot*ltot)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu)));

		//Find the orbits apo and pericentric distances
		tmporbitdata.rperi_calc = (ltot*ltot)/((1+tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

		//Find the orbits apo and pericentric distances
		tmporbitdata.rapo_calc = (ltot*ltot)/((1-tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);


		//Now done the interpolation can check if this is the closest approach so far which needs to be done in comoving
		rcomove = r * Cosmo.h / tmporbitdata.scalefactor;
		if(rcomove<orbitprops.closestapproach)
			orbitprops.closestapproach = rcomove;
			orbitprops.closestapproachscalefactor = tmporbitdata.scalefactor;

		//Put into the output data and convert back to physical
		tmporbitdata.closestapproach = orbitprops.closestapproach * orbitprops.closestapproachscalefactor/ Cosmo.h;
		tmporbitdata.closestapproachscalefactor = orbitprops.closestapproachscalefactor;

		/* Calculate various properties to be outputted if the halo is marked as orbiting */

		if(orbitprops.orbitingflag){
			// Calculate the orbit eccentricity
			if(vrad>0)
				tmporbitdata.orbitecc_ratio = (orbitprops.prevpassagercomove-rcomove)/(orbitprops.prevpassagercomove+rcomove);
			else
				tmporbitdata.orbitecc_ratio = (rcomove-orbitprops.prevpassagercomove)/(rcomove+orbitprops.prevpassagercomove);

			// //Calculate the orbit period as 2x the previous passage
			tmporbitdata.orbitperiod = 2.0* (tmporbitdata.uniage - orbitprops.prevpassagetime);

			// semiMajor = (prevpassager+r)/2.0; //sqrt((rx-orbitprops.prevpassagepos[0])*(rx-orbitprops.prevpassagepos[0]) + (ry-orbitprops.prevpassagepos[1])*(ry-orbitprops.prevpassagepos[1]) + (rz-orbitprops.prevpassagepos[2])*(rz-orbitprops.prevpassagepos[2]));

			//Remove any passages that happen within 2 snapshots since these are not sampled correctly
			if((currentsnap - orbitprops.prevpassagesnap) < 3){

				//Add the previous passage to the remove indexes
				branchorbitdata.erase(branchorbitdata.begin()+orbitprops.prevpassageindex);

				//Reset the num_entrytype for this passage and the previous passage
				num_entrytypes[0]--;
				num_entrytypes[1]--;

				// If the object has only undergone one orbit then lets remove the passages and reset so the halo is 
				// no longer set to be orbiting, otherwise update the previous passage quantities
				if(tmporbitdata.numorbits==1.0){
					orbitprops.orbitingflag=false;

					//Reset the total angular momentum in the orbit props to zero
					orbitprops.lx = 0.0;
					orbitprops.ly = 0.0;
					orbitprops.lz = 0.0;
					orbitprops.E = 0.0;
					orbitprops.mu = 0.0;
					orbitprops.hostlx = 0.0;
					orbitprops.hostly = 0.0;
					orbitprops.hostlz = 0.0;
					orbitprops.masslossrate = 0.0;
					orbitprops.phi=0.0;
					orbitprops.numorbits=0.0;
				}
				else{

					//Add the orbit average qantities to the previous passage orbit props
					passagesorbitprops.back().lx += orbitprops.lx;
					passagesorbitprops.back().ly += orbitprops.ly;
					passagesorbitprops.back().lz += orbitprops.lz;
					passagesorbitprops.back().E += orbitprops.E;
					passagesorbitprops.back().mu += orbitprops.mu;
					passagesorbitprops.back().hostlx += orbitprops.hostlx;
					passagesorbitprops.back().hostly += orbitprops.hostlx;
					passagesorbitprops.back().hostlz += orbitprops.hostlx;
					passagesorbitprops.back().masslossrate += orbitprops.masslossrate;
					passagesorbitprops.back().phi += orbitprops.phi;

					//Keep the closest approach, the scalefactor it happened at and the minmum Rscale and Rmax
					passagesorbitprops.back().closestapproach = orbitprops.closestapproach;
					passagesorbitprops.back().closestapproachscalefactor = orbitprops.closestapproachscalefactor;
					passagesorbitprops.back().minrmax = orbitprops.minrmax;
					passagesorbitprops.back().minrscale = orbitprops.minrscale;


					orbitprops = passagesorbitprops.back();

					//Delete the last entry in the passagesorbitprops
					passagesorbitprops.pop_back();

					//Remove 0.5 from the total number of orbits so far
					orbitprops.numorbits-=0.5;
				}

				return;
			}

			//Find the average energy and reduced mass
			tmporbitdata.orbitalenergy_ave = orbitprops.E / (double)(currentsnap - orbitprops.prevpassagesnap);;

			//The average mass loss rate
			tmporbitdata.masslossrate_ave = orbitprops.masslossrate / (double)(currentsnap - orbitprops.prevpassagesnap);

			//Lets output the average angular momentum since the last passage
			tmporbitdata.lxrel_ave = orbitprops.lx / (double)(currentsnap - orbitprops.prevpassagesnap);
			tmporbitdata.lyrel_ave = orbitprops.ly / (double)(currentsnap - orbitprops.prevpassagesnap);
			tmporbitdata.lzrel_ave = orbitprops.lz / (double)(currentsnap - orbitprops.prevpassagesnap);

			//Find the average angular momentum of the host halo
			hostlx = orbitprops.hostlx / (double)(currentsnap - orbitprops.prevpassagesnap);
			hostly = orbitprops.hostly / (double)(currentsnap - orbitprops.prevpassagesnap);
			hostlz = orbitprops.hostlz / (double)(currentsnap - orbitprops.prevpassagesnap);

			//Now we have the average angular momentum,the alignment with the host angular momentum can be computed
			tmporbitdata.hostalignment = acos(((hostlx*tmporbitdata.lxrel_ave) + (hostly*tmporbitdata.lyrel_ave)	+ (hostlz*tmporbitdata.lzrel_ave))/
				(sqrt(hostlx*hostlx + hostly*hostly + hostlz*hostlz) * sqrt(tmporbitdata.lxrel_ave*tmporbitdata.lxrel_ave + tmporbitdata.lyrel_ave*tmporbitdata.lyrel_ave + tmporbitdata.lzrel_ave*tmporbitdata.lzrel_ave)));

			/* Compute all the angles */


			//If the reference angles hasn't been set then lets set it
			if(accumulate(orbitprops.refangles.begin(),orbitprops.refangles.end(),0) == 0){
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

			//Store the index and snapshot of this passgae
			orbitprops.passageindex = branchorbitdata.size();
			orbitprops.passagesnap = currentsnap;

			//Add the current orbit props to the passage orbit props vector
			passagesorbitprops.push_back(orbitprops);

			//Reset the total angular momentum in the orbit props to zero
			orbitprops.lx = 0.0;
			orbitprops.ly = 0.0;
			orbitprops.lz = 0.0;
			orbitprops.E = 0.0;
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

		//Store the radial vector for this passage and the comoving radius
		orbitprops.prevpassagepos[0] = rx;
		orbitprops.prevpassagepos[1] = ry;
		orbitprops.prevpassagepos[2] = rz;
		orbitprops.prevpassagercomove = rcomove;

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
	if((orbitinghalo.x - hosthalo.x)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.x-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.y - hosthalo.y)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.y-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.z - hosthalo.z)>0.5*snapdata[currentsnap].physboxsize) orbitinghalo.z-=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.x - hosthalo.x)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.x+=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.y - hosthalo.y)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.y+=snapdata[currentsnap].physboxsize;
	if((orbitinghalo.z - hosthalo.z)<-0.5*snapdata[currentsnap].physboxsize) orbitinghalo.z+=snapdata[currentsnap].physboxsize;

	//This is where all the orbital properties are calculate for the halo at this snapshot
	double rx,ry,rz,vrx,vry,vrz,r,rcomove, vrel, vcirc, jcirc, ltot, lx, ly, lz, vcomp, vradx, vrady, vradz, vtanx, vtany, vtanz, mu, deltat;


	// Find the orbitinghalos distance to the hosthalo and its orbiting vector
	rx = orbitinghalo.x - hosthalo.x;
	ry = orbitinghalo.y - hosthalo.y;
	rz = orbitinghalo.z - hosthalo.z;
	r = sqrt(rx * rx + ry * ry + rz * rz);

	// Calculate the velocity relative to the host including the hubble flow
	vrx = orbitinghalo.vx - hosthalo.vx + r * snapdata[currentsnap].Hz;
	vry = orbitinghalo.vy - hosthalo.vy + r * snapdata[currentsnap].Hz;
	vrz = orbitinghalo.vz - hosthalo.vz + r * snapdata[currentsnap].Hz;
	vrel = sqrt(vrx*vrx + vry*vry + vrz*vrz);


	//Store the peak vmax in the orbiting halos history
	if(orbitinghalo.vmax>tmporbitdata.vmaxpeak)
		tmporbitdata.vmaxpeak = orbitinghalo.vmax;

	//Can see if this is the closest approach, this needs to be compared in comoving
	rcomove = r * Cosmo.h / snapdata[currentsnap].scalefactor;
	if(rcomove<orbitprops.closestapproach)
		orbitprops.closestapproach = rcomove;
		orbitprops.closestapproachscalefactor = snapdata[currentsnap].scalefactor;

	//Put into the output data and convert back to physical
	tmporbitdata.closestapproach = orbitprops.closestapproach * orbitprops.closestapproachscalefactor/ Cosmo.h;
	tmporbitdata.closestapproachscalefactor = orbitprops.closestapproachscalefactor;

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

	//The original Root Descendant for the orbiting halo
	tmporbitdata.halorootdescenID = orbitinghalo.origrootdescendant;

	//The host halo
	tmporbitdata.orbitedhaloID = hosthalo.origid;

	//The original Root Progenitor of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootprogenID = hosthalo.origrootprogenitor;

	//The original Root Descendant of the orbited halo from the halo catalog
	tmporbitdata.orbitedhaloorigrootdescenID = hosthalo.origrootdescendant;

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
	ltot = mu*sqrt(lx*lx + ly*ly + lz*lz);
	tmporbitdata.lxrel_inst += mu * lx;
	tmporbitdata.lyrel_inst += mu * ly;
	tmporbitdata.lzrel_inst += mu * lz;

	//Find the energy of the orbit
	tmporbitdata.orbitalenergy_inst = 0.5 * mu  * (vrel*vrel) - (Cosmo.G * orbitinghalo.mass*hosthalo.mass)/r;

	// only compute Rcirc and Jcirc if on a bound orbit (E<0)
	if(tmporbitdata.orbitalenergy_inst<0){

		//Find the radius of a circular orbit with the same energy
		tmporbitdata.rcirc = abs((Cosmo.G * orbitinghalo.mass * Menc(r,hosthalo.mass,hosthalo.cnfw,hosthalo.rvir))/(2.0 * tmporbitdata.orbitalenergy_inst));


		//Find the circular velocity of a circular orbit with the same energy
		vcirc = sqrt((- 2.0 *tmporbitdata.orbitalenergy_inst)/mu);

		//Can use these to find the orbital angular momentum of a circular orbit
		jcirc = mu * tmporbitdata.rcirc * vcirc;

		//Compute the circularity for the orbit (eta)
		tmporbitdata.eta = ltot/jcirc;
	}

	//Any additional properties to be calculated here

	//The halos orbital eccentricity calculated from energy and angular momentum
	tmporbitdata.orbitecc_calc = sqrt(1.0 + ((2.0 * tmporbitdata.orbitalenergy_inst * ltot*ltot)/((Cosmo.G * orbitinghalo.mass * hosthalo.mass)*(Cosmo.G * orbitinghalo.mass * hosthalo.mass) * mu)));

	//Find the orbits apo and pericentric distances
	tmporbitdata.rperi_calc = (ltot*ltot)/((1+tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

	//Find the orbits apo and pericentric distances
	tmporbitdata.rapo_calc = (ltot*ltot)/((1-tmporbitdata.orbitecc_calc) * Cosmo.G * orbitinghalo.mass * hosthalo.mass  * mu);

	//Set the orbit period and eccentricty as -1.0
	tmporbitdata.orbitperiod = -1.0;
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
	OrbitProps orbitprops, vector<OrbitProps> passagesorbitprops,
	vector <int> num_entrytypes){

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

	//If a halo exist less than 3 snapshots in the orbit catalogue then lets remove it
	if(halosnaps.size()<3) return;

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
				num_entrytypes, orbitprops, passagesorbitprops,
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

	//Reset all the data to zero
	tmporbitdata={};

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

	//If merged then set the MergedFlag == tur
	if(mergedflag) for(int i = 0; i<branchorbitdata.size(); i++) branchorbitdata[i].mergedflag=true;

	//Clean the orbits of "wobbles"
	CleanOrbits(opt,branchorbitdata,passagesorbitprops);

	if(branchorbitdata.size()==1) return;

	//Now finished with this branches orbital calculations so it can be added
	//into the orbitdata vector that contains all halos, it only needs to be
	//moved rather than copied
	orbitdata.insert(orbitdata.end(),make_move_iterator(branchorbitdata.begin()),make_move_iterator(branchorbitdata.end()));
}

void ProcessOrbits(Options &opt, vector<SnapData> &snapdata, vector<OrbitData> &orbitdata){

	unsigned long long orbitID = 0;
	vector<int> num_entrytypes (opt.totnumtypeofentries);

	//Keep track of the previous halos halodata and orbitdata
	HaloData prevorbitinghalo;
	HaloData prevhosthalo;
	vector<OrbitData> branchorbitdata;
	OrbitData tmporbitdata;

	//Keep track of the properties of each orbit
	OrbitProps orbitprops;
	vector<OrbitProps> passagesorbitprops;

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
			fill(num_entrytypes.begin(),num_entrytypes.end(),0);
			prevorbitinghalo = {};
			prevhosthalo = {};
			branchorbitdata.clear();
			tmporbitdata = {};
			orbitprops = {};
			passagesorbitprops.clear();

			ProcessHalo(opt,orbitID,snap,i,
				snapdata,orbitdata,
				prevorbitinghalo, prevhosthalo,
				branchorbitdata, tmporbitdata,
				orbitprops,passagesorbitprops,
				num_entrytypes);

			orbitID++;

		}
		if(opt.iverbose) cout<<"Done processing snap "<<snap<<endl;
	}
}
