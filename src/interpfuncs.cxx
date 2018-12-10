
#include "orbweaver.h"

double LogInterp(double prevdata, double nextdata, double f){
	return pow(nextdata,f) * pow(prevdata,1-f);
}
double LinInterp(double prevdata, double nextdata, double f){
	return prevdata + (nextdata - prevdata)*f;
}

void SetupPosVelInterpFunctions(vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<SnapData> &snapdata, SplineFuncs &splinefuncs){

	Int_t nhalo = halosnaps.size();
	double boundry=0, nextpos;
	double halouniages[nhalo], x[nhalo], y[nhalo], z[nhalo], vx[nhalo], vy[nhalo], vz[nhalo];


	//Lets extract the uniage for the interpolation routine for the snapshots that the halo is present and keep track
	for(int i = 0;i<nhalo;i++)
		halouniages[i] =  snapdata[halosnaps[i]].uniage;

	/* xpos */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		x[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].x+boundry;
		if(i<nhalo-1){
			nextpos=snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].x+boundry;
			if(nextpos-x[i]>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(nextpos-x[i]<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.x, halouniages, x, nhalo);


	/*  y-pos  */

	//Reset the boundary
	boundry=0;

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		y[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].y+boundry;
		if(i<nhalo-1){
			nextpos=snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].y+boundry;
			if(nextpos-y[i]>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(nextpos-y[i]<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.y, halouniages, y, nhalo);


	/*  z-pos  */

	//Reset the boundary
	boundry=0;

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		z[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].z+boundry;
		if(i<nhalo-1){
			nextpos=snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].z+boundry;
			if(nextpos-z[i]>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(nextpos-z[i]<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.z, halouniages, z, nhalo);

	/*  x-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vx[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vx;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vx, halouniages, vx, nhalo);

	/*  y-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vy[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vy;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vy, halouniages, vy, nhalo);

	/*  z-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vz[i] = snapdata[halosnaps[i]].Halo[haloindexes[i]].vz;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vz, halouniages, vz, nhalo);
}


void SetupPosVelInterpFunctionsHost(vector<Int_t> &halosnaps, vector<Int_t> &hostindexes, vector<Int_t> &haloindexes, vector<SnapData> &snapdata, SplineFuncs &splinefuncs){


	//Need to interpolate the host's position and do periodicity corrections based on the oribiting halo
	Int_t nhalo = halosnaps.size();
	double boundry=0, nextpos;
	double halouniages[nhalo], x[nhalo], y[nhalo], z[nhalo], vx[nhalo], vy[nhalo], vz[nhalo];
	double tmppos, orbitinghalopos, orbitinghalonextpos;


	//Lets extract the uniage for the interpolation routine for the snapshots that the halo is present and keep track
	for(int i = 0;i<nhalo;i++)
		halouniages[i] =  snapdata[halosnaps[i]].uniage;

	/* xpos */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		tmppos = snapdata[halosnaps[i]].Halo[hostindexes[i]].x;
		orbitinghalopos = snapdata[halosnaps[i]].Halo[haloindexes[i]].x+boundry;

		if(tmppos-orbitinghalopos>0.5*Cosmo.boxsize) tmppos -=Cosmo.boxsize;
		if(tmppos-orbitinghalopos<-0.5*Cosmo.boxsize) tmppos +=Cosmo.boxsize;

		x[i] = tmppos;
		if(i<nhalo-1){
			orbitinghalonextpos = snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].x+boundry;
			if(orbitinghalonextpos-orbitinghalopos>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(orbitinghalonextpos-orbitinghalopos<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.x, halouniages, x, nhalo);


	/*  y-pos  */

	//Reset the boundary
	boundry=0;

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		tmppos = snapdata[halosnaps[i]].Halo[hostindexes[i]].y;
		orbitinghalopos = snapdata[halosnaps[i]].Halo[haloindexes[i]].y+boundry;

		if(tmppos-orbitinghalopos>0.5*Cosmo.boxsize) tmppos -=Cosmo.boxsize;
		if(tmppos-orbitinghalopos<-0.5*Cosmo.boxsize) tmppos +=Cosmo.boxsize;

		y[i] = tmppos;
		if(i<nhalo-1){
			orbitinghalonextpos = snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].y+boundry;
			if(orbitinghalonextpos-orbitinghalopos>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(orbitinghalonextpos-orbitinghalopos<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}


	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.y, halouniages, y, nhalo);


	/*  z-pos  */

	//Reset the boundary
	boundry=0;

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++){
		tmppos = snapdata[halosnaps[i]].Halo[hostindexes[i]].z;
		orbitinghalopos = snapdata[halosnaps[i]].Halo[haloindexes[i]].z+boundry;

		if(tmppos-orbitinghalopos>0.5*Cosmo.boxsize) tmppos -=Cosmo.boxsize;
		if(tmppos-orbitinghalopos<-0.5*Cosmo.boxsize) tmppos +=Cosmo.boxsize;

		z[i] = tmppos;
		if(i<nhalo-1){
			orbitinghalonextpos = snapdata[halosnaps[i+1]].Halo[haloindexes[i+1]].z+boundry;
			if(orbitinghalonextpos-orbitinghalopos>0.5*Cosmo.boxsize) boundry-=Cosmo.boxsize;
			else if(orbitinghalonextpos-orbitinghalopos<-0.5*Cosmo.boxsize) boundry+=Cosmo.boxsize;
		}
	}


	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.z, halouniages, z, nhalo);

	/*  x-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vx[i] = snapdata[halosnaps[i]].Halo[hostindexes[i]].vx;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vx, halouniages, vx, nhalo);

	/*  y-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vy[i] = snapdata[halosnaps[i]].Halo[hostindexes[i]].vy;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vy, halouniages, vy, nhalo);

	/*  z-vel  */

	//Lets extract the data for the interpolation routine
	for(int i = 0;i<nhalo;i++)
		vz[i] = snapdata[halosnaps[i]].Halo[hostindexes[i]].vz;

	//Intialize the data for the spline
	gsl_spline_init (splinefuncs.vz, halouniages, vz, nhalo);
}

void InterpSingleHaloProps(double interpuniage, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata, SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs){

	double xhost, yhost, zhost, vxhost, vyhost, vzhost;
	double f = (interpuniage - prevuniage)/(currentuniage - prevuniage);

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

	//Interpolate the posistions from the spline functions
	tmporbitdata.x = gsl_spline_eval(splinefuncs.x,interpuniage,splinefuncs.xacc);
	tmporbitdata.y = gsl_spline_eval(splinefuncs.y,interpuniage,splinefuncs.yacc);
	tmporbitdata.z = gsl_spline_eval(splinefuncs.z,interpuniage,splinefuncs.zacc);
	xhost = gsl_spline_eval(hostsplinefuncs.x,interpuniage,hostsplinefuncs.xacc);
	yhost = gsl_spline_eval(hostsplinefuncs.y,interpuniage,hostsplinefuncs.yacc);
	zhost = gsl_spline_eval(hostsplinefuncs.z,interpuniage,hostsplinefuncs.zacc);

	tmporbitdata.xrel = xhost - tmporbitdata.x;
	tmporbitdata.yrel = yhost - tmporbitdata.y;
	tmporbitdata.zrel = zhost - tmporbitdata.z;

	if(tmporbitdata.xrel>0.5*Cosmo.boxsize)
		cout<<"Have a halo greater than the boxsize apart "<<tmporbitdata.orbitID<<endl;


	//Interpolate the velocities from the spline functions
	tmporbitdata.vx = gsl_spline_eval(splinefuncs.vx,interpuniage,splinefuncs.vxacc);
	tmporbitdata.vy = gsl_spline_eval(splinefuncs.vy,interpuniage,splinefuncs.vyacc);
	tmporbitdata.vz = gsl_spline_eval(splinefuncs.vz,interpuniage,splinefuncs.vzacc);
	vxhost = gsl_spline_eval(hostsplinefuncs.vx,interpuniage,hostsplinefuncs.vxacc);
	vyhost = gsl_spline_eval(hostsplinefuncs.vy,interpuniage,hostsplinefuncs.vyacc);
	vzhost = gsl_spline_eval(hostsplinefuncs.vz,interpuniage,hostsplinefuncs.vzacc);

	tmporbitdata.vxrel = vxhost - tmporbitdata.vx;
	tmporbitdata.vyrel = vyhost - tmporbitdata.vy;
	tmporbitdata.vzrel = vzhost - tmporbitdata.vz;

	// r = sqrt((tmporbitdata.x-xhost)*(tmporbitdata.x-xhost) + (tmporbitdata.y-yhost)*(tmporbitdata.y-yhost) + (tmporbitdata.z-zhost)*(tmporbitdata.z-zhost));

	// cout<<"Returned "<<r/tmporbitdata.rvirhost<<" "<<r<<" "<<tmporbitdata.rvirhost<<" "<<interpuniage<<endl;

}


double InterpCrossingHaloProps(float numrvircrossing, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata, SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs){


	//Set the number of points to sample the interpolation functions
	Int_t ninterp = 100;
	double f, stepsize, x, y, z, xhost, yhost, zhost, rvirhost, interpuniages[ninterp], interpvals[ninterp];
	int index = 0;

	//Set the stepsize
	stepsize = (currentuniage-prevuniage)/(double)(ninterp);

	//Set the first interpuniages to the previous uniage and evaluate the spline at this point
	interpuniages[0] = prevuniage;
	x = gsl_spline_eval(splinefuncs.x,prevuniage,splinefuncs.xacc);
	y = gsl_spline_eval(splinefuncs.y,prevuniage,splinefuncs.yacc);
	z = gsl_spline_eval(splinefuncs.z,prevuniage,splinefuncs.zacc);
	xhost = gsl_spline_eval(hostsplinefuncs.x,prevuniage,hostsplinefuncs.xacc);
	yhost = gsl_spline_eval(hostsplinefuncs.y,prevuniage,hostsplinefuncs.yacc);
	zhost = gsl_spline_eval(hostsplinefuncs.z,prevuniage,hostsplinefuncs.zacc);

	rvirhost = LogInterp(prevhosthalo.rvir,hosthalo.rvir,0);

	interpvals[0] = sqrt((x-xhost)*(x-xhost) + (y-yhost)*(y-yhost) + (z-zhost)*(z-zhost))/rvirhost - (double)abs(numrvircrossing);

	for(int i = 1; i<ninterp; i++){
		interpuniages[i] = interpuniages[i-1] + stepsize;

		//Finde the interpolated positions at this point.
		x = gsl_spline_eval(splinefuncs.x,interpuniages[i],splinefuncs.xacc);
		y = gsl_spline_eval(splinefuncs.y,interpuniages[i],splinefuncs.yacc);
		z = gsl_spline_eval(splinefuncs.z,interpuniages[i],splinefuncs.zacc);
		xhost = gsl_spline_eval(hostsplinefuncs.x,interpuniages[i],hostsplinefuncs.xacc);
		yhost = gsl_spline_eval(hostsplinefuncs.y,interpuniages[i],hostsplinefuncs.yacc);
		zhost = gsl_spline_eval(hostsplinefuncs.z,interpuniages[i],hostsplinefuncs.zacc);

		f = (interpuniages[i] - prevuniage)/(currentuniage - prevuniage);

		rvirhost = LogInterp(prevhosthalo.rvir,hosthalo.rvir,f);

		//Find the difference to the crossing point of interest 
		interpvals[i] = sqrt((x-xhost)*(x-xhost) + (y-yhost)*(y-yhost) + (z-zhost)*(z-zhost))/rvirhost - (double)abs(numrvircrossing);
	}


	//Now find the index of the smallest element
    for(int i = 1; i < ninterp; i++) if(abs(interpvals[i]) < abs(interpvals[index]))
            index = i;

	//Now the corresponding time can be used to interpolate the rest of the halo's properties
	InterpSingleHaloProps(interpuniages[index], currentuniage, prevuniage, orbitinghalo, hosthalo, prevorbitinghalo, prevhosthalo, tmporbitdata, snapdata, splinefuncs, hostsplinefuncs);

	return interpuniages[index];
}

HaloData InterpHaloProps(Options &opt, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<Int_t> &interpsnaps, vector<SnapData> &snapdata, SplineFuncs &splinefuncs){

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

	// InterpHaloPosVel(nhalo,ninterp,halouniages,interpuniages,halosnaps,haloindexes,snapdata,interphalos);

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

			//Lets first use the cubic spline to get the positions and velocties
			interphalos[j].x = gsl_spline_eval(splinefuncs.x,snapdata[currentsnap].uniage,splinefuncs.xacc);
			interphalos[j].y = gsl_spline_eval(splinefuncs.y,snapdata[currentsnap].uniage,splinefuncs.yacc);
			interphalos[j].z = gsl_spline_eval(splinefuncs.z,snapdata[currentsnap].uniage,splinefuncs.zacc);
			interphalos[j].vx = gsl_spline_eval(splinefuncs.vx,snapdata[currentsnap].uniage,splinefuncs.vxacc);
			interphalos[j].vy = gsl_spline_eval(splinefuncs.vy,snapdata[currentsnap].uniage,splinefuncs.vyacc);
			interphalos[j].vz = gsl_spline_eval(splinefuncs.vz,snapdata[currentsnap].uniage,splinefuncs.vzacc);

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
				cout<<"Warning: this halo is set to orbit a host halo thats at a different snapshot, have the host halos been interpolated?"<<endl;

			//Iterate the snapshot
			currentsnap++;
			j++;
		}
	}
	//Now lets add this interpolated halo into the snapsdata at each of the interpolation
	//snapshots and also add to the number of halos at each snapshot
	for(int i=0;i<ninterp;i++){
		snapdata[interpsnaps[i]].Halo.push_back(interphalos[i]);

		// //insert these halos into the halosnaps and indexes
		halosnaps.insert(halosnaps.begin()+interpsnaps[i]-halosnaps[0],interpsnaps[i]);
		haloindexes.insert(haloindexes.begin()+interpsnaps[i]-halosnaps[0],snapdata[interpsnaps[i]].numhalos);
		snapdata[interpsnaps[i]].numhalos++;
	}
}