

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
	// cout<<rx<<" "<<ry<<" "<<rz<<endl;
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



void ProcessHalo(Int_t snap, Int_t i, Options &opt, SnapData *&snapdata, vector<vector<OrbitData>> &orbitdata){

	unsigned long long descendantID = snapdata[snap].Halo[i].descendant;
	Int_t descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
	Int_t descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);
	unsigned long long haloID = snapdata[snap].Halo[i].id;
	Int_t halosnap = (Int_t)(haloID/opt.TEMPORALHALOIDVAL);
	Int_t haloindex = (Int_t)(haloID%opt.TEMPORALHALOIDVAL-1);


	//Keep track of this halo's snap and index for interpolation
	vector<Int_t> halosnaps = {halosnap};
	vector<Int_t> haloindexes = {haloindex};

	//Store the index of the halo it is orbiting
	Int_t orbitinghaloindex;

	//Set the data for the interpolation routine
	HaloData interphalo;
	Int_t interporbitinghaloindex;
	unsigned long long ihaloID,idescendantID;
	Int_t ihalosnap,ihaloindex,idescendantsnap,idescendantindex;


	//Keep track of the previous halos halodata and orbitdata
	HaloData prevorbitinghalo = snapdata[halosnap].Halo[haloindex];
	HaloData prevhosthalo;
	OrbitData prevorbitdata = {0};

	//Keep track of what snapshot
	Int_t currentsnap = snap;

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

		//Keep track of the orbiting halo's index index if need to interpolate
		interporbitinghaloindex = orbitinghaloindex;

		

		//Only if the current halosnap is greater than the final snapshot in the list of snapshots
		if(halosnap>halosnaps.back()){

			//Lets add this halo into the vectors for the interpolation routine
			halosnaps.push_back(halosnap);
			haloindexes.push_back(haloindex);
		}

		//lets move onto the descendant
		haloID = descendantID;
		halosnap = descendantsnap;
		haloindex = descendantindex;

		//Extract its descendant
		descendantID = snapdata[halosnap].Halo[haloindex].descendant;
		descendantsnap = (Int_t)(descendantID/opt.TEMPORALHALOIDVAL);
		descendantindex = (Int_t)(descendantID%opt.TEMPORALHALOIDVAL-1);

		// Move to the next snapshot
		currentsnap++;

	}


}

void ProcessOrbits(Options &opt, SnapData *&snapdata, vector<vector<OrbitData>> &orbitdata){

	int numsnaps =  opt.fsnap - opt.isnap+1;
	vector<double> scalefactors;
	scalefactors.resize(numsnaps);

	// Initilize the flag which marks the halo as being processed to false
	// and the vector to store the orbital data
	for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++){

		//Extract the scalefactor at each snapshot
		scalefactors[snap] = snapdata[snap].scalefactor;

		//Initilize a vector and resize it (this intilizes the values to 0)
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
		cout<<"Done processing snap "<<snap<<endl;
		if(done) break;
	}
}