

#include "orbweaver.h"


void CleanOrbits(Options &opt, vector<OrbitData> &branchorbitdata, vector<OrbitProps> passagesorbitprops){
	/*

	This function cleans the orbits of "wobbles" where passages seem to happen due to an interaction with another object.
	The cleaning is done by removing any passage with phi and orbitecc_ratio that are below a reigion in the Phi and
	orbiteccratio plane. This region is set by PHICLEANLIMIT and ECCCLEANLIMI set in allvars.h.

	When removing passages it also requires updating the number of orbits for the crossing points between passages and
	any points after the deleted passaged. In addition the orbit averaged quanties

	*/

	int passageindex, prevpassageindex, nextpassageindex, updatepassageindex;
	int currpassagensnap, prevpassagesnap;
	int numpassages = passagesorbitprops.size();
	double ecclimit, hostlx, hostly, hostlz, rcomove;
	int i=0,numoutputs = branchorbitdata.size();


	//Store a vector of indices to delete
	vector<int> idel;

	//Loop over all the passages in the passage orbit props
	while(i<passagesorbitprops.size()){

		//Extract the index which this is in the branchorbitdata
		passageindex = passagesorbitprops[i].passageindex;


		//Check if this passage is within CLEANRATIOHOSTRADIUS * host radius, if so then lets continue to the next passage.
		//This is so orbits which have not been properly sampled by the simulation are also included if desired
		if((branchorbitdata[passageindex].r/branchorbitdata[passageindex].rvirhost)<CLEANRATIOHOSTRADIUS){
			i++;
			continue;
		}

		//Lets see if the eccentricity and Phi are below the region as described in Poulton et at., in prep
		ecclimit = ECCCLEANLIMIT - (ECCCLEANLIMIT/PHICLEANLIMIT) * branchorbitdata[passageindex].phi;
		if((ecclimit>branchorbitdata[passageindex].orbitecc_ratio) & (ecclimit<ECCCLEANLIMIT)){

			//Now we have found it to be less than the clean limits it can be add this index to the delete vector
			idel.push_back(passageindex);

			//Initially add the previous passage as the one to be deleted
			prevpassageindex = passagesorbitprops[i].prevpassageindex;

			//So one passage is to be deleted, but another needed to be deleted since there would be two passages of the same time
			// therefore we need to check which one caused this wobble either the previous or next passage (if it exists)
			if(i+1<numpassages){

				//Extract the index for the next passage
				nextpassageindex = passagesorbitprops[i+1].passageindex;

				//Lets check if the next is worse than this one
				if((branchorbitdata[nextpassageindex].phi<branchorbitdata[passageindex].phi) & (branchorbitdata[nextpassageindex].orbitecc_ratio<branchorbitdata[passageindex].orbitecc_ratio)){

					//The next passage is wrose to add to the delete vector
					idel.push_back(nextpassageindex);

					//Lets first update the number of orbits we only need to subtract 0.5 an orbit between the passage index and the next passage index
					for(int j=passageindex; j<nextpassageindex; j++)
						branchorbitdata[j].numorbits-=0.5;

					//But we need to remove 1.0 orbit after the next passage index since two passages have been removed
					for(int j=nextpassageindex; j<branchorbitdata.size(); j++)
						branchorbitdata[j].numorbits-=1.0;

					//Lets check if there is another passage, if so then its properties needs to be updated
					if(i+2<numpassages){

						//Lets add the properties of the deleted passages to this passage
						passagesorbitprops[i+2]+=passagesorbitprops[i];
						passagesorbitprops[i+2]+=passagesorbitprops[i+1];

						//Update the previous passage index to point to the passsage before the deleted ones
						passagesorbitprops[i+2].prevpassageindex = passagesorbitprops[i].prevpassageindex;

						//Extract the index for this passage
						updatepassageindex = passagesorbitprops[i+2].passageindex;

						//Find the current and previous passages for snapshots
						currpassagensnap = passagesorbitprops[i+2].passagesnap;
						prevpassagesnap = passagesorbitprops[i].prevpassagesnap;

						//Need to update the orbit eccentricity ratio for this passages
						rcomove = branchorbitdata[updatepassageindex].r * Cosmo.h / branchorbitdata[updatepassageindex].scalefactor;
						passagesorbitprops[i+2].prevpassageindex = passagesorbitprops[i].prevpassageindex;if(branchorbitdata[updatepassageindex].entrytype==99)
							branchorbitdata[updatepassageindex].orbitecc_ratio = (passagesorbitprops[i].prevpassagercomove-rcomove)/(passagesorbitprops[i].prevpassagercomove+rcomove);
						else
							branchorbitdata[updatepassageindex].orbitecc_ratio = (rcomove-passagesorbitprops[i].prevpassagercomove)/(rcomove+passagesorbitprops[i].prevpassagercomove);

						//We need to re-calculate the orbit averaged quanties
						//Find the average energy and reduced mass
						branchorbitdata[updatepassageindex].orbitalenergy_ave = passagesorbitprops[i+2].E / (double)(currpassagensnap - prevpassagesnap);;

						//The average mass loss rate
						branchorbitdata[updatepassageindex].masslossrate_ave = passagesorbitprops[i+2].masslossrate / (double)(currpassagensnap -prevpassagesnap);

						//Lets output the average angular momentum since the last passage
						branchorbitdata[updatepassageindex].lxrel_ave = passagesorbitprops[i+2].lx / (double)(currpassagensnap - prevpassagesnap);
						branchorbitdata[updatepassageindex].lyrel_ave = passagesorbitprops[i+2].ly / (double)(currpassagensnap - prevpassagesnap);
						branchorbitdata[updatepassageindex].lzrel_ave = passagesorbitprops[i+2].lz / (double)(currpassagensnap - prevpassagesnap);

						//Find the average angular momentum of the host halo
						hostlx = passagesorbitprops[i+2].hostlx / (double)(currpassagensnap - prevpassagesnap);
						hostly = passagesorbitprops[i+2].hostly / (double)(currpassagensnap - prevpassagesnap);
						hostlz = passagesorbitprops[i+2].hostlz / (double)(currpassagensnap - prevpassagesnap);

						//Now we have the average angular momentum,the alignment with the host angular momentum can be computed
						branchorbitdata[updatepassageindex].hostalignment = acos(((hostlx*branchorbitdata[updatepassageindex].lxrel_ave) + (hostly*branchorbitdata[updatepassageindex].lyrel_ave)	+ (hostlz*branchorbitdata[updatepassageindex].lzrel_ave))/
							(sqrt(hostlx*hostlx + hostly*hostly + hostlz*hostlz) * sqrt(branchorbitdata[updatepassageindex].lxrel_ave*branchorbitdata[updatepassageindex].lxrel_ave + branchorbitdata[updatepassageindex].lyrel_ave*branchorbitdata[updatepassageindex].lyrel_ave + branchorbitdata[updatepassageindex].lzrel_ave*branchorbitdata[updatepassageindex].lzrel_ave)));

						//Update the value of phi for this orbit
						branchorbitdata[updatepassageindex].phi = passagesorbitprops[i+2].phi;

						//Now we have updated the orbit averaged quantities we also need to update the num_entrytype,
						// only 1 needs to be removed from each num_entrytype since 1 pericenter and 1 apocentre have been deleted
						for(int j=i+2; j<passagesorbitprops.size(); j++){

							//Extract the index of the passage
							passageindex = passagesorbitprops[j].passageindex;

							//Remove 1 from its num_entrytype
							branchorbitdata[passageindex].num_entrytype-=1;
						}

					}
					//Otherwise nothing needs to be done

					//Need to add to the iterator here as the next passage has already been deleted
					i++;

				}
				else{
					// If there is another passage but it the phi and eccentricity is not worse then we need to update a different passage


					//If there is no next passage then we just need to delete the previous passage
					idel.push_back(prevpassageindex);

					//Then update the number of orbits between the passages, first update between the passages and subtract 0.5 an orbit
					for(int j=prevpassageindex; j<passageindex; j++)
						branchorbitdata[j].numorbits-=0.5;

					//Next update after the passages and remove 1.0 orbit (two deleted passages)
					for(int j=passageindex; j<branchorbitdata.size(); j++)
						branchorbitdata[j].numorbits-=1.0;

					//Extract the index for this passage
					updatepassageindex = passagesorbitprops[i+1].passageindex;

					//We need to update the orbit properties for the new apsis points, this only needs to be done when i>0
					//i.e. the passage to be updated is not first passage
					if((i>0) & (idel.size()<i+2)){
						//Mark this passage as beeb deleted

						//Lets add the properties of the deleted passages to this passage
						passagesorbitprops[i+1]+=passagesorbitprops[i-1];
						passagesorbitprops[i+1]+=passagesorbitprops[i];

						//Find the current and previous passages for snapshots
						currpassagensnap = passagesorbitprops[i+1].passagesnap;
						prevpassagesnap = passagesorbitprops[i-1].prevpassagesnap;

						//Update the previous passage index to point to the passsage before the deleted ones
						passagesorbitprops[i+1].prevpassageindex = passagesorbitprops[i-1].prevpassageindex;

						//Need to update the orbit eccentricity ratio for this passages
						rcomove = branchorbitdata[updatepassageindex].r * Cosmo.h / branchorbitdata[updatepassageindex].scalefactor;
						if(branchorbitdata[updatepassageindex].entrytype==99)
							branchorbitdata[updatepassageindex].orbitecc_ratio = (passagesorbitprops[i-1].prevpassagercomove-rcomove)/(passagesorbitprops[i-1].prevpassagercomove+rcomove);
						else
							branchorbitdata[updatepassageindex].orbitecc_ratio = (rcomove-passagesorbitprops[i-1].prevpassagercomove)/(rcomove+passagesorbitprops[i-1].prevpassagercomove);

						//We need to re-calculate the orbit averaged quanties
						//Find the average energy and reduced mass
						branchorbitdata[updatepassageindex].orbitalenergy_ave = passagesorbitprops[i+1].E / (double)(currpassagensnap - prevpassagesnap);;

						//The average mass loss rate
						branchorbitdata[updatepassageindex].masslossrate_ave = passagesorbitprops[i+1].masslossrate / (double)(currpassagensnap -prevpassagesnap);

						//Lets output the average angular momentum since the last passage
						branchorbitdata[updatepassageindex].lxrel_ave = passagesorbitprops[i+1].lx / (double)(currpassagensnap - prevpassagesnap);
						branchorbitdata[updatepassageindex].lyrel_ave = passagesorbitprops[i+1].ly / (double)(currpassagensnap - prevpassagesnap);
						branchorbitdata[updatepassageindex].lzrel_ave = passagesorbitprops[i+1].lz / (double)(currpassagensnap - prevpassagesnap);

						//Find the average angular momentum of the host halo
						hostlx = passagesorbitprops[i+1].hostlx / (double)(currpassagensnap - prevpassagesnap);
						hostly = passagesorbitprops[i+1].hostly / (double)(currpassagensnap - prevpassagesnap);
						hostlz = passagesorbitprops[i+1].hostlz / (double)(currpassagensnap - prevpassagesnap);

						//Now we have the average angular momentum,the alignment with the host angular momentum can be computed
						branchorbitdata[updatepassageindex].hostalignment = acos(((hostlx*branchorbitdata[updatepassageindex].lxrel_ave) + (hostly*branchorbitdata[updatepassageindex].lyrel_ave)	+ (hostlz*branchorbitdata[updatepassageindex].lzrel_ave))/
							(sqrt(hostlx*hostlx + hostly*hostly + hostlz*hostlz) * sqrt(branchorbitdata[updatepassageindex].lxrel_ave*branchorbitdata[updatepassageindex].lxrel_ave + branchorbitdata[updatepassageindex].lyrel_ave*branchorbitdata[updatepassageindex].lyrel_ave + branchorbitdata[updatepassageindex].lzrel_ave*branchorbitdata[updatepassageindex].lzrel_ave)));

						//Update the value of phi for this orbit
						branchorbitdata[updatepassageindex].phi = passagesorbitprops[i+1].phi;
					}
					else{
						//Otherwise the passage to be updated is the first passage in the orbit so should not have the above quantities calculated
						//as there is not another passage which they can be calculated with reference to
						branchorbitdata[updatepassageindex].orbitecc_ratio = 0;
						branchorbitdata[updatepassageindex].orbitalenergy_ave = 0;
						branchorbitdata[updatepassageindex].masslossrate_ave = 0;
						branchorbitdata[updatepassageindex].lxrel_ave = 0;
						branchorbitdata[updatepassageindex].lyrel_ave = 0;
						branchorbitdata[updatepassageindex].lzrel_ave = 0;
						branchorbitdata[updatepassageindex].hostalignment = 0;
						branchorbitdata[updatepassageindex].phi = 0;
					}


					//Now we have updated the orbit averaged quantities we also need to update the num_entrytype,
					// only 1 needs to be removed from each num_entrytype since 1 pericenter and 1 apocentre have been deleted
					for(int j=i+1; j<passagesorbitprops.size(); j++){

						//Extract the index of the passage
						passageindex = passagesorbitprops[j].passageindex;

						//Remove 1 from its num_entrytype
						branchorbitdata[passageindex].num_entrytype-=1;
					}

				}

			}
			else{
				//If there is no next passage then we just need to delete the previous passage
				idel.push_back(prevpassageindex);

				//Then update the number of orbits between the passages, first update between the passages and subtract 0.5 an orbit
				for(int j=prevpassageindex; j<passageindex; j++)
					branchorbitdata[j].numorbits-=0.5;

				//Next update after the passages and remove 1.0 orbit (two deleted passages)
				for(int j=passageindex; j<branchorbitdata.size(); j++)
					branchorbitdata[j].numorbits-=1.0;

			}
		}

		//Add to the iterator
		i++;
	}

	//Now lets sort the delete vector so in sorted order
	sort(idel.begin(),idel.end());

	//Check if there is any indexes to delete
	if(idel.size()>0){
		for(int i=idel.size()-1;i>=0;i--){
			branchorbitdata.erase(branchorbitdata.begin()+idel[i]);
		}
	}

	//Remove all the values from the delete vector
	idel.clear();

}