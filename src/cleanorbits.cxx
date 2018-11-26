

#include "orbweaver.h"


void CleanOrbits(vector<OrbitData> &branchorbitdata){

	//Store data and the number of passages
	int prevpassageindex, nextpassageindex;
	int numpassages=0;
	float rx,ry,rz, r, tmprx, tmpry, tmprz, tmpr, prevphi, nextphi;
	int numentries = branchorbitdata.size();

	//Store a vector of indices to delete
	vector<int> idel;

	//Lets clean the passage points by check
	for(int i=0;i<branchorbitdata.size();i++){
		if(abs(branchorbitdata[i].entrytype)==99){

			//Skip the first passage by checking if the orbital period is 0
			if(numpassages>0){

				//Calculate the angle between this passage and the previous
				rx = branchorbitdata[i].xrel;
				ry = branchorbitdata[i].yrel;
				rz = branchorbitdata[i].zrel;
				r = sqrt(rx*rx + ry*ry + rz*rz);

				tmprx = branchorbitdata[prevpassageindex].xrel;
				tmpry = branchorbitdata[prevpassageindex].yrel;
				tmprz = branchorbitdata[prevpassageindex].zrel;
				tmpr = sqrt(tmprx*tmprx + tmpry*tmpry + tmprz*tmprz);

				//Now calculate the size of the angle
				prevphi = acos((rx*tmprx + ry*tmpry + rz*tmprz)/(r*tmpr));

				if(prevphi<M_PI/2.0){

					//Lets delete this passage
					idel.push_back(i);
					numpassages--;

					//Now we need to check which of the surrounding passages to remove
					//based on how much the halo has moved around its host
					nextpassageindex = i+1;

					//Move to the next passage point
					while((nextpassageindex<numentries) & (abs(branchorbitdata[nextpassageindex].entrytype)!=99)) nextpassageindex++;

					//If at the end of the entries the delete the previous passage and
					//the clean can be stopeed
					if(nextpassageindex==numentries){
						idel.push_back(prevpassageindex);
						numpassages--;
						break;
					}

					//Otherwise lets calculate the size of the angle between this passage and the previous
					tmprx = branchorbitdata[nextpassageindex].xrel;
					tmpry = branchorbitdata[nextpassageindex].yrel;
					tmprz = branchorbitdata[nextpassageindex].zrel;
					tmpr = sqrt(tmprx*tmprx + tmpry*tmpry + tmprz*tmprz);

					//Now calculate the size of the angle
					nextphi = acos((rx*tmprx + ry*tmpry + rz*tmprz)/(r*tmpr));

					//Check which angle between the passages is larger
					if(prevphi>nextphi){

						//In this case lets remove the next passage as it has moved through a smaller angel
						idel.push_back(nextpassageindex);
						numpassages--;

						//We can now skip the next passage as it will be deleted
						while((i<numentries) & (i<nextpassageindex+1)) i++;
					}
					else{
						//Other wise remove the previous passage
						idel.push_back(prevpassageindex);
						numpassages--;

						//Now need to update the previous pasage index since it is being deleted
						//if numpassages>0, otherwise if numpassages is 0 then all the reference
						// angles needs to be updated
						if(numpassages>0){

							//Find the index of the previous passage
							prevpassageindex--;
							while(abs(branchorbitdata[prevpassageindex].entrytype)!=99) prevpassageindex--;
						}
					}
					continue;
				}
			}

			//Store the index of the previous passage
			prevpassageindex = i;

			//Add one to the number of passages
			numpassages++;
		}
	}

	//Now lets sort the delete vector so in sorted order
	sort(idel.begin(),idel.end());

	//Check if there is any indexes to delete
	if(idel.size()>0){
		for(int i=idel.size()-1;i>=0;i--)
			branchorbitdata.erase(branchorbitdata.begin()+idel[i]);
	}

}