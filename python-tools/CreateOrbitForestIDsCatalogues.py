import VELOCIraptor_Python_Tools.velociraptor_python_tools as VPT 
import SetOrbitForestID as SOF 
import numpy as np
import scipy.spatial as spatial
import time
from scipy.interpolate import interp1d

np.set_printoptions(threshold=10)
desiredfields = ["ID","Head","Tail","RootTail","RootHead","hostHaloID","npart","Mass_200crit","R_200crit","Xc","Yc","Zc","VXc","VYc","VZc","Vmax","Rmax","cNFW","Mass_tot","Mass_FOF","Lx","Ly","Lz"]
atime,tree,numhalos,halodata,cosmodata,unitdata = VPT.ReadUnifiedTreeandHaloCatalog("/mnt/su3ctm/rpoulton/waves/analysis/waves_40_512/VELOCIraptor.tree.t4.unifiedhalotree",desiredfields=desiredfields)

atime = atime[::-1]
numhalos = numhalos[::-1]
halodata = halodata[::-1]

#Extract all the datatypes from the numpy arrays
datatypes = {field:halodata[0][field].dtype for field in desiredfields}

buildparams = {}
buildparams["VERSION"] = 0.10
buildparams["numsnaps"] = 200
buildparams["numRvirSearch"] = 4
buildparams["NpartLimHost"] = 100000
buildparams["MinSnapExist"] = 20
buildparams["TEMPORALHALOIDVAL"] = 1000000000000
buildparams["numOrbitalForestPerFile"] = 2000
iverbose = 1
outfilebasename = "/mnt/su3ctm/rpoulton/orbitdata/cat"

# Build a new data stucture to contain the information for this file
treefields = ["origID","ID","Head","Tail","OrbitingHaloID","hostFlag"]
orbitalfields = ["Mass_200crit","Mass_FOF","R_200crit","npart","Xc","Yc","Zc","VXc","VYc","VZc","Vmax","Rmax","cNFW","Lx","Ly","Lz"]
orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(buildparams["numsnaps"])]

#Add in the extra datatypes for the extra orbit fields
datatypes["origID"] =np.dtype("int64")
datatypes["hostFlag"] = np.dtype("bool")
datatypes["OrbitingHaloID"] = np.dtype("int64")

#initialize the dictionaries
for snap in range(buildparams["numsnaps"]):
	#Set a done flag for the halos who's orbits around have already be analysed
	halodata[snap]["doneFlag"] = np.zeros(numhalos[snap],dtype=bool)

# built KD tree to quickly search for near neighbours. only build if not passed.
start=time.clock()
boxsize=cosmodata['BoxSize']
hval=cosmodata['Hubble_param']
pos=[[]for j in range(buildparams["numsnaps"])]
pos_tree=[[]for j in range(buildparams["numsnaps"])]
start=time.clock()
if (iverbose): print("KD tree build")
for snap in range(buildparams["numsnaps"]-1,-1,-1):
	if (numhalos[snap]>0):
		boxval=boxsize*atime[snap]/hval
		pos[snap]=np.transpose(np.asarray([halodata[snap]["Xc"],halodata[snap]["Yc"],halodata[snap]["Zc"]]))
		pos_tree[snap]=spatial.cKDTree(pos[snap],boxsize=boxval)
if (iverbose): print("done building in",time.clock()-start)


#Now walk backwards in time along the halos history finding any unique
#branch that comes within numRvirSearch * R_200crit and set it to 
#be a member of this orbital forest ID
orbitforestidval=0
start=time.clock()
inumForest = 0
prevorbitforestidval=0
for j in range(buildparams["numsnaps"]-1,-1,-1):
	start2=time.clock()
	if (numhalos[j]==0): continue
	#First define halos of interest, intially just do it based on mass and how long the halo has existed for
	haloIndexes = np.where((halodata[j]["npart"]>buildparams["NpartLimHost"]) & ((halodata[j]["RootHead"]/buildparams["TEMPORALHALOIDVAL"]-halodata[j]["RootTail"]/buildparams["TEMPORALHALOIDVAL"]).astype(int)>=buildparams["MinSnapExist"]))[0]
	#Loop over all the interesting halos
	for indx in haloIndexes:

		#Skip if this halo's orbits have already been extracted
		if(halodata[j]["doneFlag"][indx]): continue

		start3 = time.clock()
		#Set the OrbitalForestID
		print("On oribital forest",orbitforestidval)
		SOF.SetOrbitalForestID(buildparams["numsnaps"],numhalos,halodata,halodata[j]["ID"][indx],orbitforestidval,orbitdata,atime,treefields,orbitalfields,pos_tree,cosmodata,buildparams["TEMPORALHALOIDVAL"],buildparams["numRvirSearch"])

		#Keep track of the current number of forest
		inumForest+=1

		if(inumForest==buildparams["numOrbitalForestPerFile"]):
			SOF.OutputOrbitalForestIDFile(buildparams["numsnaps"],outfilebasename,orbitdata,datatypes,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata,buildparams)

			inumForest = 0
			prevorbitforestidval = orbitforestidval +1

			#Reset the orbitdat
			orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(buildparams["numsnaps"])]

		#Iterate the orbitforestidval
		orbitforestidval+=1


	if (iverbose): print("Done snap",j,time.clock()-start2)

if(orbitforestidval!=prevorbitforestidval -1):
	SOF.OutputOrbitalForestIDFile(buildparams["numsnaps"],outfilebasename,orbitdata,datatypes,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata,buildparams)
print("Done generating forest",time.clock()-start)


# #get the size of each forest
# OrbitForestSize=np.zeros(orbitforestidval,dtype=np.int64)
# for j in range(buildparams["numsnaps"]):
# 	if (numhalos[j]==0): continue
# 	uniqueforest,counts=np.unique(halodata[j]['OrbitForestID'],return_counts=True)
# 	for icount in range(len(uniqueforest)):
# 		OrbitForestSize[uniqueforest[icount]-1]+=counts[icount]
# 	# if (iverbose): print("Finished processing forest size for snap",j)

# for i in range(orbitforestidval):
# 	print("Orbit Forest ID",i,"has size",OrbitForestSize[i])

