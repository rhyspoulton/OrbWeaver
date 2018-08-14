import VELOCIraptor_Python_Tools.velociraptor_python_tools as VPT 
import SetOrbitForestID as SOF 
import numpy as np
import scipy.spatial as spatial
import time
from scipy.interpolate import interp1d

np.set_printoptions(threshold=10)

atime,tree,numhalos,halodata,cosmodata,unitdata = VPT.ReadUnifiedTreeandHaloCatalog("/mnt/su3ctm/pelahi/waves/analysis/waves_40_512/oldtree/VELOCIraptor.tree.t4.unifiedhalotree",desiredfields=["ID","Head","Tail","RootTail","RootHead","npart","Mass_200crit","R_200crit","Xc","Yc","Zc","VXc","VYc","VZc","Vmax","Rmax"])

atime = atime[::-1]
numhalos = numhalos[::-1]
halodata = halodata[::-1]


numsnaps = 200

numRvirSearch = 4

NpartLim = 100000
MinSnapExist = 20
TEMPORALHALOIDVAL = 1000000000000
iverbose = 1
numOrbitalForestPerFile = 50
outfilebasename = "/mnt/su3ctm/rpoulton/orbitdata/cat"

# Build a new data stucture to contain the information for this file
treefields = ["origID","ID","Head","Tail","OrbitingHaloID","hostHaloID"]
orbitalfields = ["Mass_200crit","R_200crit","npart","Xc","Yc","Zc","VXc","VYc","VZc","Vmax","Rmax"]
orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(numsnaps)]

#initialize the dictionaries
for snap in range(numsnaps):
	#store id and snap and mass of last major merger and while we're at it, store number of major mergers
	halodata[snap]["OrbitingHaloID"] = np.zeros(numhalos[snap],dtype=np.int64)

# built KD tree to quickly search for near neighbours. only build if not passed.
start=time.clock()
boxsize=cosmodata['BoxSize']
hval=cosmodata['Hubble_param']
pos=[[]for j in range(numsnaps)]
pos_tree=[[]for j in range(numsnaps)]
start=time.clock()
if (iverbose): print("KD tree build")
for snap in range(numsnaps-1,-1,-1):
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
for j in range(numsnaps-1,-1,-1):
	start2=time.clock()
	if (numhalos[j]==0): continue
	#First define halos of interest, intially just do it based on mass and how long the halo has existed for
	haloIndexes = np.where((halodata[j]["npart"]>NpartLim) & ((halodata[j]["RootHead"]/TEMPORALHALOIDVAL-halodata[j]["RootTail"]/TEMPORALHALOIDVAL).astype(int)>=MinSnapExist))[0]
	#Loop over all the interesting halos
	for indx in haloIndexes:
		#Skip if we have already set this halos OrbitForestID
		if(halodata[j]["OrbitingHaloID"][indx]!=0): continue

		start3 = time.clock()
		#Set the OrbitalForestID
		print("On forest",orbitforestidval)
		SOF.SetOrbitalForestID(j,numsnaps,halodata,halodata[j]["ID"][indx],orbitforestidval,orbitdata,treefields,orbitalfields,pos_tree)

		
		inumForest+=1

		if(inumForest==numOrbitalForestPerFile):
			SOF.OutputOrbitalForestIDFile(numsnaps,outfilebasename,orbitdata,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata)

			inumForest = 0
			prevorbitforestidval = orbitforestidval +1

			#Reset the orbitdat
			orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(numsnaps)]

		#Iterate the orbitforestidval
		orbitforestidval+=1


	if (iverbose): print("Done snap",j,time.clock()-start2)

if(orbitforestidval!=prevorbitforestidval -1):
	SOF.OutputOrbitalForestIDFile(numsnaps,outfilebasename,orbitdata,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata)
print("Done generating forest",time.clock()-start)


# #get the size of each forest
# OrbitForestSize=np.zeros(orbitforestidval,dtype=np.int64)
# for j in range(numsnaps):
# 	if (numhalos[j]==0): continue
# 	uniqueforest,counts=np.unique(halodata[j]['OrbitForestID'],return_counts=True)
# 	for icount in range(len(uniqueforest)):
# 		OrbitForestSize[uniqueforest[icount]-1]+=counts[icount]
# 	# if (iverbose): print("Finished processing forest size for snap",j)

# for i in range(orbitforestidval):
# 	print("Orbit Forest ID",i,"has size",OrbitForestSize[i])

