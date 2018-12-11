import VELOCIraptor_Python_Tools.velociraptor_python_tools as VPT 
import SetOrbitForestID as SOF 
from ui import Options
import numpy as np
import scipy.spatial as spatial
import time
from scipy.interpolate import interp1d
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-c",action="store",dest="configfile",help="Configuration file (wherewolf.cfg)",required=True)
parser.add_argument("-t",action="store",dest="inputtree",help="The VELOCIraptor walkabletree file",required=True)
parser.add_argument("-i",action="store",dest="inputhalobbasename",help="The base name for the VELOCIraptor catalog (/path/to/halo/catalog/snapshot_)",required=True)
parser.add_argument("-o",action="store",dest="outfilebasename",help="The base name for the output catalogs",required=True)
tmpOpt = parser.parse_args()

# np.set_printoptions(threshold=10)
desiredfields = ["hostHaloID","npart","Mass_200crit","R_200crit","Xc","Yc","Zc","VXc","VYc","VZc","Vmax","Rmax","cNFW","Mass_tot","Mass_FOF","Lx","Ly","Lz"]
# # atime,tree,numhalos,halodata,cosmodata,unitdata = VPT.ReadUnifiedTreeandHaloCatalog("/mnt/su3ctm/rpoulton/waves/analysis/waves_40_512/VELOCIraptor.tree.t4.unifiedhalotree",desiredfields=desiredfields)

print("Reading the walkable tree")
tree,numsnaps = VPT.ReadWalkableHDFTree(tmpOpt.inputtree,False)
print("Done reading the walkable tree")

print("Reading in the halo catalog")
numhalos=np.zeros(numsnaps,dtype=np.uint64)
halodata=[dict() for i in range(numsnaps)]
atime=np.zeros(numsnaps)
for i in range(numsnaps):
    halodata[i],numhalos[i]=VPT.ReadPropertyFile(tmpOpt.inputhalobbasename+'%03d.VELOCIraptor'%i, 2, 0, 0, desiredfields)
    atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']
    for key in halodata[i].keys():
        if (key == 'SimulationInfo' or key == 'UnitInfo'): continue
        if (halodata[i][key].dtype==np.float64):
            halodata[i][key] = np.array(halodata[i][key],dtype=np.float32)
print('Finished reading halo properties')

opt = Options(tmpOpt,numsnaps)

unitdata = halodata[0]["UnitInfo"]
cosmodata = halodata[0]["SimulationInfo"]

#Extract all the datatypes from the numpy arrays
datatypes = {field:halodata[0][field].dtype for field in desiredfields}
datatypes.update({field:tree[0][field].dtype for field in tree[0].keys()})


# Build a new data stucture to contain the information for this file
treefields = ["origID","ID","Head","Tail","OrbitingHaloID","hostFlag"]
orbitalfields = [field for field in desiredfields if field!="hostHaloID"]
orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(opt.numsnaps)]

#Add in the extra datatypes for the extra orbit fields
datatypes["origID"] =np.dtype("int64")
datatypes["hostFlag"] = np.dtype("bool")
datatypes["OrbitingHaloID"] = np.dtype("int64")

#initialize the dictionaries
for snap in range(opt.numsnaps):
	#Set a done flag for the halos who's orbits around have already be analysed
	halodata[snap]["doneFlag"] = np.zeros(numhalos[snap],dtype=bool)

VPT.AdjustforPeriod(opt.numsnaps, numhalos, halodata, icomove=0)

# built KD tree to quickly search for near neighbours. only build if not passed.
start=time.clock()
pos_tree=[[]for j in range(opt.numsnaps)]
start=time.clock()
if (opt.iverbose): print("KD tree build")
for snap in range(opt.numsnaps-1,-1,-1):
	if (numhalos[snap]>0):
		pos=np.transpose(np.asarray([halodata[snap]["Xc"],halodata[snap]["Yc"],halodata[snap]["Zc"]]))
		pos_tree[snap]=spatial.cKDTree(pos,boxsize=halodata[snap]["SimulationInfo"]["Period"])
if (opt.iverbose): print("done building in",time.clock()-start)


#Now walk backwards in time along the halos history finding any unique
#branch that comes within numRvirSearch * R_200crit and set it to 
#be a member of this orbital forest ID
orbitforestidval=0
start=time.clock()
inumForest = 0
prevorbitforestidval=0
for j in range(opt.numsnaps-1,-1,-1):
	start2=time.clock()
	if (numhalos[j]==0): continue
	#First define halos of interest, intially just do it based on mass and how long the halo has existed for
	haloIndexes = np.where((halodata[j]["npart"]>opt.NpartLimHost) & ((tree[j]["RootHead"]/opt.TEMPORALHALOIDVAL-tree[j]["RootTail"]/opt.TEMPORALHALOIDVAL).astype(int)>=opt.MinSnapExist))[0]
	#Loop over all the interesting halos
	for indx in haloIndexes:

		#Skip if this halo's orbits have already been extracted
		if(halodata[j]["doneFlag"][indx]): continue

		start3 = time.clock()
		#Set the OrbitalForestID
		print("On oribital forest",orbitforestidval)
		SOF.SetOrbitalForestID(opt,numhalos,halodata,tree,tree[j]["ID"][indx],orbitforestidval,orbitdata,atime,treefields,orbitalfields,pos_tree,cosmodata)

		#Keep track of the current number of forest
		inumForest+=1

		if(inumForest==opt.numOrbitalForestPerFile):
			SOF.OutputOrbitalForestIDFile(opt,outfilebasename,orbitdata,datatypes,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata)

			inumForest = 0
			prevorbitforestidval = orbitforestidval +1

			#Reset the orbitdat
			orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(opt.numsnaps)]

		#Iterate the orbitforestidval
		orbitforestidval+=1


	if (opt.iverbose): print("Done snap",j,time.clock()-start2)

if(orbitforestidval!=prevorbitforestidval -1):
	SOF.OutputOrbitalForestIDFile(opt,orbitdata,datatypes,atime,prevorbitforestidval,orbitforestidval,cosmodata,unitdata)
print("Done generating forest",time.clock()-start)


# #get the size of each forest
# OrbitForestSize=np.zeros(orbitforestidval,dtype=np.int64)
# for j in range(opt.numsnaps):
# 	if (numhalos[j]==0): continue
# 	uniqueforest,counts=np.unique(halodata[j]['OrbitForestID'],return_counts=True)
# 	for icount in range(len(uniqueforest)):
# 		OrbitForestSize[uniqueforest[icount]-1]+=counts[icount]
# 	# if (iverbose): print("Finished processing forest size for snap",j)

# for i in range(orbitforestidval):
# 	print("Orbit Forest ID",i,"has size",OrbitForestSize[i])

