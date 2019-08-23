import sys, os, glob
import numpy as np
import scipy.spatial as spatial
import time
from scipy.interpolate import interp1d
import argparse

#Load the other routines
scriptpath=os.path.abspath(__file__)
basecodedir=scriptpath.split('CreateOrbitCatalogs.py')[0]
sys.path.append(basecodedir+'/src/')

from MakeOrbitForest import CreateOrbitForest
from ui import GetDatasetNames,Options
from orbio import ReadVELOCIraptorTreeandHalodata,OutputOrbitCatalog


parser = argparse.ArgumentParser()
parser.add_argument("-c",action="store",dest="configfile",help="Configuration file (orbweaver.cfg)",required=True)
parser.add_argument("-t",action="store",dest="inputtreefilename",help="The VELOCIraptor walkabletree file",required=True)
parser.add_argument("-i",action="store",dest="inputhalofilelistname",help="The filelist containing a list of the base filenames for VELOCIraptor output (e.g. /path/to/VELOCIraptor/output/snapshot_###.VELOCIraptor)",required=True)
parser.add_argument("-o",action="store",dest="outfilebasename",help="The base name for the output catalogs",required=True)
tmpOpt = parser.parse_args()

#store number of snaps
opt = Options(tmpOpt)

#Get the name of the datasets for the desired input format
inputfields = GetDatasetNames(opt,basecodedir)

#Load in the VELOCIraptor halodata and tree
atime, numhalos, halodata, tree, unitdata, cosmodata = ReadVELOCIraptorTreeandHalodata(opt,inputfields)


# built KD tree to quickly search for near neighbours. only build if not passed.
start=time.clock()
pos_tree=[[] for j in range(opt.numsnaps)]
start=time.clock()
if (opt.iverbose): print("KD tree build")
for snap in range(opt.numsnaps-1,-1,-1):
	if (numhalos[snap]>0):
		start1 = time.clock()
		if(opt.iverbose>1): print('Snapshot', snap, 'producing spatial tree')
		pos=np.transpose(np.asarray([halodata[snap]["X"],halodata[snap]["Y"],halodata[snap]["Z"]]))
		pos_tree[snap]=spatial.cKDTree(pos,boxsize=halodata[snap]["SimulationInfo"]["Period"])
		if(opt.iverbose>1): print('Done',snap,'in',time.clock()-start1)
if (opt.iverbose): print("Done building in",time.clock()-start)
sys.stdout.flush()


#Extract all the datatypes from the numpy arrays
datatypes = {field:halodata[0][field].dtype for field in inputfields.keys()}
datatypes.update({field:tree[0][field].dtype for field in tree[0].keys()})

# Build a new data stucture to contain the information for this file
treefields = ["OrigID","OrigRootProgenID","OrigRootDescenID","ID","Head","Tail","OrbitedHaloID","FieldHalo","RatioOfMassinSubsStruct","hostMerges"]
orbitalfields = [field for field in inputfields.keys() if field!="host_id"]
orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(opt.numsnaps)]

orbithalocount = np.zeros(opt.numOrbitForestPerFile, dtype=np.int64)
orbitforestfields = ['Number_of_halos', 'Main_branch_length']
orbitforestdata = {orbitforestfield:[] for orbitforestfield in orbitforestfields}


#Add in the extra datatypes for the extra orbit fields
datatypes["OrigID"] =np.dtype("uint64")
datatypes["OrigRootProgenID"] = np.dtype("int64")
datatypes["OrigRootDescenID"] = np.dtype("int64")
datatypes["FieldHalo"] = np.dtype("bool")
datatypes["OrbitedHaloID"] = np.dtype("int64")
datatypes["RatioOfMassinSubsStruct"] = np.dtype("float32")
datatypes["hostMerges"] = np.dtype("bool")

#initialize the dictionaries and also generate the HeadIndex for the halos as this offers a significant speed up if computed here
for snap in range(opt.numsnaps):
	#Set a done flag for the halos who's orbits around have already be analysed
	halodata[snap]["doneFlag"] = np.zeros(numhalos[snap],dtype=bool)
	halodata[snap]["MassinSubStruct"] = np.zeros(numhalos[snap],dtype=np.float32)

#Lets pre-compute the mass in substructre for all halos
for snap in range(opt.numsnaps):
	indexes = np.where(halodata[snap]["host_id"]>-1)[0]
	hostIndxes = np.array(halodata[snap]["host_id"][indexes]%opt.TEMPORALHALOIDVAL-1,dtype=np.int64)
	for i in range(indexes.size):
		halodata[snap]["MassinSubStruct"][hostIndxes[i]]+=halodata[snap]["Mass"][indexes[i]]


if(opt.iverbose): print("Building the orbit forests")

#Now walk backwards in time along the halos history finding any unique
#branch that comes within numRvirSearch * Radius and set it to 
#be a member of this orbital forest ID
orbitforestidval=0
start=time.clock()
inumForest = 0
prevorbitforestidval=0
ifileno=0
for j in range(opt.numsnaps-1,-1,-1):
	start2=time.clock()
	if (numhalos[j]==0): continue
	#First define halos of interest, intially just do it based on mass and how long the halo has existed for
	haloIndexes = np.where((halodata[j]["npart"]>opt.NpartLimHost) & ((tree[j]["RootHeadSnap"]-tree[j]["RootTailSnap"])>=opt.MinNumSnapExistHost))[0]
	if(opt.iverbose>1): print('Snapshot',j,' containing initial set of ', haloIndexes.size, 'orbital forest candidates') 
	sys.stdout.flush()

	#Loop over all the interesting halos
	for indx in haloIndexes:

		#Skip if this halo's orbits have already been extracted
		if(halodata[j]["doneFlag"][indx]): continue

		start3 = time.clock()
		#Set the OrbitalForestID
		if (opt.iverbose > 1):
			print("On oribital forest", orbitforestidval)
			sys.stdout.flush()

		orbithalocount[inumForest] = CreateOrbitForest(opt,numhalos,halodata,tree,tree[j]["ID"][indx],orbitforestidval,orbitdata,atime,treefields,orbitalfields,pos_tree,cosmodata)

		#If there were zero entries in this halos orbit forest
		if(orbithalocount[inumForest]==0):
			continue

		if (opt.iverbose > 1):
			print("Done orbital forest", orbitforestidval, time.clock()-start3)
			sys.stdout.flush()

		#Keep track of the current number of forest
		inumForest+=1

		if(inumForest == opt.numOrbitForestPerFile):
			#get orbit forest statistic
			orbitforestdata['Number_of_halos']=np.array(orbithalocount, dtype=np.int64)

			OutputOrbitCatalog(opt,orbitdata,datatypes,
				prevorbitforestidval,orbitforestidval,inumForest,ifileno,
				atime,cosmodata,unitdata,
				orbitforestdata)

			inumForest = 0
			prevorbitforestidval = orbitforestidval +1
			ifileno+=1

			#Reset the orbitdat
			orbitdata = [{field:[] for field in orbitalfields+treefields} for snap in range(opt.numsnaps)]

			#Reset the orbit forsest info 
			orbitforestdata = {orbitforestfield:[] for orbitforestfield in orbitforestfields}

		#Iterate the orbitforestidval
		orbitforestidval+=1


	if (opt.iverbose>1):
		print("Done snap",j,time.clock()-start2)
		sys.stdout.flush()


if(orbitforestidval!=prevorbitforestidval -1):
	orbitforestdata['Number_of_halos']=np.array(orbithalocount[:inumForest], dtype=np.int64)
	OutputOrbitCatalog(opt,orbitdata,datatypes,prevorbitforestidval,orbitforestidval,inumForest,ifileno,atime,cosmodata,unitdata,orbitforestdata)
print("Done generating orbit forest",time.clock()-start)
sys.stdout.flush()


#Create a filelist containing all the filenames
ifile = open(opt.outfilebasename + ".orbweaver.filelist.txt","w")

#Output the number of files
ifile.write("%i\n" %(ifileno+1))

for i in range(ifileno+1):

	#If on the final fileno then output a line without a newline
	if(i==ifileno):
		ifile.write(opt.outfilebasename + ".%i" %i)
	else:
		ifile.write(opt.outfilebasename + ".%i\n" %i)

ifile.close()

print("The file containing the list of the catalogues is here:\n\t",opt.outfilebasename + ".orbweaver.filelist.txt")

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

