import os, sys, glob
import h5py
import numpy as np 
import time


#load local python routines
scriptpath=os.path.abspath(__file__)
basecodedir=scriptpath.split('orbio.py')[0]
sys.path.append(basecodedir+'/VELOCIraptor_Python_Tools/')

#load the cythonized code if compiled
if (len(glob.glob(basecodedir+'/VELOCIraptor_Python_Tools/velociraptor_python_tools_cython.*.so'))==1):
	print('using cython VR+TF toolkit')
	import velociraptor_python_tools_cython as VPT
else:
	print('using python VR+TF toolkit')
	import velociraptor_python_tools as VPT


def ReadVELOCIraptorTreeandHalodata(opt,inputfields):

	start=time.time()
	if(opt.iverbose): print("Reading the walkable tree")
	sys.stdout.flush()
	tree,numsnaps = VPT.ReadWalkableHDFTree(opt.inputtreefilename,False)
	if(opt.iverbose): print("Done reading the walkable tree in",time.time()-start)
	sys.stdout.flush()

	#Update the number of snapshots from the tree
	opt.update_numsnaps(numsnaps)

	#The fields that needs to be converted to physical
	physconvertfields = [inputfields[field] for field in inputfields.keys() if field in ["X","Y","Z","Radius","Lx","Ly","Lz"]]

	#Open up the filelist and extract all the filenames
	snapfilelist = open(opt.inputhalofilelistname,"r")
	snapfilenames = [line.strip() for line in snapfilelist]
	snapfilelist.close()

	#Check the number of files in the VELOCIraptor filelist is the same as whats reported by the walkable tree
	if(len(snapfilenames)!=opt.numsnaps):
		raise IOError("The number of input snapshots as reported by the walkable tree is not the same as the\n number of files in the input halo filelist, please correct this")

	#extract the desiredfields from the inputfields
	desiredfields = list(inputfields.values())

	start=time.time()
	if(opt.iverbose): print("Reading in the halo catalog")
	sys.stdout.flush()
	numhalos=np.zeros(opt.numsnaps,dtype=np.uint64)
	halodata=[dict() for i in range(opt.numsnaps)]
	atime=np.zeros(opt.numsnaps)
	for i in range(opt.numsnaps):
		start1 = time.time()
		halodata[i],numhalos[i]=VPT.ReadPropertyFile(snapfilenames[i], 2, 0, 0, desiredfields)
		atime[i]=halodata[i]['SimulationInfo']['ScaleFactor']

		for key in halodata[i].keys():
			if (key == 'SimulationInfo' or key == 'UnitInfo' or key == 'ConfigurationInfo'): continue

			#Reduce the precision of the datasets since OrbWeaver is all in float32
			if (halodata[i][key].dtype==np.float64):
				halodata[i][key].astype(np.float32,casting="same_kind",copy=False)

			#Lets convert fields to physical if required
			if((key in physconvertfields) & (halodata[i]['UnitInfo']["Comoving_or_Physical"]==1)):
				halodata[i][key] *= atime[i]/halodata[i]["SimulationInfo"]["h_val"]

		#Set the flag to physical
		halodata[i]['UnitInfo']["Comoving_or_Physical"] = 0

		if(opt.iverbose > 1): print('Snapshot', i,'done in', time.time()-start1)
		sys.stdout.flush()

	if(opt.iverbose): print('Finished reading halo properties in', time.time()-start)
	sys.stdout.flush()

	VPT.AdjustforPeriod(opt.numsnaps, numhalos, halodata,tree)

	unitdata = halodata[0]["UnitInfo"]
	cosmodata = halodata[0]["SimulationInfo"]
	cosmodata["ComovingBoxSize"] = np.round(cosmodata["Period"]*cosmodata["h_val"]/cosmodata["ScaleFactor"],1)

	#Now it has been loaded the dictionary keys needs to be updated in the halodata so they are compatible with OrbWeaver
	for i in range(opt.numsnaps):
		for field in inputfields.keys():
			halodata[i][field] = halodata[i].pop(inputfields[field])

	return atime, numhalos, halodata, tree, unitdata, cosmodata






def OutputOrbitCatalog(opt,
	orbitdata,datatypes,
	orbitForestIDStart,orbitForestIDEnd,inumForest,fileno,
	atime,cosmodata,unitdata,
	orbitforestdata):

	#Set the filename for the catalog
	filename = opt.outfilebasename + ".%i.orbweaver.preprocessed.hdf" %(fileno)

	print("Outputting data to file",filename)

	# Open up the hdf5 file
	hdffile = h5py.File(filename,"w")

	hdrgrp=hdffile.create_group("Header")

	hdrgrp.attrs["NSnaps"]=np.int32(opt.numsnaps)
	tothalos = 0 
	for i in range(opt.numsnaps):
		tothalos += len(orbitdata[i]['ID'])
	hdrgrp.attrs["Fileno"] = np.int32(fileno)
	hdrgrp.attrs["Total_number_of_halos"] = np.uint64(tothalos)
	hdrgrp.attrs["Start_orbit_forest_ID"] = np.uint64(orbitForestIDStart)
	hdrgrp.attrs["End_orbit_forest_ID"] = np.uint64(orbitForestIDEnd)
	hdrgrp.attrs["TEMPORALHALOIDVAL"] = np.uint64(opt.TEMPORALHALOIDVAL)
	hdrgrp.attrs["numRvirSearch"] = opt.numRvirSearch

	unitgrp=hdrgrp.create_group("Units")

	unitgrp.attrs["Comoving_or_Physical"] = 0
	unitgrp.attrs["Length_unit_to_kpc"] = unitdata["Length_unit_to_kpc"]
	unitgrp.attrs["Mass_unit_to_solarmass"] = unitdata["Mass_unit_to_solarmass"]
	unitgrp.attrs["Velocity_unit_to_kms"] = unitdata["Velocity_unit_to_kms"]

	cosmogrp=hdrgrp.create_group("Cosmology")

	cosmogrp.attrs["BoxSize"] = cosmodata["ComovingBoxSize"]
	cosmogrp.attrs["Hubble_param"] = cosmodata["h_val"]
	cosmogrp.attrs["Gravity"] = cosmodata["Gravity"]
	cosmogrp.attrs["Omega_Lambda"] = cosmodata["Omega_Lambda"]
	cosmogrp.attrs["Omega_m"] = cosmodata["Omega_m"]
	cosmogrp.attrs["Omega_b"] = cosmodata["Omega_b"]

	#Add the other omega's only if they exist otherwise set it to zero
	if("Omega_r" in cosmodata):
		cosmogrp.attrs["Omega_r"] = cosmodata["Omega_r"]
	else:
		cosmogrp.attrs["Omega_r"] = 0.0

	if("Omega_k" in cosmodata):
		cosmogrp.attrs["Omega_k"] = cosmodata["Omega_k"]
	else:
		cosmogrp.attrs["Omega_k"] = 0.0

	#Add the catalogue build params
	buildgrp = hdrgrp.create_group("Build_Params")
	buildparams = vars(opt)
	for field in buildparams.keys():
		buildgrp.attrs[field] = buildparams[field]


	# Put the data per snapshot
	for snap in range(opt.numsnaps):

		#Open up the snapgroup
		snapgrp=hdffile.create_group("Snap_%03d"%(snap))

		snapgrp.attrs["Snapnum"] = snap
		snapgrp.attrs["scalefactor"] = atime[snap]
		snapgrp.attrs["NHalos"] = np.int64(len(orbitdata[snap]["ID"]))

		# Output each of the datasets to the file
		for field in orbitdata[snap].keys():

			snapgrp.create_dataset(field,data=np.asarray(orbitdata[snap][field],dtype = datatypes[field]),compression='gzip', compression_opts=6)


	orbitgrp=hdrgrp.create_group("OrbitInfo")
	orbitgrp.attrs["Number_of_Orbital_Forest_IDs"] = inumForest
	orbitgrp.create_dataset('Number_of_Halos_In_Forest',data=orbitforestdata['Number_of_halos'], dtype = np.int64, compression='gzip', compression_opts=6)
	orbitgrp.create_dataset('OrbitForestID', data = np.arange(orbitForestIDStart, orbitForestIDEnd, dtype=np.int64), compression='gzip', compression_opts=6)

	hdffile.close()