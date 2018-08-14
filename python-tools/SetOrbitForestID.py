import numpy as np 
from scipy.interpolate import interp1d
import h5py


def SetOrbitalForestID(snap,numsnaps,halodata,HaloID,orbitforestid,orbitdata,treefields,orbitalfields,pos_tree,
	TEMPORALHALOIDVAL = 1000000000000,searchSnapLim=5,numRvirSearch=4,ireversesnaporder=False):
	"""
	Sets the orbital forestID by finding any halos which come with numRvirSearch x Rvir of
	the branch of interest over all the snapshots that it exists in

	"""

	#Lets first walk down extracting the haloIDs for this branch of interest
	mainOrbitHaloIDs = -1 * np.ones(numsnaps,dtype=np.int64)
	snap = int(HaloID/TEMPORALHALOIDVAL)
	index = int(HaloID%TEMPORALHALOIDVAL-1)

	# Set the snapshot for which this branch first comes into existence
	mainRootTailSnap = int(halodata[snap]["RootTail"][index]/TEMPORALHALOIDVAL)
	mainRootHeadSnap = int(halodata[snap]["RootHead"][index]/TEMPORALHALOIDVAL)
	# Lets walk up this branch while it exists
	ID = halodata[snap]["RootTail"][index]
	haloIndex = int(ID%TEMPORALHALOIDVAL-1)
	haloSnap = mainRootTailSnap

	#Lets extract the full history of this halo
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):


		#Extract all the information for this halo if it exists
		if(ID!=0):
			# halodata[snap]["OrbitForestID"][haloIndex]=orbitforestid
			for field in orbitalfields:
				orbitdata[snap][field].append(halodata[snap][field][haloIndex])

		else: # Otherwise lets interpolate its properties

			for field in orbitalfields:

				#Set the interpolation data
				f = interp1d([haloSnap,headSnap],[halodata[haloSnap][field][haloIndex],halodata[headSnap][field][headIndex]])

				#Do the interpolation
				orbitdata[snap][field].append(f(snap))

		#Set the re-mapped tree information for this branch
		orbitdata[snap]["origID"].append(np.uint64(ID))
		mainOrbitHaloID = np.uint64(snap * TEMPORALHALOIDVAL + len(orbitdata[snap]["ID"]) + 1)
		orbitdata[snap]["ID"].append(mainOrbitHaloID)

		# #Set it orbiting forest ID and its Orbiting halo ID
		# orbitdata[snap]["OrbitForestID"].append(orbitforestid)
		orbitdata[snap]["OrbitingHaloID"].append(-1)

		halodata[snap]["OrbitingHaloID"][haloIndex]=mainOrbitHaloID


		#Store the mainOrbitID
		mainOrbitHaloIDs[snap] = mainOrbitHaloID

		# Lets re-map the tree information
		if((snap==mainRootTailSnap) & (snap==mainRootHeadSnap)):
			orbitdata[snap]["Tail"].append(mainOrbitHaloID)
			orbitdata[snap]["Head"].append(mainOrbitHaloID)
		elif(snap==mainRootTailSnap):
			orbitdata[snap]["Tail"].append(mainOrbitHaloID)
			orbitdata[snap]["Head"].append(np.uint64((snap+1) * TEMPORALHALOIDVAL + len(orbitdata[snap+1]["Head"]) + 1))
		elif(snap==mainRootHeadSnap):
			orbitdata[snap]["Tail"].append(np.uint64((snap-1) * TEMPORALHALOIDVAL + len(orbitdata[snap-1]["Tail"])))
			orbitdata[snap]["Head"].append(mainOrbitHaloID)
		else:
			orbitdata[snap]["Tail"].append(np.uint64((snap-1) * TEMPORALHALOIDVAL + len(orbitdata[snap-1]["Tail"])))
			orbitdata[snap]["Head"].append(np.uint64((snap+1) * TEMPORALHALOIDVAL + len(orbitdata[snap+1]["Head"]) + 1))

		# orbitdata[snap]["RootTail"].append(np.uint64(mainRootTailSnap * TEMPORALHALOIDVAL +1))
		# orbitdata[snap]["RootHead"].append(np.uint64(mainRootHeadSnap * TEMPORALHALOIDVAL +1))

		if(ID!=0):
			#Extract its head
			head = halodata[snap]["Head"][haloIndex]
			headSnap = int(head/TEMPORALHALOIDVAL)
			headIndex = int(head%TEMPORALHALOIDVAL-1)

		#Lets check if the head is at the next snapshot
		if(headSnap==snap+1):

			#If it is at the next snapshot then move to the Descedant
			ID = head
			haloSnap = headSnap
			haloIndex = headIndex

		else: #Otherwise lets set its ID==0 showing that its properties needs to be interpolated
			ID = 0


	#List to keep track of which indexes are needed to be extracted per snapshot
	extractIndexes = [[] for i in range(numsnaps)]
	# Now for the full history of the main branch lets find which halos come within Rvir of this halo
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):
		#Extract the ID of the halo at this snapshot and its index
		mainOrbitHaloID = mainOrbitHaloIDs[snap]
		index = int(mainOrbitHaloID%TEMPORALHALOIDVAL-1)
		mainOrbitHaloNpart = orbitdata[snap]["npart"][index]

		#Lets find any halos that are within numRvirSearch of this halo
		indexes = pos_tree[snap].query_ball_point([orbitdata[snap]["Xc"][index],orbitdata[snap]["Yc"][index],orbitdata[snap]["Zc"][index]],r = numRvirSearch * orbitdata[snap]["R_200crit"][index])


		# Walk along the branches of the halos within numRvirSearch
		for iindex in indexes:
			#Skip this halo if its orbital halo ID has already been set to this one or the halo is greater than the number of particles in the halo currently being orbited
			if((halodata[snap]["OrbitingHaloID"][iindex]==mainOrbitHaloID) | (halodata[snap]["npart"][iindex]>mainOrbitHaloNpart)):
				continue

			# Note: no interpolation is done here for the orbiting branches as this will be done by OrbWeaver

			# Go straight to the branches roottail 
			ID = halodata[snap]["RootTail"][iindex]
			haloSnap = int(ID/TEMPORALHALOIDVAL)
			haloIndex = int(ID%TEMPORALHALOIDVAL-1)

			#Start to walk up the branch
			head = halodata[haloSnap]["Head"][haloIndex]
			headSnap = int(head/TEMPORALHALOIDVAL)
			headIndex = int(head%TEMPORALHALOIDVAL-1)

			# #Set the re-mapped RootTails and RootHead for this branch

			# #If this branch RootTailSnap is before the mainRootTailSnap then set this branch to only exits after mainRootTailSnap
			# if(haloSnap<mainRootTailSnap):
			# 	reMappedRootTail = np.uint64(mainRootTailSnap * TEMPORALHALOIDVAL + len(orbitdata[mainRootTailSnap]["ID"]) + 1)
			# else: #Otherwise it starts at its RootTailSnap
			# 	reMappedRootTail = np.uint64(haloSnap * TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)

			# RootHeadSnap = int(halodata[haloSnap]["RootHead"][iindex]/TEMPORALHALOIDVAL)
			# #Do the same for the RootHeadSnap if it exist after mainRootHeadSnap
			# if(RootHeadSnap>mainRootHeadSnap):
			# 	reMappedRootHead = np.uint64(mainRootHeadSnap * TEMPORALHALOIDVAL + len(orbitdata[mainRootHeadSnap]["ID"]) + 1)
			# else: #Otherwise it starts at its RootTailSnap
			# 	reMappedRootHead = np.uint64(RootHeadSnap * TEMPORALHALOIDVAL + len(orbitdata[RootHeadSnap]["ID"]) + 1)

			# Now go from the base of the branch up
			while(True):

				# Only set to be part of this OrbitForestID if the main branch exists
				if(haloSnap>=mainRootTailSnap):

					#Lets also store its TailSanp
					tail = halodata[haloSnap]["Tail"][haloIndex]
					tailSnap = int(tail/TEMPORALHALOIDVAL)

					#Check we are still on the main branch
					headTail = halodata[headSnap]["Tail"][headIndex]


					# halodata[haloSnap]["OrbitForestID"][haloIndex] = orbitforestid
					halodata[haloSnap]["OrbitingHaloID"][haloIndex] = mainOrbitHaloIDs[haloSnap]
					#Lets set all its orbital properties
					# for field in orbitalfields:
					# 	orbitdata[haloSnap][field].append(halodata[haloSnap][field][haloIndex])

					extractIndexes[haloSnap].append(haloIndex)

					# Set this halos orbital forest ID and the halo it is orbiting
					# orbitdata[haloSnap]["OrbitForestID"].append(orbitforestid)
					orbitdata[haloSnap]["OrbitingHaloID"].append(mainOrbitHaloIDs[haloSnap])

					# Re-map its tree properties
					orbitdata[haloSnap]["origID"].append(np.uint64(ID))
					orbitingHaloID = np.uint64(haloSnap * TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)
					orbitdata[haloSnap]["ID"].append(orbitingHaloID)

					# Mark the end of the branch if it has merged with something else or end of its existence or at the end of the main orbiting halo's existence
					if((headTail!=ID) | (head==ID) | (haloSnap==mainRootHeadSnap)):
						orbitdata[haloSnap]["Tail"].append(np.uint64(tailSnap * TEMPORALHALOIDVAL + len(orbitdata[tailSnap]["Tail"])))
						orbitdata[haloSnap]["Head"].append(orbitingHaloID)
						break

					elif((haloSnap==tailSnap) | (tailSnap<=mainRootTailSnap)): # If at the base of the branch or when the main orbiting halo formed
						orbitdata[haloSnap]["Tail"].append(orbitingHaloID)
						orbitdata[haloSnap]["Head"].append(np.uint64(headSnap * TEMPORALHALOIDVAL + len(orbitdata[headSnap]["Head"]) + 1))

					else: #Otherwise set its head/ tail as normal
						orbitdata[haloSnap]["Tail"].append(np.uint64(tailSnap * TEMPORALHALOIDVAL + len(orbitdata[tailSnap]["Tail"])))
						orbitdata[haloSnap]["Head"].append(np.uint64(headSnap * TEMPORALHALOIDVAL + len(orbitdata[headSnap]["Head"]) + 1))


				#Stop walking if this branch if it merges, or at the end of the simulation or the main orbit branch no longer exits
				# if((headTail!=ID) | (head==ID) | (haloSnap==mainRootHeadSnap)):
				# 	print("Got here")
				# 	break

				# Move to its head
				ID = head 
				haloSnap = headSnap
				haloIndex = headIndex

				#Extract its descendant
				head = halodata[haloSnap]["Head"][haloIndex]
				headSnap = int(head/TEMPORALHALOIDVAL)
				headIndex = int(head%TEMPORALHALOIDVAL-1)

	#Lets set all its orbital properties
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):
		for field in orbitalfields:
			orbitdata[snap][field].extend(halodata[snap][field][extractIndexes[snap]].tolist())

def OutputOrbitalForestIDFile(numsnaps,basefilename,orbitdata,atime,orbitForestIDStart,orbitForestIDEnd,cosmodata,unitdata):

	#Set the filename for the catalog
	filename = basefilename + ".orbweaver.orbitForestIDs.%05d-%05d.hdf" %(orbitForestIDStart,orbitForestIDEnd) 

	# Open up the hdf5 file
	hdffile = h5py.File(filename,"w")

	hdrgrp=hdffile.create_group("Header")

	unitgrp=hdrgrp.create_group("Units")

	unitgrp.attrs["Comoving_or_Physical"] = False
	unitgrp.attrs["Length_unit_to_kpc"] = 1000
	unitgrp.attrs["Mass_unit_to_solarmass"] = 1.0
	unitgrp.attrs["Velocity_unit_to_kms"] = 1.0

	cosmogrp=hdrgrp.create_group("Cosmology")

	cosmogrp.attrs["BoxSize"] = cosmodata["BoxSize"]
	cosmogrp.attrs["Hubble_param"] = cosmodata["Hubble_param"]

	# Put the data per snapshot
	for snap in range(numsnaps):

		#Open up the snapgroup
		snapgrp=hdffile.create_group("Snap_%03d"%(snap))

		snapgrp.attrs["scalefactor"] = atime[snap]
		snapgrp.attrs["NHalos"] = np.int64(len(orbitdata[snap]["ID"]))

		# Output each of the datasets to the file
		for field in orbitdata[snap].keys():

			snapgrp.create_dataset(field,data=orbitdata[snap][field],compression='gzip', compression_opts=6)


	hdffile.close()





