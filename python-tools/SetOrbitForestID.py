import numpy as np 
from scipy.interpolate import interp1d
import h5py

def SetOrbitalForestID(opt,numhalos,halodata,tree,HaloID,orbitforestid,orbitdata,atime,treefields,orbitalfields,pos_tree,cosmodata):
	"""
	Sets the orbital forestID by finding any halos which come with opt.numRvirSearch x Rvir of
	the branch of interest over all the snapshots that it exists in

	"""

	#Set a dataset to mark a halo as done within this orbitID
	processedFlag = [[] for i in range(opt.numsnaps)]
	for i in range(opt.numsnaps):
		processedFlag[i] = np.zeros(numhalos[i],dtype=bool)

	#Lets first walk down extracting the haloIDs for this branch of interest
	mainHaloIDs = -1 * np.ones(opt.numsnaps,dtype=np.int64)
	mainOrbitHaloIDs = -1 * np.ones(opt.numsnaps,dtype=np.int64)
	snap = int(HaloID/opt.TEMPORALHALOIDVAL)
	index = int(HaloID%opt.TEMPORALHALOIDVAL-1)

	# Set the snapshot for which this branch first comes into existence
	mainRootTailSnap = int(tree[snap]["RootTail"][index]/opt.TEMPORALHALOIDVAL)
	mainRootHeadSnap = int(tree[snap]["RootHead"][index]/opt.TEMPORALHALOIDVAL)
	# Lets walk up this branch while it exists
	ID = tree[snap]["RootTail"][index]
	haloIndex = int(ID%opt.TEMPORALHALOIDVAL-1)
	haloSnap = mainRootTailSnap

    #store number of halos in orbital forest. Start with number of halos in main branch
    localhalocount = mainRootHeadSnap + 1 - mainRootTailSnap

	#Lets extract the full history of this halo
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):

		#Set the re-mapped tree information for this branch
		mainOrbitHaloID = np.uint64(snap * opt.TEMPORALHALOIDVAL + len(orbitdata[snap]["ID"]) + 1)
		orbitdata[snap]["ID"].append(mainOrbitHaloID)


		# #Set it orbiting forest ID and its Orbiting halo ID
		# orbitdata[snap]["OrbitForestID"].append(orbitforestid)
		orbitdata[snap]["OrbitingHaloID"].append(-1)

		#Store the mainOrbitID
		mainOrbitHaloIDs[snap] = mainOrbitHaloID
		mainHaloIDs[snap]=ID

		#Extract all the information for this halo if it exists
		if(ID!=0):
			#Store the orginal haloID
			orbitdata[snap]["origID"].append(np.uint64(ID))

			#Set a boolean if this halo is a host halo or not
			#orbitdata[snap]["hostFlag"].append(True if halodata[snap]["hostHaloID"][haloIndex]==-1 else False)
			orbitdata[snap]["hostFlag"].append(halodata[snap]["hostHaloID"][haloIndex])

			for field in orbitalfields:
				orbitdata[snap][field].append(halodata[snap][field][haloIndex])

		else: # Otherwise lets interpolate its properties

			#If the halo is interoplated set ist origID to -1
			orbitdata[snap]["origID"].append(-1)

			#Set a boolean if this halo if it a host halo or not based on the surrounding snapshots
			#orbitdata[snap]["hostFlag"].append(True if((halodata[haloSnap]["hostHaloID"][haloIndex]==-1) & (halodata[headSnap]["hostHaloID"][headIndex]==-1)) else False)
			orbitdata[snap]["hostFlag"].append(halodata[haloSnap]["hostHaloID"][haloIndex] & halodata[headSnap]["hostHaloID"][headIndex])

			for field in orbitalfields:

				#If the field is positional lets see if a periodicity correction needs to be done
				if((field=="Xc") | (field=="Yc") |(field=="Zc")):

					halopos=halodata[haloSnap][field][haloIndex]
					headpos=halodata[headSnap][field][headIndex]

					if(headpos-halopos>0.5*halodata[snap]["SimulationInfo"]["Period"]): halopos+=halodata[snap]["SimulationInfo"]["Period"]
					elif(headpos-halopos<-0.5*halodata[snap]["SimulationInfo"]["Period"]): halopos-=halodata[snap]["SimulationInfo"]["Period"]

					f = interp1d([haloSnap,headSnap],[halopos,headpos])

				else:

					#Set the interpolation data
					f = interp1d([haloSnap,headSnap],[halodata[haloSnap][field][haloIndex],halodata[headSnap][field][headIndex]])

				#Do the interpolation
				orbitdata[snap][field].append(f(snap))

		# Lets re-map the tree information
		if((snap==mainRootTailSnap) & (snap==mainRootHeadSnap)):
			orbitdata[snap]["Tail"].append(mainOrbitHaloID)
			orbitdata[snap]["Head"].append(mainOrbitHaloID)
		elif(snap==mainRootTailSnap):
			orbitdata[snap]["Tail"].append(mainOrbitHaloID)
			orbitdata[snap]["Head"].append(np.uint64((snap+1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap+1]["Head"]) + 1))
		elif(snap==mainRootHeadSnap):
			orbitdata[snap]["Tail"].append(np.uint64((snap-1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap-1]["Tail"])))
			orbitdata[snap]["Head"].append(mainOrbitHaloID)
		else:
			orbitdata[snap]["Tail"].append(np.uint64((snap-1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap-1]["Tail"])))
			orbitdata[snap]["Head"].append(np.uint64((snap+1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap+1]["Head"]) + 1))


		if(ID!=0):
			#Extract its head
			head = tree[snap]["Head"][haloIndex]
			headSnap = int(head/opt.TEMPORALHALOIDVAL)
			headIndex = int(head%opt.TEMPORALHALOIDVAL-1)

			#Mark this halo as done in the global
			halodata[snap]["doneFlag"][haloIndex]=True

			#Mark this halo as done so it is not walked again in this orbital forest ID
			processedFlag[snap][haloIndex]=True

			#Extract the head tail to check if this branch 
			headTail=tree[headSnap]["Tail"][headIndex]

			#Lets check if this halo merges in the next snapshot, if so then set its head to the current halo
			if(headTail!=ID):
				mainRootHeadSnap=haloSnap
				orbitdata[snap]["Head"][-1]=mainOrbitHaloID
				break

		#Lets check if the head is at the next snapshot
		if(headSnap==snap+1):

			#If it is at the next snapshot then move to the Descedant
			ID = head
			haloSnap = headSnap
			haloIndex = headIndex

		else: #Otherwise lets set its ID==0 showing that its properties needs to be interpolated
			ID = 0

        # likely unecessary
		# orbitdata[snap]["RootTail"].append(np.uint64(mainRootTailSnap * opt.TEMPORALHALOIDVAL +1))
		# orbitdata[snap]["RootHead"].append(np.uint64(mainRootHeadSnap * opt.TEMPORALHALOIDVAL +1))

	#List to keep track of which indexes are needed to be extracted per snapshot
	extractIndexes = [[] for i in range(opt.numsnaps)]
	# Now for the full history of the main branch lets find which halos come within Rvir of this halo
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):
		#Extract the ID of the halo at this snapshot and its index
		mainOrbitHaloID = mainOrbitHaloIDs[snap]
		index = int(mainOrbitHaloID%opt.TEMPORALHALOIDVAL-1)
		mainOrbitHaloNpart = orbitdata[snap]["npart"][index]

		#Lets find any halos that are within opt.numRvirSearch of this halo
		indexes = pos_tree[snap].query_ball_point([orbitdata[snap]["Xc"][index],orbitdata[snap]["Yc"][index],orbitdata[snap]["Zc"][index]],r = opt.numRvirSearch * orbitdata[snap]["R_200crit"][index])

		# Walk along the branches of the halos within opt.numRvirSearch
		for iindex in indexes:
			#Skip this halo if its orbital halo ID has already been set to this one or the halo is greater than the number of particles in the halo currently being orbited
			if((processedFlag[snap][iindex]) | (halodata[snap]["npart"][iindex]>mainOrbitHaloNpart)):
				continue

			# Note: no interpolation is done here for the orbiting branches as this will be done by OrbWeaver

			# Go straight to the branches  tail 
			ID = tree[snap]["Tail"][iindex]
			haloSnap = int(ID/opt.TEMPORALHALOIDVAL)
			haloIndex = int(ID%opt.TEMPORALHALOIDVAL-1)

			#Start to walk up the branch
			head = tree[haloSnap]["Head"][haloIndex]
			headSnap = int(head/opt.TEMPORALHALOIDVAL)
			headIndex = int(head%opt.TEMPORALHALOIDVAL-1)

			#Lets also store its TailSanp
			tail = tree[haloSnap]["Tail"][haloIndex]
			tailSnap = int(tail/opt.TEMPORALHALOIDVAL)

			#Check we are still on the main branch
			headTail = tree[headSnap]["Tail"][headIndex]

			# #Set the re-mapped RootTails and RootHead for this branch

			# #If this branch RootTailSnap is before the mainRootTailSnap then set this branch to only exits after mainRootTailSnap
			# if(haloSnap<mainRootTailSnap):
			# 	reMappedRootTail = np.uint64(mainRootTailSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[mainRootTailSnap]["ID"]) + 1)
			# else: #Otherwise it starts at its RootTailSnap
			# 	reMappedRootTail = np.uint64(haloSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)

			# RootHeadSnap = int(halodata[haloSnap]["RootHead"][iindex]/opt.TEMPORALHALOIDVAL)
			# #Do the same for the RootHeadSnap if it exist after mainRootHeadSnap
			# if(RootHeadSnap>mainRootHeadSnap):
			# 	reMappedRootHead = np.uint64(mainRootHeadSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[mainRootHeadSnap]["ID"]) + 1)
			# else: #Otherwise it starts at its RootTailSnap
			# 	reMappedRootHead = np.uint64(RootHeadSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[RootHeadSnap]["ID"]) + 1)

			# Now go from the base of the branch up
			while(True):

				# Only set to be part of this OrbitForestID if the main branch exists
				if(haloSnap>=mainRootTailSnap):

					#Have a list of all the extracted indexes per
					extractIndexes[haloSnap].append(haloIndex)

					#Lets also store its TailSanp
					tail = tree[haloSnap]["Tail"][haloIndex]
					tailSnap = int(tail/opt.TEMPORALHALOIDVAL)

					#Check we are still on the main branch
					headTail = tree[headSnap]["Tail"][headIndex]

					# Set this halos orbital forest ID and the halo it is orbiting
					# orbitdata[haloSnap]["OrbitForestID"].append(orbitforestid)
					orbitdata[haloSnap]["OrbitingHaloID"].append(mainOrbitHaloIDs[haloSnap])

					#Set a boolean if this halo is a host halo or not
					#orbitdata[haloSnap]["hostFlag"].append(True if halodata[haloSnap]["hostHaloID"][haloIndex]==-1 else False)
					orbitdata[haloSnap]["hostFlag"].append(halodata[haloSnap]["hostHaloID"][haloIndex])

					# Re-map its tree properties
					orbitdata[haloSnap]["origID"].append(np.uint64(ID))
					orbitingHaloID = np.uint64(haloSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)
					orbitdata[haloSnap]["ID"].append(orbitingHaloID)

					#Mark this halo as done so it is not walked again in this orbital forest ID
					# halodata[haloSnap]["OrbitForestID"][haloIndex] = orbitforestid
					processedFlag[haloSnap][haloIndex] = True

                    # increment local halo counter
                    localhalocount += 1

					# Mark the end of the branch if we have reached end of its existence or at the end of the main orbiting halo's existence
					if((head==ID) | (haloSnap==mainRootHeadSnap) | (headSnap>mainRootHeadSnap)):
						orbitdata[haloSnap]["Tail"].append(np.uint64(tailSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[tailSnap]["Tail"])))
						orbitdata[haloSnap]["Head"].append(orbitingHaloID)
						break
					elif(headTail!=ID): # If its heads tail is not pointing back to itself then it has merged with something

						orbitdata[haloSnap]["Tail"].append(np.uint64(tailSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[tailSnap]["Tail"])))

						#See if this halo merges with the halo its is orbiting if so then set it to merge with it or have it point back to itself
						sel = np.where(mainOrbitHaloIDs==head)[0]
						if(sel):
							orbitdata[haloSnap]["Head"].append(mainOrbitHaloIDs[sel])
						else:
							orbitdata[haloSnap]["Head"].append(orbitingHaloID)

						break

					elif((haloSnap==tailSnap) | (tailSnap<=mainRootTailSnap)): # If at the base of the branch or when the main orbiting halo formed
						orbitdata[haloSnap]["Tail"].append(orbitingHaloID)
						orbitdata[haloSnap]["Head"].append(np.uint64(headSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[headSnap]["Head"]) + 1))

					else: #Otherwise set its head/ tail as normal
						orbitdata[haloSnap]["Tail"].append(np.uint64(tailSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[tailSnap]["Tail"])))
						orbitdata[haloSnap]["Head"].append(np.uint64(headSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[headSnap]["Head"]) + 1))


				#Stop walking if this branch if it merges, or at the end of the simulation or the main orbit branch no longer exits
				# if((headTail!=ID) | (head==ID) | (haloSnap==mainRootHeadSnap)):
				# 	print("Got here")
				# 	break

				# Move to its head
				ID = head 
				haloSnap = headSnap
				haloIndex = headIndex

				#Extract its descendant
				head = tree[haloSnap]["Head"][haloIndex]
				headSnap = int(head/opt.TEMPORALHALOIDVAL)
				headIndex = int(head%opt.TEMPORALHALOIDVAL-1)
	#Lets set all its orbital properties
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):
		for field in orbitalfields:
			orbitdata[snap][field].extend(halodata[snap][field][extractIndexes[snap]].tolist())

    return localhalocount


def OutputOrbitalForestIDFile(opt,
    orbitdata,datatypes,
    orbitForestIDStart,orbitForestIDEnd,
    atime,cosmodata,unitdata,
    orbitforestdata):

	#Set the filename for the catalog
	filename = opt.outfilebasename + ".orbweaver.orbitForestIDs.%09d-%09d.hdf" %(orbitForestIDStart,orbitForestIDEnd) 

	print("Outputting data to file",filename)

	# Open up the hdf5 file
	hdffile = h5py.File(filename,"w")

	hdrgrp=hdffile.create_group("Header")

	hdrgrp.attrs["NSnaps"]=opt.numsnaps
    tothalos = 0 
    for i in range(opt.numsnaps):
        tothalos += orbitdata[i]['ID'].size
	hdrgrp.attrs["Total_number_of_Halos"]=np.int64(tothalos)

	unitgrp=hdrgrp.create_group("Units")

	unitgrp.attrs["Comoving_or_Physical"] = unitdata["Comoving_or_Physical"]
	unitgrp.attrs["Length_unit_to_kpc"] = unitdata["Length_unit_to_kpc"]
	unitgrp.attrs["Mass_unit_to_solarmass"] = unitdata["Mass_unit_to_solarmass"]
	unitgrp.attrs["Velocity_unit_to_kms"] = unitdata["Velocity_unit_to_kms"]

	cosmogrp=hdrgrp.create_group("Cosmology")

	cosmogrp.attrs["BoxSize"] = np.round(cosmodata["Period"]*cosmodata["h_val"]/cosmodata["ScaleFactor"],1)
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
	orbitgrp.attrs["Number_of_Orbital_Forest_IDs"]=orbitalForestIDEnd - orbitalForestIDStart
    orbitgrp.create_dataset('Number_of_Halos_In_Forest',data=orbitforestdata['Number_of_halos'], dtype = np.int64, compression='gzip', compression_opts=6)
    orbitgrp.create_dataset('OrbitForestID', data = np.arange(orbitalForestIDStart, orbitalForestIDEnd, dtype=np.int64), compression='gzip', compression_opts=6)

	hdffile.close()
