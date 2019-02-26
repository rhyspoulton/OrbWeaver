import numpy as np 
from scipy.interpolate import InterpolatedUnivariateSpline
import h5py

def LogInterp(prevdata,nextdata,f):
	inputdtype = prevdata.dtype
	prevdata = np.float64(prevdata)
	nextdata = np.float64(nextdata)
	interpoutput = (nextdata**f) * (prevdata**(1-f))
	return interpoutput.astype(inputdtype)

def LinInterp(prevdata,nextdata,f):
	inputdtype = prevdata.dtype
	prevdata = np.float64(prevdata)
	nextdata = np.float64(nextdata)
	interpoutput = prevdata + (nextdata - prevdata)*f
	return interpoutput.astype(inputdtype)


def CreateOrbitForest(opt,numhalos,halodata,tree,HaloID,orbitforestid,orbitdata,atime,treefields,orbitalfields,pos_tree,cosmodata):
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
	mainRootTailSnap = tree[snap]["RootTailSnap"][index]
	mainRootHeadSnap = tree[snap]["RootHeadSnap"][index]

	# Lets walk up this branch while it exists
	ID = tree[snap]["RootTail"][index]
	haloIndex = int(ID%opt.TEMPORALHALOIDVAL-1)
	haloSnap = mainRootTailSnap

	#Flag to keep track if the host has merged
	hostMerges = False

	#Walk this halo and see if it merges with anything
	while(True):

		headID = tree[haloSnap]["Head"][haloIndex]
		headSnap = tree[haloSnap]["HeadSnap"][haloIndex]
		headIndex = tree[haloSnap]["HeadIndex"][haloIndex]

		headTail = tree[headSnap]["Tail"][headIndex]

		if(haloSnap==mainRootHeadSnap):
			break
		elif(headTail!=ID):

			#If this halo has existed less than 20 snapshots before merging then don't include it in the orbit catalog
			if((haloSnap-mainRootTailSnap)<opt.MinSnapExist):
				return 0
			else:
				mainRootHeadSnap = haloSnap
				hostMerges = True
				break

		ID = headID
		haloIndex = headIndex
		haloSnap = headSnap

	# Lets walk up this branch while it exists
	ID = tree[snap]["RootTail"][index]
	haloIndex = int(ID%opt.TEMPORALHALOIDVAL-1)
	haloSnap = mainRootTailSnap

	#Store the snapshots to be interpolated
	interpsnaps = []
	splinefields = ["Xc","Yc","Zc","VXc","VYc","VZc","Lx","Ly","Lz"]
	loginterpfields = ["Mass_200crit","R_200crit","Vmax","Rmax","Mass_tot","Mass_FOF"]
	lininterpfields = ["npart","cNFW","numSubStruct"]

	#Store the orbital data
	tmporbitdata = {field:-1 * np.ones(opt.numsnaps,dtype=halodata[0][field].dtype) for field in splinefields}

	#store number of halos in orbital forest. Start with number of halos in main branch
	localhalocount = mainRootHeadSnap + 1 - mainRootTailSnap

	#Lets extract the full history of this halo
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):

		if(snap==mainRootTailSnap):
			#At the first snapshot set the tailID = haloID
			mainOrbitHaloID = np.int64(haloSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)
			mainOrbitTailID = mainOrbitHaloID
			if(mainRootHeadSnap>snap+1):
				mainOrbitHeadID = np.int64((snap+1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap+1]["ID"]) + 1)
		elif(snap==mainRootHeadSnap):
			#If at the final snapshot set headID==haloID
			mainOrbitHaloID = mainOrbitHeadID
		else:
			#Set the re-mapped tree information for this branch
			mainOrbitHaloID = mainOrbitHeadID
			mainOrbitHeadID = np.int64((snap+1) * opt.TEMPORALHALOIDVAL + len(orbitdata[snap+1]["ID"]) + 1)

		# Lets re-map the tree information
		orbitdata[snap]["ID"].append(mainOrbitHaloID)
		orbitdata[snap]["Tail"].append(mainOrbitTailID)
		orbitdata[snap]["Head"].append(mainOrbitHeadID)


		#Set its OrbitingHaloID to -1 as it is a host halo
		orbitdata[snap]["OrbitedHaloID"].append(-1)

		#Store the mainOrbitID
		mainOrbitHaloIDs[snap] = mainOrbitHaloID
		mainHaloIDs[snap]=ID

		#Extract all the information for this halo if it exists
		if(ID!=0):
			#Store the orginal haloID
			orbitdata[snap]["OrigID"].append(np.uint64(ID))

			#Set a boolean if this halo is a host halo or not
			orbitdata[snap]["FieldHalo"].append(halodata[snap]["hostHaloID"][haloIndex]==-1)

			#Append if the orbit host has merged
			orbitdata[snap]["hostMerges"].append(hostMerges)

			#Find the ratio of mass in substructure only if these if sub structure
			prevMassinSubStructure =halodata[snap]["MassinSubStruct"][haloIndex]/halodata[snap]["Mass_200crit"][haloIndex]

			orbitdata[snap]["RatioOfMassinSubsStruct"].append(prevMassinSubStructure)

			for field in orbitalfields:
				orbitdata[snap][field].append(halodata[snap][field][haloIndex])

				#Store the data for the spline functions
				if(field in splinefields):
					tmporbitdata[field][snap] = halodata[snap][field][haloIndex]

		else: # Otherwise lets interpolate its properties

			#Store The snpashots to be interpolated
			interpsnaps.append(snap)

			#If the halo is interpolated set ist origID to 0
			orbitdata[snap]["OrigID"].append(0)

			#Set a boolean if this halo if it a host halo or not based on the surrounding snapshots
			orbitdata[snap]["FieldHalo"].append(halodata[haloSnap]["hostHaloID"][haloIndex]==-1 & halodata[headSnap]["hostHaloID"][headIndex]==-1)

			#Append if the orbit host has merged
			orbitdata[snap]["hostMerges"].append(hostMerges)

			#Find the f needed to interpolate
			f = (snap - haloSnap)/(headSnap - haloSnap)

			#Find the ratio of mass in substructure
			# nextSubStructIndexes = np.where(halodata[snap]["hostHaloID"]==ID)[0]
			nextMassinSubStructure =halodata[headSnap]["MassinSubStruct"][headIndex]/halodata[headSnap]["Mass_200crit"][headIndex]

			orbitdata[snap]["RatioOfMassinSubsStruct"].append(LinInterp(prevMassinSubStructure,nextMassinSubStructure,f))

			#Do a log interpolation of the fields
			for field in loginterpfields:
				orbitdata[snap][field].append(LogInterp(halodata[haloSnap][field][haloIndex],halodata[headSnap][field][headIndex],f))

			#Do a linear interpolation of the fields
			for field in lininterpfields:
				orbitdata[snap][field].append(LinInterp(halodata[haloSnap][field][haloIndex],halodata[headSnap][field][headIndex],f))

		#Store the current haloID as the next tail
		mainOrbitTailID = mainOrbitHaloID

		#Check if non-zero haloID if so then extract its head information
		if(ID!=0):
			#Extract its head
			head = tree[snap]["Head"][haloIndex]
			headSnap = tree[snap]["HeadSnap"][haloIndex]
			headIndex = tree[snap]["HeadIndex"][haloIndex]

			#Mark this halo as done in the global
			halodata[snap]["doneFlag"][haloIndex]=True

			#Mark this halo as done so it is not walked again in this orbital forest ID
			processedFlag[snap][haloIndex]=True

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


	if(len(interpsnaps)):
		for field in splinefields:

			halosnaps = np.where(tmporbitdata[field]!=-1)[0]
			nhalo = halosnaps.size

			#If the field is positional lets see if a periodicity correction needs to be done
			if((field=="Xc") | (field=="Yc") |(field=="Zc")):
				boundry = 0

				for i in range(nhalo):
					tmporbitdata[field][halosnaps[i]] = tmporbitdata[field][halosnaps[i]]+boundry;
					if(i<nhalo-1):
						nextpos=tmporbitdata[field][halosnaps[i+1]]+boundry
						if(nextpos-tmporbitdata[field][halosnaps[i]]>0.5*cosmodata["ComovingBoxSize"]): boundry-=cosmodata["ComovingBoxSize"]
						elif(nextpos-tmporbitdata[field][halosnaps[i]]<-0.5*cosmodata["ComovingBoxSize"]): boundry+=cosmodata["ComovingBoxSize"]

			#Lets setup the interpolation routine
			f = InterpolatedUnivariateSpline(halosnaps,tmporbitdata[field][halosnaps])

			#Interpolate the data
			for snap in interpsnaps:
				orbitdata[snap][field].append(f(snap))

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

			# Extract the haloID
			ID = tree[snap]["ID"][iindex]
			haloSnap = snap
			haloIndex = iindex

			#Start to walk up the branch
			head = tree[haloSnap]["Head"][haloIndex]
			headSnap = tree[haloSnap]["HeadSnap"][haloIndex]
			headIndex = tree[haloSnap]["HeadIndex"][haloIndex]

			#Check we are still on the main branch
			headTail = tree[headSnap]["Tail"][headIndex]

			#Lets check if its the halo doesn't point to its own head and so exists in the next snapshot
			if((ID==head) | (headSnap>mainRootHeadSnap) | (headTail!=ID)):
				continue

			#Store the halos original ID
			orbitdata[haloSnap]["OrigID"].append(np.uint64(ID))

			#At the first snapshot lets set the tail == haloID and its head
			orbitingHaloID = np.int64(haloSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[haloSnap]["ID"]) + 1)
			orbitingTailID = orbitingHaloID
			orbitingHeadID = np.int64(headSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[headSnap]["ID"]) + 1)

			#Add to the orbit catalog
			orbitdata[haloSnap]["ID"].append(orbitingHaloID)
			orbitdata[haloSnap]["Tail"].append(orbitingHaloID)
			orbitdata[haloSnap]["Head"].append(orbitingHeadID)

			#Update a list of all the extracted indexes per
			extractIndexes[haloSnap].append(haloIndex)

			# Set this halos orbital forest ID and the halo it is orbiting
			orbitdata[haloSnap]["OrbitedHaloID"].append(mainOrbitHaloIDs[haloSnap])

			#Set a boolean if this halo is a host halo or not
			orbitdata[haloSnap]["FieldHalo"].append(halodata[haloSnap]["hostHaloID"][haloIndex]==-1)

			#Append if this halo's orbit host has merged
			orbitdata[haloSnap]["hostMerges"].append(hostMerges)

			#Find the ratio of mass in substructure only if these if sub structure
			orbitdata[haloSnap]["RatioOfMassinSubsStruct"].append(halodata[haloSnap]["MassinSubStruct"][haloIndex]/halodata[haloSnap]["Mass_200crit"][haloIndex])


			# increment local halo counter
			localhalocount += 1

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

				# Move to its head
				ID = head
				haloSnap = headSnap
				haloIndex = headIndex

				#Extract its descendant
				head = tree[haloSnap]["Head"][haloIndex]
				headSnap = tree[haloSnap]["HeadSnap"][haloIndex]
				headIndex = tree[haloSnap]["HeadIndex"][haloIndex]

				#Have a list of all the extracted indexes per
				extractIndexes[haloSnap].append(haloIndex)

				#Check we are still on the main branch
				headTail = tree[headSnap]["Tail"][headIndex]

				# Set this halos orbital forest ID and the halo it is orbiting
				orbitdata[haloSnap]["OrbitedHaloID"].append(mainOrbitHaloIDs[haloSnap])

				#Set a boolean if this halo is a host halo or not
				orbitdata[haloSnap]["FieldHalo"].append(halodata[haloSnap]["hostHaloID"][haloIndex]==-1)

				#Append if this halo's orbit host has merged
				orbitdata[haloSnap]["hostMerges"].append(hostMerges)

				# Re-map its tree properties
				orbitdata[haloSnap]["OrigID"].append(np.uint64(ID))
				orbitingHaloID = orbitingHeadID
				orbitingHeadID = np.int64(headSnap * opt.TEMPORALHALOIDVAL + len(orbitdata[headSnap]["ID"]) + 1)
				orbitdata[haloSnap]["ID"].append(orbitingHaloID)

				#Find the ratio of mass in substructure only if these if sub structure
				orbitdata[haloSnap]["RatioOfMassinSubsStruct"].append(halodata[haloSnap]["MassinSubStruct"][haloIndex]/halodata[haloSnap]["Mass_200crit"][haloIndex])

				#Mark this halo as done so it is not walked again in this orbital forest ID
				# halodata[haloSnap]["OrbitForestID"][haloIndex] = orbitforestid
				processedFlag[haloSnap][haloIndex] = True

				# increment local halo counter
				localhalocount += 1

				# Mark the end of the branch if we have reached end of its existence or at the end of the main orbiting halo's existence
				if((head==ID) | (headSnap>mainRootHeadSnap)):
					orbitdata[haloSnap]["Tail"].append(orbitingTailID)
					orbitdata[haloSnap]["Head"].append(orbitingHaloID)
					break

				elif(headTail!=ID): # If its heads tail is not pointing back to itself then it has merged with something

					orbitdata[haloSnap]["Tail"].append(orbitingTailID)

					#See if this halo merges with the halo its is orbiting if so then set it to merge with it or have it point back to itself
					sel = np.where(mainHaloIDs==head)[0]
					if(sel.size):
						orbitdata[haloSnap]["Head"].append(mainOrbitHaloIDs[sel])
					else:
						orbitdata[haloSnap]["Head"].append(orbitingHaloID)

					break

				else: #Otherwise set its head/ tail as normal
					orbitdata[haloSnap]["Tail"].append(orbitingTailID)
					orbitdata[haloSnap]["Head"].append(orbitingHeadID)

				orbitingTailID = orbitingHaloID



	#Lets set all its orbital properties
	for snap in range(mainRootTailSnap,mainRootHeadSnap+1):
		for field in orbitalfields:
			orbitdata[snap][field].extend(halodata[snap][field][extractIndexes[snap]].tolist())

	return localhalocount
