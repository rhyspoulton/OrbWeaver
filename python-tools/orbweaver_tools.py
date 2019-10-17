import numpy as np 
import h5py
import os
import time


def ReadOrbitData(filenamelist,iFileno=False,apsispoints=True,crossingpoints=True,endpoints=True,desiredfields=[]):

	"""

	Function to read in the data from the .orbweaver.orbitdata.hdf files

	Parameters
	----------

	filenamelist : str
		The file containing the basenames for the orbweaver catalogue (can use the file generated from the creation of the (preprocessed) orbit catalogue)

	iFileno : bool, optional
		If True, the file number where the each of the entries came from is to be outputted.

	apsispoints : bool, optional
		If False, a boolean selection is done on the data to be loaded so the apsispoints are not loaded in. It is done so it also reduces the memory used. But this means the reading takes significatly longer.

	crossingpoints : bool, optional
		If False, a boolean selection is done on the data to be loaded so the crossingpoints are not loaded in. It is done so it also reduces the memory used. But this means the reading takes significatly longer.

	endpoints : bool, optional
		If False, a boolean selection is done on the data to be loaded so the endpoints are not loaded in. It is done so it also reduces the memory used. But this means the reading takes significatly longer.

	desiredfields : list, optional
		A list containing the desired fields to put returned, please see the FieldsDescriptions.md for the list of fields availible. If not supplied then all fields are returned

	Returns
	-------

	orbitdata : dict
		Dictionary of the fields to be outputted, where each field is a ndarray .

	"""
	start =  time.time()

	#See if any of the desired datasets are false
	createselection=False
	if((apsispoints==False) | (crossingpoints==False) | (endpoints==False)):
		createselection  = True

	#First see if the file exits
	if(os.path.isfile(filenamelist)==False):
		raise IOError("The filelist",filenamelist,"does not exist")

	filelist = open(filenamelist,"r")

	#Try and read the first line as int
	try:
		numfiles = int(filelist.readline())
	except ValueError:
		raise IOError("The first line of the filelist (which says the number of files), cannot be interpreted as a integer")

	if(iFileno): fileno = np.zeros(numfiles,dtype=np.int32)
	numentries = np.zeros(numfiles,dtype=np.uint64)
	maxorbitIDs = np.zeros(numfiles,dtype=np.uint64)
	prevmaxorbitID = np.uint64(0)
	filenames = [""]*numfiles
	orbitdatatypes = {}
	orbitdatakeys = None
	if(createselection): selArray = [[] for i in range(numfiles)]

	#Loop through the filelist and read the header of each file
	for i in range(numfiles):

		#Extract the filename
		filename = filelist.readline().strip("\n")
		filename+=".orbweaver.orbitdata.hdf"

		if(os.path.isfile(filename)==False):
			raise IOError("The file",filename,"does not exits")

		#Open up the file
		hdffile = h5py.File(filename,"r")

		#Read the header information
		if(iFileno): fileno[i] = np.int32(hdffile.attrs["Fileno"][...])
		numentries[i] =  np.uint64(hdffile.attrs["Number_of_entries"][...])
		maxorbitIDs[i] = prevmaxorbitID
		prevmaxorbitID += np.uint64(hdffile.attrs["Max_orbitID"][...])

		#Use the entrytype dataset to find points to be extracted
		if(createselection):

			#Load in the entrytype dataset to create the selection
			ientrytype = np.asarray(hdffile["entrytype"],dtype=np.float64)

			#Create an array the size of the number of entries to load in the data
			sel = np.zeros(numentries[i],dtype=bool)

			#If want apsis points
			if(apsispoints):
				sel = (np.round(np.abs(ientrytype),1) == 99.0)

			#If crossing points are also desired
			if(crossingpoints):
				sel = ((np.round(np.abs(ientrytype),1) != 99.0) & (np.round(np.abs(ientrytype),1) != 0.0)) | sel

			#The final endpoint for the orbiting halo
			if(endpoints):
				sel = (np.round(np.abs(ientrytype),1) == 0.0) | sel

			selArray[i] = sel

			#Update the number of entries based on the selection
			numentries[i] = np.sum(sel,dtype = np.uint64)

		#If the first file then file then find the dataset names and their datatypes
		if(i==0):
			if(len(desiredfields)>0):
				orbitdatakeys = desiredfields
			else:
				orbitdatakeys = list(hdffile.keys())

			for key in orbitdatakeys:
				orbitdatatypes[key] = hdffile[key].dtype

		hdffile.close()

		#Add this filename to the filename list
		filenames[i] = filename

	#Now can initilize the array to contain the data
	totnumentries = np.sum(numentries)
	orbitdata = {key:np.zeros(totnumentries,dtype = orbitdatatypes[key]) for key in orbitdatakeys}

	if(iFileno):
		#Add a field to contain the file number
		orbitdata["Fileno"] = np.zeros(totnumentries,dtype=np.int32)

	ioffset = np.uint64(0)

	#Now read in the data
	for i in range(numfiles):

		filename=filenames[i]

		print("Reading orbitdata from",filename)

		#Open up the file
		hdffile = h5py.File(filename,"r")

		#Get the start and end index
		startindex = ioffset
		endindex = np.uint64(ioffset+numentries[i])

		#Read the datasets
		for key in orbitdatakeys:
			if(createselection):
				orbitdata[key][startindex:endindex] = np.asarray(hdffile[key][selArray[i]],dtype=orbitdatatypes[key])
			else:
				orbitdata[key][startindex:endindex] = np.asarray(hdffile[key],dtype=orbitdatatypes[key])

		if("OrbitID" in orbitdatakeys):
			#Lets offset the orbitID to make it unique across all the data
			orbitdata["OrbitID"][startindex:endindex]+= np.uint64(maxorbitIDs[i] + i)

		if(iFileno):
			#Set the fileno that this came from
			orbitdata["Fileno"][startindex:endindex]=fileno[i]

		hdffile.close()

		ioffset+=numentries[i]

	print("Done reading in",time.time()-start)

	return orbitdata
