import numpy as np 
import h5py
import os


def ReadOrbitData(filenamelist):


	#First see if the file exits
	if(os.path.isfile(filenamelist)==False):
		raise IOError("The filelist",filenamelist,"does not exist")

	filelist = open(filenamelist,"r")

	#Try and read the first line as int
	try:
		numfiles = int(filelist.readline())
	except ValueError:
		raise IOError("The first line of the filelist (which says the number of files), cannot be interpreted as a integer")

	fileno = np.zeros(numfiles,dtype=np.int32)
	numentries = np.zeros(numfiles,dtype=np.uint64)
	maxorbitIDs = np.zeros(numfiles,dtype=np.uint64)
	prevmaxorbitID = np.uint64(0)
	filenames = [""]*numfiles
	orbitdatatypes = {}
	orbitdatakeys = None

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
		fileno[i] = np.int32(hdffile.attrs["Fileno"][...])
		numentries[i] =  np.uint64(hdffile.attrs["Number_of_entries"][...])
		maxorbitIDs[i] = prevmaxorbitID
		prevmaxorbitID += np.uint64(hdffile.attrs["Max_orbitID"][...])

		#If the first file then file then find the dataset names and their datatypes
		if(i==0):
			orbitdatakeys = list(hdffile.keys())
			for key in orbitdatakeys:
				orbitdatatypes[key] = hdffile[key].dtype

		hdffile.close()

		#Add this filename to the filename list
		filenames[i] = filename

	#Now can initilize the array to contain the data
	totnumentries = np.sum(numentries)
	orbitdata = {key:np.zeros(totnumentries,dtype = orbitdatatypes[key]) for key in orbitdatakeys}

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
			orbitdata[key][startindex:endindex] = np.asarray(hdffile[key],dtype=orbitdatatypes[key])

		#Lets offset the orbitID to make it unique across all the data
		orbitdata["OrbitID"][startindex:endindex]+=maxorbitIDs[i]

		#Set the fileno that this came from
		orbitdata["Fileno"][startindex:endindex]=fileno[i]

		hdffile.close()

		ioffset+=numentries[i]

	return orbitdata
