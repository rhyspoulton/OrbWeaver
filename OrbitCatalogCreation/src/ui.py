import numpy as np

def GetDatasetNames(opt,basecodedir):

	orbitfields = {}

	#Find the name of the datasetfile given the desired InputFormat
	if(opt.InputFormat==0):

		datasetfile = basecodedir+"/example_inputs/input_VELOCIraptor_catalog.txt"

		print("Loading in the halos from VELOCIraptor and merger tree from TreeFrog")
		print("Extracting the datasets names from",datasetfile)

	else:
		raise NotImplementedError(opt.InputFormat,"is currently not supported. If you are interested in having OrbWeaver support it please contact the developer")

	#Extract the dataset names from the file
	with open(datasetfile,"r") as f:

		for line in f:

			if(line[0]=="#"): continue

			line = line.replace(" ","")

			line = line.strip()

			if(not line): continue

			line = line.split(":")

			if(len(line)!=2):
				IOError(line,"has no corresponding halo catalog field  please update this")

			orbitfields[line[0]] = line[1]

	return orbitfields

class Options(object):


	def __init__(self,tmpOpt):

		self.VERSION = 0.10
		self.configfile = tmpOpt.configfile
		self.inputtreefilename = tmpOpt.inputtreefilename
		self.inputhalofilelistname = tmpOpt.inputhalofilelistname
		self.outfilebasename = tmpOpt.outfilebasename
		self.numsnaps=100
		self.InputFormat = 0
		self.numRvirSearch = 4
		self.NpartLimHost = 10000
		self.MinNumSnapExistHost = 20
		self.MinNumSnapExistSat = 10
		self.TEMPORALHALOIDVAL = 1000000000000
		self.numOrbitForestPerFile = 2000
		self.iverbose = 1

		with open(tmpOpt.configfile,"r") as f:

			for line in f:

				if(line[0]=="#"): continue

				line = line.replace(" ","")

				line = line.strip()

				if(not line): continue

				line = line.split("=")

				if(line[0]=="InputFormat"):
					self.InputFormat=int(line[1])

				elif(line[0]=="numRvirSearch"):
					self.numRvirSearch=np.float32(line[1])

				elif(line[0]=="NpartLimHost"):
					self.NpartLimHost=np.int64(line[1])

				elif(line[0]=="MinNumSnapExistHost"):
					self.MinNumSnapExistHost=int(line[1])

				elif(line[0]=="MinNumSnapExistSat"):
					self.MinNumSnapExistSat=int(line[1])

				elif(line[0]=="TEMPORALHALOIDVAL"):
					self.TEMPORALHALOIDVAL=np.uint64(line[1])

				elif(line[0]=="numOrbitForestPerFile"):
					self.numOrbitForestPerFile=np.int64(line[1])

				elif(line[0]=="iverbose"):
					self.iverbose=int(line[1])

				else:
					print("Invalid config option %s, please only use the options in the sample config file" %line[0])

	def update_numsnaps(self,numsnaps):
		self.numsnaps = numsnaps