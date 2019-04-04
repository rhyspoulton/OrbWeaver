import numpy as np



class Options(object):


	def __init__(self,tmpOpt):

		self.VERSION = 0.10
		self.configfile = tmpOpt.configfile
		self.inputtree = tmpOpt.inputtree
		self.inputhalobbasename = tmpOpt.inputhalobbasename
		self.outfilebasename = tmpOpt.outfilebasename
		self.numsnaps=100
		self.numRvirSearch = 4
		self.NpartLimHost = 10000
		self.MinSnapExist = 20
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

				if(line[0]=="numRvirSearch"):
					self.numRvirSearch=np.float32(line[1])

				elif(line[0]=="NpartLimHost"):
					self.NpartLimHost=np.int64(line[1])

				elif(line[0]=="MinSnapExist"):
					self.MinSnapExist=int(line[1])

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