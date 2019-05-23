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

				elif(line[0]=="HostSelectionFile"):
					self.HostSelectionFile = line[1]

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



def GetHostSelection(opt,fields):

	snapshotselexpr = ""
	hostselexpr  = ""
	operators = "=><!()"

	with open(opt.HostSelectionFile,"r") as file:

		for i,line in enumerate(file):

			if(line[0]=="#"): continue

			line = line.replace(" ","")

			line = line.strip()

			if(not line): continue

			ifound = False

			#See if there has been a desired range to find hosts
			if(line.find("redshift")!=-1):

				#Replace inf with np.inf
				line = line.replace("inf","np.inf")

				#Lets test the line is valid
				redshift = 0
				try:
					eval(line)
				except SyntaxError:
					raise SyntaxError("Please check line %i (%s) for the redshift range in %s is valid with python syntax" %(i,line,opt.HostSelectionFile))

				if(snapshotselexpr==""):
					snapshotselexpr += "(" +line + ")"
				else:
					snapshotselexpr += " & " + "(" +line + ")"

			elif(line.find("MinSnapExist")!=-1): #See if there is a desired number of snapshost that the host exists for

				line = line.split("=")

				opt.MinSnapExist = int(line[1])

			else:
				#Store where the fieldname starts and the fieldname
				fieldstarts = []
				fieldnames = []

				#Lets see if this dataset is in the halodata
				for field in fields:

					#The offset if more than one occurance if the string is found
					offset = 0

					#Loop over number of occurances of the string
					for inum in range(line.count(field)):

						indx = line.find(field,offset)

						#See if the string is surrounded by at least one operator(s)
						if(indx==0):
							if(line[indx+len(field)] in operators):
								fieldstarts.append(indx)
								fieldnames.append(field)

						elif(indx+len(field)==len(line)):
							if(line[indx-1] in operators):
								fieldstarts.append(indx)
								fieldnames.append(field)

						elif(line[indx-1] in operators)  & (line[indx+len(field)] in operators):
								fieldstarts.append(indx)
								fieldnames.append(field)

						#Lets make sure this has not been already been set
						if(field not in locals()):
							locals()[field] = 0

						offset = indx + len(field)

				#Lets test that the line is valid
				try:
					eval(line)
				except SyntaxError:
					raise SyntaxError("Please check line %i (%s) in %s is valid with python syntax" %(i,line,opt.HostSelectionFile))

				#Check that a valid datasets have been found
				if(len(fieldstarts)==0):
					raise IOError("Issue on line %i (%s) in the configuration file please check the dataset exists \ncheck: https://github.com/ICRAR/VELOCIraptor-STF/blob/master/doc/output.rst" % (i,line))


				#Lets set the line to point to the correct dataset int the halodata
				for indx,field in zip(fieldstarts,fieldnames):
					tmp = [line[:indx],line[indx+len(field):]]
					line = "halodata[j]['" + field + "']"
					line = line.join(tmp)

				#Add in the space around the or operator
				line = line.split("|")
				line = " | ".join(line)

				#Add this into the host selection string
				if(hostselexpr==""):
					hostselexpr += "(" +line + ")"
				else:
					hostselexpr += " & " + "(" +line + ")"

	return hostselexpr, snapshotselexpr
