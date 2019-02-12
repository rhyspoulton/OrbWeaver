import subprocess
import argparse
import os
import sys

scriptpath=os.path.abspath(__file__)
pythontoolsdir=scriptpath.split('runorbweaver.py')[0]
sys.path.append(pythontoolsdir)
baseorbweaverdir  = pythontoolsdir.split('python-tools')[0]

#Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i",action="store",dest="inputfilelist",help="file list containing the names of the preprocessed catalogues",required=True)
parser.add_argument("-s",action="store",default='None',dest="schedulertype",help="The type of scheduler avalible either Slurm, PBS or None, if None then python's multiprocessing will be used to run orbweaver concurrently",required=True)
parser.add_argument("-c",action="store",type=float,default=0.0,dest="fracrvircross",help="Configuration file (orbweaver.cfg)",required=False)
parser.add_argument("-v",action="store",type=int,default=0,dest="iverbose",help="How verbose the code is 0=None, 1=talkative",required=False)
tmpOpt = parser.parse_args()


if(tmpOpt.fracrvircross<1e-5):
	print("The fraction of the host's rvir at which crossing points are to be outputted (-c) has not been supplied, defaulting to 0.5 * Rvir_host")
	tmpOpt.fracrvircross=0.5
elif(tmpOpt.fracrvircross<0.1):
	raise ValueError("The fraction of the host's rvir at which crossing points are to be outputted is set below 0.1, which is most likely smaller than the simulation output times.\nPlease input a value 0.1<=c<=3.0")
elif(tmpOpt.fracrvircross>3.0):
	raise ValueError("The fraction of the host's rvir at which crossing points are to be outputted is set above 3.0, which is outside where orbit properties are calculated.\nPlease input a value 0.1<=c<=3.0")

with open(tmpOpt.inputfilelist,"r") as filenamelist:

	#Try and read the first line as int
	try:
		numfiles = int(filenamelist.readline())
	except ValueError:
		raise IOError("The first line of the filelist (which says the number of files), cannot be interpreted as a integer")

	#Loop over all the filename in the list
	for i,basefilename in enumerate(filenamelist):

		#Remove any whitespace
		basefilename=basefilename.strip()

		#Check if the string has a value
		if(basefilename==""):
			continue

		#Set the name of the input file
		inputfilename = basefilename + ".orbweaver.preprocessed.hdf"

		#Check if either scheduler type avalible has been selected
		if((tmpOpt.schedulertype=="PBS") | (tmpOpt.schedulertype=="Slurm")):

			#Lets create a folder for all the submit files
			os.makedirs(pythontoolsdir + "/runscripts", exist_ok=True)

			#Get the base submit file names
			if(tmpOpt.schedulertype=="PBS"): basesubmitfilename = pythontoolsdir+"examples/qsub.runorbweaver.base.sh"
			else: basesubmitfilename = pythontoolsdir+"examples/runorbweaver.base.sbatch"

			#Open up the base submit file
			with open(basesubmitfilename,"r") as f:
				
				#Update the file contents
				fc = f.read()
				fc = fc.replace("JOBNAME","OrbWeaver-%i"%i)
				fc = fc + "\ncd " + baseorbweaverdir + "/build"
				fc = fc + "\n./orbweaver -i " + inputfilename + " -o " + basefilename + " -c " + str(tmpOpt.fracrvircross) + " -v " + str(tmpOpt.iverbose)
				
			#Open up a new submit file with all the updated contents
			if(tmpOpt.schedulertype=="PBS"): submitfilename = pythontoolsdir + "/runscripts/qsub.runorbweaver.%i.sh" %i 
			else: submitfilename = pythontoolsdir + "/runscripts/runorbweaver.%i.sbatch" %i 
			submitfile = open(submitfilename,"w")
			submitfile.write(fc)
			submitfile.close()

			#Submit the job
			if(tmpOpt.schedulertype=="PBS"): returncode = subprocess.call(["qsub",submitfilename],shell=True)
			else: returncode = subprocess.call("sbatch "+submitfilename,shell=True)

			#Check if it returned a non-zero error on submission
			if(returncode!=0):
				raise OSError("Excution of job submission failed, are you sure",tmpOpt.schedulertype,"is avalible if not please set -s None")
			else:
				print("OrbWeaver has been submitted for file",inputfilename)

		else:

			#This is where multiprocessing is to be done 
			print("Python multiprocessing is currently not implemented")