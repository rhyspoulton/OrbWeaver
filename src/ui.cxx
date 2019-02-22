
#include "orbweaver.h"



// void GetParamFile(Options &opt)
// {
// 	string line,sep="=";
// 	string tag,val;
// 	char buff[1024],*pbuff,tbuff[1024],vbuff[1024],fname[1024];
// 	fstream paramfile,cfgfile;
// 	if (!FileExists(opt.configname)){
// 			cerr<<"Config file does not exist or can't be read, terminating"<<endl;
// 			EXIT_CODE = 1;
// 			throw exception();
// 	}
// 	paramfile.open(opt.configname, ios::in);
// 	unsigned j,k;
// 	if (paramfile.is_open())
// 	{
// 		while (paramfile.good()){
// 			getline(paramfile,line);
// 			//if line is not commented out or empty
// 			if (line[0]!='#'&&line.length()!=0) {
// 				if (j=line.find(sep)){
// 					//clean up string
// 					tag=line.substr(0,j);
// 					strcpy(buff, tag.c_str());
// 					pbuff=strtok(buff," ");
// 					strcpy(tbuff, pbuff);
// 					val=line.substr(j+1,line.length()-(j+1));
// 					strcpy(buff, val.c_str());
// 					pbuff=strtok(buff," ");
// 					if (pbuff==NULL) continue;
// 					strcpy(vbuff, pbuff);
// 					//number of snapshot

// 					if (strcmp(tbuff, "Temporal_haloidval")==0) {
// 						opt.TEMPORALHALOIDVAL = atol(vbuff);
// 					}

// 					else if (strcmp(tbuff, "iverbose")==0) {
// 						opt.iverbose = atoi(vbuff);
// 					}

// 				}
// 			}
// 		}
// 		paramfile.close();
// 	}
// 	return;
// }



void GetArgs(int argc, char *argv[], Options &opt){

	int option;
	int itotnumtypeofentries;
	// int configflag=0;
	while ((option = getopt(argc, argv, "i:c:o:v:")) != EOF)
	{

		switch(option){
			case 'i':
				opt.fname = optarg;
				break;

			case 'c':
				opt.fracrvircross = atof(optarg);
				break;

			case 'o':
				opt.outputbasename = optarg;
				break;

			case 'v':
				opt.iverbose = atoi(optarg);
				break;

		}

	}

	// if(configflag){
	// 	cout<<"Reading config file"<<endl;
	// 	GetParamFile(opt);
	// }
	// else {
	// 	cout<<"NO CONFIG FILE PASSED! Using default values"<<endl;
	// }

	ConfigCheck(opt);

	//The config check has passed can now calculate the number of types of crossing entry to expect
	for(float i = 3.0;i>0.0;i-=opt.fracrvircross)
		opt.numtypeofcrossingentries++;

	//Calculate the total number of entries
	itotnumtypeofentries = opt.numtypeofcrossingentries;

	//This need to be 2x since we can also have negative (outward going) crossing points
	itotnumtypeofentries*=2;

	//Then add two more for the apsis points
	itotnumtypeofentries+=2;

	opt.totnumtypeofentries = itotnumtypeofentries;

	return;
}

void ConfigCheck(Options &opt){

	if(opt.fname==NULL){

		cerr<<"Must provide input file, usage: \n	orbweaver -i <input catalogue> -o <base output name>  [-c <fraction host rvir crossing point rvir 0.1<=fracrvircross<=3.0> -v <verbose 0 = none, 1 = talkative>]"<<endl;

		EXIT_CODE = 1;
		throw exception();

	}

	if(opt.fracrvircross<1e-5){

		cout<<"The fraction of the host's rvir at which crossing points are to be outputted has not been supplied, defaulting to 0.5 * Rvir_host"<<endl;
		opt.fracrvircross=0.5;
	}

	if(opt.fracrvircross<0.1){
		cerr<<"The fraction of the host's rvir at which crossing points are to be outputted is set below 0.1, which is most likely smaller than the simulation outputs.\nPlease input a value 0.1<=c<=3.0"<<endl;
		EXIT_CODE = 1;
		throw exception();
	}

	if(opt.fracrvircross>3.0){
		cerr<<"The fraction of the host's rvir at which crossing points are to be outputted is set above 3.0, which is outside where orbit properties are calculated.\nPlease input a value 0.1<=c<=3.0"<<endl;
		EXIT_CODE = 1;
		throw exception();
	}

	if(opt.outputbasename==NULL){

		cerr<<"Must provide output base filename, usage: \n	orbweaver -i <input catalogue> -o <base output name> [-c <fraction host rvir crossing point rvir 0.1<=fracrvircross<=3.0> -v <verbose 0 = none, 1 = talkative>]"<<endl;

		EXIT_CODE = 1;
		throw exception();

	}
	return;
}
