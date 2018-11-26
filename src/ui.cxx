
#include "orbweaver.h"



void GetParamFile(Options &opt)
{
	string line,sep="=";
	string tag,val;
	char buff[1024],*pbuff,tbuff[1024],vbuff[1024],fname[1024];
	fstream paramfile,cfgfile;
	if (!FileExists(opt.configname)){
			cerr<<"Config file does not exist or can't be read, terminating"<<endl;
			EXIT_CODE = 1;
			throw exception();
	}
	paramfile.open(opt.configname, ios::in);
	unsigned j,k;
	if (paramfile.is_open())
	{
		while (paramfile.good()){
			getline(paramfile,line);
			//if line is not commented out or empty
			if (line[0]!='#'&&line.length()!=0) {
				if (j=line.find(sep)){
					//clean up string
					tag=line.substr(0,j);
					strcpy(buff, tag.c_str());
					pbuff=strtok(buff," ");
					strcpy(tbuff, pbuff);
					val=line.substr(j+1,line.length()-(j+1));
					strcpy(buff, val.c_str());
					pbuff=strtok(buff," ");
					if (pbuff==NULL) continue;
					strcpy(vbuff, pbuff);
					//number of snapshot
					if (strcmp(tbuff, "isnap")==0) {
						opt.isnap = atoi(vbuff);
					}
					else if (strcmp(tbuff, "fsnap")==0) {
						opt.fsnap = atoi(vbuff);
					}

					else if (strcmp(tbuff, "Temporal_haloidval")==0) {
						opt.TEMPORALHALOIDVAL = atol(vbuff);
					}

					else if (strcmp(tbuff, "iverbose")==0) {
						opt.iverbose = atoi(vbuff);
					}

				}
			}
		}
		paramfile.close();
	}
	return;
}



void GetArgs(int argc, char *argv[], Options &opt){

	int option;
	int NumArgs = 0;
	int configflag=0;
	while ((option = getopt(argc, argv, "i:c:o:")) != EOF)
	{

		switch(option){
			case 'i':
				opt.fname = optarg;
				NumArgs+=2;
				break;

			case 'c':
				opt.configname = optarg;
				NumArgs+=2;
				configflag = 1;
				break;
			case 'o':
				opt.outputbasename = optarg;
				NumArgs+=2;
				break;

		}

	}

	if(configflag){
		cout<<"Reading config file"<<endl;
		GetParamFile(opt);
	}
	else {
		cout<<"NO CONFIG FILE PASSED! Using default values"<<endl;
	}

	ConfigCheck(opt);


	return;
}

void ConfigCheck(Options &opt){

	if(opt.fname==NULL){

		cerr<<"Must provide input file, usage: \n	orbweaver -i <input catalogue> -o <base output name> [-C <config file>]"<<endl;

		EXIT_CODE = 1;
		throw exception();

	}

	if(opt.outputbasename==NULL){

		cerr<<"Must provide output base filename, usage: \n	orbweaver -i <input catalogue> -o <base output name> [-C <config file>]"<<endl;

		EXIT_CODE = 1;
		throw exception();

	}


	if(opt.isnap>opt.fsnap){

		cerr<<"isnap > fsnap please adjust these values"<<endl;

		EXIT_CODE = 1;
		throw exception();
	}

	opt.numsnaps=opt.fsnap-opt.isnap+1;

	if(opt.numsnaps<2){

		cerr<<"Need 2 or more snapshots to process orbits please adjust the isnap and fsnap value"<<endl;

		EXIT_CODE = 1;
		throw exception();
	}
	return;
}
