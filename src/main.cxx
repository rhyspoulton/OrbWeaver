#include "orbweaver.h"

using namespace std;

int main(int argc,char **argv)
{
    int nthreads;
    int itotnumtypeofentries;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#else
    nthreads=1;
#endif


#ifdef USEOPENMP
    cout<<"orbweaver running with OpenMP. Number of openmp threads: "<<nthreads<<endl;
#endif

    // Put everything in a try - catch  for graceful exit
    try{

        // Load in the options from the config file
        Options opt;
        GetArgs(argc,argv,opt);
        cout.precision(10);

        vector<SnapData> snapdata;

        // First need to load in the catalogue data from VELOCIraptor and TreeFrog catalogues
        ReadData(opt,snapdata);

        //The options and data has been readcan now calculate the number of types of crossing entry to expect
        for(float i = opt.numRvirSearch;i>0.0;i-=opt.fracrvircross)
            opt.numtypeofcrossingentries++;

        //Calculate the total number of entries
        itotnumtypeofentries = opt.numtypeofcrossingentries;

        //This need to be 2x since we can also have negative (outward going) crossing points
        itotnumtypeofentries*=2;

        //Then add two more for the apsis points
        itotnumtypeofentries+=2;
        opt.totnumtypeofentries = itotnumtypeofentries;


        //Find the ages of the universe from the scalefactors in the snapdata
        for(int snap=0;snap<opt.numsnaps;snap++)
            snapdata[snap].uniage = GetUniverseAge(snapdata[snap].scalefactor);

        //Declare the orbitdata as a 2d vector so the interpolated halos can be added
        vector<OrbitData> orbitdata;

        ProcessOrbits(opt, snapdata,orbitdata);

        WriteOrbitData(opt,orbitdata);

    }
    catch(const exception&){

        return EXIT_CODE;

    }


    return EXIT_CODE;
}