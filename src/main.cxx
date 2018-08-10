#include "orbweaver.h"

using namespace std;
using namespace Math;
using namespace NBody;

int main(int argc,char **argv)
{
    int nthreads;
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

        SnapData *snapdata;
        snapdata = new SnapData[opt.numsnaps];

        // First need to load in the catalogue data from VELOCIraptor and TreeFrog catalogues
        ReadData(opt,snapdata);

        //Find the ages of the universe from the scalefactors in the snapdata
        for(Int_t snap=opt.isnap;snap<=opt.fsnap;snap++)
            snapdata[snap].uniage = GetUniverseAge(snapdata[snap].scalefactor);

        //Declare the orbitdata as a 2d vector so the interpolated halos can be added
        vector<OrbitData> orbitdata;

        ProcessOrbits(opt, snapdata,orbitdata);

    }
    catch(const exception&){

        return EXIT_CODE;

    }
    




    return EXIT_CODE;
}