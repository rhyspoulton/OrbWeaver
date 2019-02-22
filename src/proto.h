
#ifndef ORBPROTO_H
#define ORBPROTO_H

#include "allvars.h"

bool FileExists(const char *fname);
void GetArgs(int argc, char *argv[], Options &opt);
void ConfigCheck(Options &opt);
void ReadHeader(Options &opt, H5File *Fhdf,HDFCatalogNames hdfnames);
void ReadData(Options &opt, vector<SnapData> &snapdata);
vector<HaloData> ReadSnapshotData(int snap, int i,Group snapgroup, Options &opt, vector<SnapData> &snapdata, HDFCatalogNames hdfnames);

double GetUniverseAge(double scalefactor);
double GetScaleFactor(double uniage);

double LogInterp(double prevdata, double nextdata, double f);
double LinInterp(double prevdata, double nextdata, double f);
void SetupPosVelInterpFunctions(vector<int> &halosnaps, vector<unsigned long long> &haloindexes, vector<SnapData> &snapdata, SplineFuncs& splinefuncs);
void SetupPosVelInterpFunctionsHost(vector<int> &halosnaps, vector<unsigned long long> &hostindexes, vector<unsigned long long> &haloindexes, vector<SnapData> &snapdata, SplineFuncs &splinefuncs);
void InterpSingleHaloProps(double interpuniage, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata, SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs);
HaloData InterpHaloProps(Options &opt, vector<int> &halosnaps, vector<unsigned long long> &haloindexes, vector<int> &interpsnaps, vector<SnapData> &snapdata, SplineFuncs &splinefuncs);
double InterpCrossingHaloProps(float numrvircrossing, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata, SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs);

void CleanOrbits(vector<OrbitData> &branchorbitdata, double simtime);

void ProcessOrbits(Options &opt, vector<SnapData> &snapdata, vector<OrbitData> &orbitdata);
void ProcessHalo(Options &opt, int snap, unsigned long long i, vector<SnapData> &snapdata, OrbitData &orbitdata);
void CalcOrbitProps(Options &opt, unsigned long long orbitID, int currentsnap, int prevsnap, unsigned long long descendantProgenID, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &branchorbitdata, OrbitData &tmporbitdata, vector<SnapData> &snapdata, int* num_entrytypes, OrbitProps &orbitprops, SplineFuncs &splinefuncs, SplineFuncs &hostsplinefuncs);

void WriteOrbitData(Options &opt, vector<OrbitData> &orbitdata);
#endif //ifndef ORBPROTO_H