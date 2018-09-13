
#ifndef ORBPROTO_H
#define ORBPROTO_H

#include "allvars.h"

bool FileExists(const char *fname);
void GetArgs(int argc, char *argv[], Options &opt);
void ConfigCheck(Options &opt);
void ReadHeader(H5File *Fhdf,HDFCatalogNames hdfnames);
void ReadData(Options &opt, vector<SnapData> &snapdata);
vector<HaloData> ReadSnapshotData(int snap, int i,Group snapgroup, Options &opt, vector<SnapData> &snapdata, HDFCatalogNames hdfnames);

double GetUniverseAge(double scalefactor);

double LogInterp(double prevdata, double nextdata, double f);
double LinInterp(double prevdata, double nextdata, double f);
void InterpHaloPosVel(int nhalo, int ninterp, double *halouniages, double *interpuniages, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<SnapData> &snapdata, vector<HaloData> &interphalos);
void InterpPassageHaloProps(double interpuniage, double currentuniage, double prevuniage, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, OrbitData &tmporbitdata, vector<SnapData> &snapdata);
HaloData InterpHaloProps(Options &opt, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, vector<Int_t> &interpsnaps, vector<SnapData> &snapdata);
void InterpPassagePoints(vector<Int_t> halosnaps,vector<Int_t> haloindexes,vector<Int_t> hostindexes, vector<SnapData> &snapdata, vector<OrbitData> &branchorbitdata, OrbitProps &orbitprops);

void ProcessOrbits(Options &opt, vector<SnapData> &snapdata, vector<OrbitData> &orbitdata);
void ProcessHalo(Int_t snap, Int_t i, Options &opt, vector<SnapData> &snapdata, OrbitData &orbitdata);
void CalcOrbitProps(Int_t orbitID, int currentsnap, int prevsnap, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &orbitdata, OrbitData &tmporbitdata, vector<SnapData> &snapdata, OrbitProps &orbitprops);

void WriteOrbitData(Options &opt, vector<OrbitData> &orbitdata);
#endif //ifndef ORBPROTO_H