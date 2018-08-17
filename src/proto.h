
#ifndef ORBPROTO_H
#define ORBPROTO_H

#include "allvars.h"

bool FileExists(const char *fname);
void GetArgs(int argc, char *argv[], Options &opt);
void ConfigCheck(Options &opt);
void ReadHeader(H5File *Fhdf,HDFCatalogNames hdfnames);
void ReadData(Options &opt, SnapData *&snapdata);
vector<HaloData> ReadSnapshotData(int snap, int i,Group snapgroup, Options &opt, SnapData *&snapdata, HDFCatalogNames hdfnames);

double GetUniverseAge(double scalefactor);

HaloData InterpHaloProps(Options &opt, vector<Int_t> &halosnaps, vector<Int_t> &haloindexes, SnapData *&snapdata);
void InterpPassagePoints(vector<Int_t> halosnaps,vector<Int_t> haloindexes,vector<Int_t> hostindexes, SnapData *&snapdata, vector<OrbitData> &branchorbitdata);
void ProcessOrbits(Options &opt, SnapData *&snapdata, vector<OrbitData> &orbitdata);
void ProcessHalo(Int_t snap, Int_t i, Options &opt, SnapData *&snapdata, OrbitData &orbitdata);
void CalcOrbitProps(Int_t orbitID, int currentsnap, int prevsnap, HaloData &orbitinghalo, HaloData &hosthalo, HaloData &prevorbitinghalo, HaloData &prevhosthalo, vector<OrbitData> &orbitdata, OrbitData &tmporbitdata, SnapData *&snapdata, OrbitProps &orbitprops);

void WriteOrbitData(Options &opt, vector<OrbitData> &orbitdata);
#endif //ifndef ORBPROTO_H