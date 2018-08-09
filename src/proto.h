
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
void ProcessOrbits(Options &opt, SnapData *&snapdata, vector< vector<OrbitData>> &orbitdata);
void ProcessHalo(Int_t snap, Int_t i, Options &opt, SnapData *&snapdata, vector<vector<OrbitData>> &orbitdata);
OrbitData CalcOrbitProps(HaloData &orbitinghalo, HaloData &hosthalo);
#endif //ifndef ORBPROTO_H