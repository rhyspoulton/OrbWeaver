// #ifndef HDF5_H
// #define HDF5_H

// #include "orbweaver.h"

// //The number of members(fields) in OrbitData
// #define NFIELDS 26


// /* Calculate the size and the offsets of the struct members in memory */
// OrbitData  tmp;
// size_t dataset_size = sizeof(OrbitData);
// size_t dataset_offsets[NFIELDS] = {
// 	HOFFSET(OrbitData, orbitid),
// 	HOFFSET(OrbitData, entrytype),
// 	HOFFSET(OrbitData, numorbits),
// 	HOFFSET(OrbitData, orbitperiod),
// 	HOFFSET(OrbitData, closestapproach),
// 	HOFFSET(OrbitData, masslossrate),
// 	HOFFSET(OrbitData, scalefactor),
// 	HOFFSET(OrbitData, x),
// 	HOFFSET(OrbitData, y),
// 	HOFFSET(OrbitData, z),
// 	HOFFSET(OrbitData, vx),
// 	HOFFSET(OrbitData, vy),
// 	HOFFSET(OrbitData, vz),
// 	HOFFSET(OrbitData, mvir),
// 	HOFFSET(OrbitData, vmax),
// 	HOFFSET(OrbitData, rmax),
// 	HOFFSET(OrbitData, xrel),
// 	HOFFSET(OrbitData, yrel),
// 	HOFFSET(OrbitData, zrel),
// 	HOFFSET(OrbitData, vxrel),
// 	HOFFSET(OrbitData, vyrel),
// 	HOFFSET(OrbitData, vzrel),
// 	HOFFSET(OrbitData, rvirhost),
// 	HOFFSET(OrbitData, mvirhost),
// 	HOFFSET(OrbitData, vmaxhost),
// 	HOFFSET(OrbitData, rmaxhost)}; 
// size_t dataset_sizes[NFIELDS] = {
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.entrytype),
// 	sizeoff(tmp.numorbits),
// 	sizeoff(tmp.orbitperiod),
// 	sizeoff(tmp.closestapproach),
// 	sizeoff(tmp.masslossrate),
// 	sizeoff(tmp.scalefactor),
// 	sizeoff(tmp.x),
// 	sizeoff(tmp.y),
// 	sizeoff(tmp.z),
// 	sizeoff(tmp.vx),
// 	sizeoff(tmp.vy),
// 	sizeoff(tmp.vz),
// 	sizeoff(tmp.mvir),
// 	sizeoff(tmp.vmax),
// 	sizeoff(tmp.rmax),
// 	sizeoff(tmp.xrel),
// 	sizeoff(tmp.yrel),
// 	sizeoff(tmp.zrel),
// 	sizeoff(tmp.vxrel),
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.orbitid),
// 	sizeoff(tmp.orbitid)};











// #endif // HDF5_H
