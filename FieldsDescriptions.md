| Field | Units | Output point (Apsis or Crossing point or Both) | Description |
|------------------------------|----------------------------|-----------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| OrbitID | N/A | Both | Unique ID to identify this orbit in each file |
| HaloID_orbweaver | N/A | Both | Unique ID to identify this halo in the pre-processed catalog in each file. The 0 values represent interpolated halos that do not exist in the pre-processed catalog |
| HaloID_orig | N/A | Both | Unique ID to identify this halo in the orignal catlogue. 0 values means it is a interpolated halo that does not exist in the original catalog |
| HaloRootProgen_orig | N/A | Both | Unique ID to indentify the root progenitor for the halo from the original halo catlog |
| HaloRootDescen_orig | N/A | Both | Unique ID to indentify the root descendant for the halo from the original halo catlog |
| OrbitedHaloID_orig | N/A | Both | Unique ID to identify the orbit host halo in the orginal catalogue. 0 values means it is a interpolated halo that does not exist in the original catalog |
| OrbitedHaloRootProgen_orig | N/A | Both | Unique ID to identify the original root progenitor of the orbited halo in the halo catalog. This can also be used to find any halos that have orbited this object, by finding the halos that have the same OrbitedHaloRootProgen_orig ID. |
| OrbitedHaloRootDescen_orig | N/A | Both | Unique ID to identify the original root descendant of the orbited halo in the halo catalog. This can be used to extract the super set of halo that orbit this host and also the halos that orbit any object that merges with this host, by finding the halos that have the same OrbitedHaloRootDescen_orig ID. |
| entrytype | N/A | Both | This value states if this entry is either: <br> **-99** = Apo-center, <br> **99** = pericenter, <br> **0** = endpoint of a orbit (the halo has either; terminated/merged with another halo, merged with its host [see MergedFlag] or host has terminated/ merged [see hostMergedFlag] ) <br> **All other values** = the fraction of host crossing i.e. entrytype * R\_200crit\_host crossing (positive if infalling and negative if outfalling). To get the correct number of R\_200crit\_host, this dataset will need to be rounded to the desired number of decimals.  |
| num_entrytype | N/A | Both | This values tells you the number of each entry type so far in the orbit, such that if you want to extract the first crossing of rvir you can query the whole dataset using: entrytype==1.0 & num_entrytype==1 |
| numorbits | N/A | Both | Number of orbits the halo has completed since its first pericentric passage |
| totnumorbits | N/A | Both | The total number of orbits that the halo has completed in the simulation |
| orbitalperiod | Gyr | Apsis | Current period of its orbit |
| orbitecc_ratio | N/A | Apsis | The orbital eccetricity found from the peri/apo-centric distances in the simulation |
| orbitalenergy_inst | solarmasses km^2 / s^2 | Both | The instantaneous energy of the orbit |
| orbitalenergy_ave | solarmasses km^2 / s^2 | Apsis | The average energy of the orbit since infall or last passage, only outputted at apsis points |
| R_circ | phys Mpc | Both | The radius of a circular orbit with the same orbital energy, calculated from Khochfar and Burkert 2006 |
| V_circ | km/s | Both | The velocity of a circular orbit with the same orbital energy, calculated from Khochfar and Burkert 2006 |
| L_circ | solar masses phys Mpc km/s | Both | The orbital angular momentum of a circular orbit with the same orbital energy, calculated from Khochfar and Burkert 2006 |
| Eta | N/A | Both | The ratio of the (instantaneous) orbital angular momentum to the orbital angularmomentum of a circular orbit with the same orbital energy (L\_circ). This is useful to identify the type of orbit the halo is on, where 0 is a highly radial orbit and 1 is a circular orbit. L\_circ is calculated from Khochfar and Burkert 2006 |
| Rperi_calc | phys Mpc | Both | The calculated peri-centric distance from Wetzel, 2011. |
| Rapo_calc | phys Mpc | Both | The calculated apo-centric distance from Wetzel, 2011. |
| orbitecc_calc | N/A | Both | The calculated eccentricity of its orbit from Khochfar and Burkert 2006. |
| closestapproach | phys Mpc | Both | Closest approach the halo has had to is host up to this point in its orbit |
| closestapproachscalefactor | N/A | Both | The scalefactor which the closest approach occured at |
| masslossrate_inst | solar masses/Gyr | Both | The instantaneous rate at which the halo is losing mass (negative values means mass has been accreted) |
| masslossrate_ave | solar masses/Gyr | Apsis | The average rate at which the halo is losing mass, this is only calculated at apsis points so will be average since its last passage (negative values means mass has been accreted) |
| LongAscNode | Radian | Apsis | The angle of longitude of the ascending node with respect to the intial orbital plane defined here. |
| Inclination | Radian | Apsis | The inclination of the halos orbit with respect to the intial orbital plane |
| ArgPeriap | Radian | Apsis | The argument of periapsis with respect to the intial orbital plane |
| HostAlignment | Radian | Apsis | The alignment of the orbital angular momentum vector with the host halo's angular momentum vector |
| Phi | Radian | Apsis | The angle moved through since last passage |
| scalefactor | N/A | Both | Scalefactor of this entry |
| uniage | Gyr | Both | Age of the universe of this entry |
| X | phys Mpc | Both | X position of the halo in the simulation |
| Y | phys Mpc | Both | Y position of the halo in the simulation |
| Z | phys Mpc | Both | Z position of the halo in the simulation |
| VX | km/s | Both | X component of the halos velocity in the simulation |
| VY | km/s | Both | Y component of the halos velocity in the simulation |
| VZ | km/s | Both | Z component of the halos velocity in the simulation |
| npart | N/A | Both | Number of particles in the orbiting halo |
| Mass | solar masses | Both | Mass of the halo (depends on mass definition given) |
| Radius | phys Mpc | Both | Radius of the halo (depends on mass definition given) |
| min_Rscale | phys Mpc | Both | Minimum of the scale radius in this orbit history, used to see if the halo can be disrupted artificially due to inadequate force softening |
| min_Rmax | phys Mpc | Both | Minimum of the radius that the maximu circular velocity is at in this orbit history, used to see if the halo can be disrupted artificially due to inadequate force softening. |
| Rmax | phys Mpc | Both | Radial distance of Vmax |
| Vmax | km/s | Both | Maxiumum circular velocity of the halo |
| Vmaxpeak | km/s | Both | The peak Vmax has had in its existence up to the current entry time |
| cNFW | N/A | Both | Concentration of the halo |
| fieldHalo | N/A | Both | Flag if this orbiting halo is top of its spatial hierachy (not a subhalo) ,where: 0 = No, 1 = Yes |
| numSubStruct | N/A | Both | The number of substructure that this halo contains |
| RatioOfMassinSubsStruct | N/A | Both | The ratio of how much this halo's mass is in substructure |
| MergerTimeScale | Gyr | Crossing Point (see description) | How long the halo takes to merge once crossing 1.0 Rvir of its host halo, this is set the first time the orbiting halo crosses 1.0 Rvir (entrytype==1.0 & num_entrytype==1)|
| MergedFlag | N/A | Both | A flag identifying if this halo merges with its orbit host |
| Vrad | km/s | Both | The radial velocity of the orbiting halo with respect to its host |
| Vtan | km/s | Both | The tangential velocity of the orbiting halo with respect to its host |
| R | phys Mpc | Both | The distance of the halo relative to its host |
| Xrel | phys Mpc | Both | X position of the halo relative to its host |
| Yrel | phys Mpc | Both | Y position of the halo relative to its host |
| Zrel | phys Mpc | Both | Z position of the halo relative to its host |
| Vrel | km/s | Both | The halos velocity relative to its host |
| VXrel | km/s | Both | X component of the halos velocity relative to its host |
| VYrel | km/s | Both | Y component of the halos velocity relative to its host |
| Vzrel | km/s | Both | Z component of the halos velocity relative to its host |
| Lrel_inst | solar masses phys Mpc km/s | Both | The instantaneous orbital angular momentum vector |
| LXrel_inst | solar masses phys Mpc km/s | Both | The instantaneous x component of the orbital angular momentum vector |
| LYrel_inst | solar masses phys Mpc km/s | Both | The instantaneous y component of the orbital angular momentum vector |
| LZrel_inst | solar masses phys Mpc km/s | Both | The instantaneous z component of the orbital angular momentum vector |
| Lrel_ave | solar masses phys Mpc km/s | Apsis | The average orbital angular momentum vector since last apsis point |
| LXrel_ave | solar masses phys Mpc km/s | Apsis | The average x component of the orbital angular momentum vector since last apsis point |
| LYrel_ave | solar masses phys pc km/s | Apsis | The average y component of the orbital angular momentum vector since last apsis point |
| LZrel_ave | solar masses phys Mpc km/s | Apsis | The average z component of the orbital angular momentum vector since last apsis point |
| Npart_host | N/A | Both | Number of particle in its host halo |
| Mass_host | solar masses | Both | Mass of its host (depends on mass definition given) |
| Radius_host | phys Mpc | Both | Radius of its host (depends on mass definition given) |
| Rmax_host | phys Mpc | Both | Radial distance of Vmax for its host |
| Vmax_host | km/s | Both | Maximum circular velocity of its host |
| cNFW_host | N/A | Both | Concentration of its host |
| fieldhalo_host | N/A | Both | If this halo host is the top of its spatial hierachy (not a subhalo), where: 0 = No, 1 = Yes |
| hostMerges | N/A | Both | A flag to indicate if the orbit host merges with another halo or not where 1 = merges, 0 = does not merge |
| numSubStruct_host | N/A | Both | The number of substructure that its host contains |
| RatioOfMassinSubsStruct_host | N/A | Both | The ratio of how much of its host's mass is in substructure |