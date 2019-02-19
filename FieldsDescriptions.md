| Field  | Units  | Description |  
|---|---|---|
| OrbitID  | N/A | Unique ID to identify this orbit in each file |
| HaloID_orbweaver | N/A | Unique ID to identify this halo in the pre-processed catalog in each file. The 0 values represent interpolated halos that do not exist in the pre-processed catalog |
| HaloID_orig  | N/A | Unique ID to identify this halo in the orignal catlogue. 0 values means it is a interpolated halo that does not exist in the original catalog |
| OrbitedHaloID_orig  | N/A | Unique ID to identify the orbit host halo in the orginal catalogue. 0 values means it is a interpolated halo that does not exist in the original catalog  |
| entrytype  | N/A | This value tell if this entry is either -99 = Apo-center, 99 = pericenter, 0 = mergerpoint or the entrytype * Rvir_host crossing (positive if infalling and negative if outfalling)  |
| numorbits  | N/A | Number of orbits the halo has completed since its first pericentric passage |
| orbitalperiod  | Gyr | Current period of its orbit |
| orbitecc_Wetzel2011  | N/A | The calculated eccentricity of its orbit from [Wetzel, 2011](https://doi.org/10.1111/j.1365-2966.2010.17877.x). |
| orbiteccratio | N/A | The orbital eccetricity found from the peri/apo-centric distances in the simulation |
| orbitalEnergy | solarmasses km^2 / s^2 | The average energy of the orbit since infall or last passage |
| Rperi_Wetzel2011 | Mpc | The calculated peri-centric distance from [Wetzel, 2011](https://doi.org/10.1111/j.1365-2966.2010.17877.x). |
| Rapo_Wetzel2011 | Mpc | The calculated apo-centric distance from [Wetzel, 2011](https://doi.org/10.1111/j.1365-2966.2010.17877.x). |
| closestapproach  | Mpc | Closest approach the halo has had to is host |
| masslossrate_inst  | solar masses/Gyr | The instantaneous rate at which the halo is losing mass (negative values means mass has been accreted) |
| masslossrate_ave  | solar masses/Gyr | The average rate at which the halo is losing mass, this is only calculated at apsis points so will be average since its last passage  (negative values means mass has been accreted) |
| LongAscNode | Radian | The angle of longitude of the ascending node with respect to the intial orbital plane defined [here](https://en.wikipedia.org/wiki/Orbital_elements). |
| Inclination | Radian | The inclination of the halos orbit with respect to the intial orbital plane |
| ArgPeriap | Radian | The argument of periapsis with respect to the intial orbital plane |
| HostAlignment | Radian | The alignment of the orbital angular momentum vector with the host halo's angular momentum vector |
| Phi | Radian | The angle moved through since last passage |
| MergerTimeScale | Gyr | How long the halo takes to merge once crossing 1.0 Rvir of its host halo, this is set the first time the orbiting halo crosses 1.0 Rvir |
| scalefactor  | N/A | Scalefactor of this entry |
| uniage | Gyr | Age of the universe of this entry |
| X  | Mpc | X position of the halo in the simulation |
| Y  | Mpc | Y position of the halo in the simulation |
| Z  | Mpc | Z position of the halo in the simulation |
| VX  | km/s | X component of the halos velocity in the simulation |
| VY  | km/s | Y component of the halos velocity in the simulation |
| VZ  | km/s | Z component of the halos velocity in the simulation |
| Npart | N/A | Number of particles in the orbiting halo |
| Mass  | solar masses | Mass of the halo (depends on mass definition given) |
| Rmax  | Mpc | Radial distance of Vmax |
| Vmax  | km/s | Maxiumum circular velocity of the halo |
| Vmaxpeak | km/s | The peak Vmax has had in its existence up to the current entry time |
| cNFW | N/A  | Concentration of the halo |
| fieldHalo | N/A | Flag if this orbiting halo is top of its spatial hierachy (not a subhalo) ,where: 0 = No, 1 = Yes |
| Vrad | km/s | The radial velocity of the orbiting halo with respect to its host |
| Vtan | km/s | The tangential velocity of the orbiting halo with respect to its host |
| Xrel  | Mpc | X position of the halo relative to its host |
| Yrel  | Mpc | Y position of the halo relative to its host |
| Zrel  | Mpc | Z position of the halo relative to its host |
| VXrel  | km/s | X component of the halos velocity relative to its host |
| VYrel  | km/s | Y component of the halos velocity relative to its host |
| Vzrel  | km/s | Z component of the halos velocity relative to its host |
| Lxrel | solar masses Mpc km/s | The x comoponent of average  orbital angular momentum vector |
| Lyrel | solar masses Mpc km/s | The x comoponent of average  orbital angular momentum vector |
| Lzrel | solar masses Mpc km/s | The x comoponent of average orbital angular momentum vector |
|Npart_host | N/A | Number of particle in its host halo |
| Mass_host  | solar masses | Mass of its host (depends on mass definition given) |
| R\_200crit_host  | Mpc | Virial radius of its host |
| Rmax_host  | Mpc | Radial distance of Vmax for its host |
| Vmax_host  | km/s | Maximum circular velocity of its host
| cNFW_host | N/A  | Concentration of its host |
|fieldhalo_host | N/A | If this halo host is the top of its spatial hierachy (not a subhalo), where: 0 = No, 1 = Yes |

