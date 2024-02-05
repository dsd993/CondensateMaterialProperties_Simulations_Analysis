#importing the libraries
import sys, os, numpy as np
import hoomd, hoomd.md as md
from hoomd import azplugins
import gsd, gsd.hoomd, gsd.pygsd
nsteps = 300000000 #simulation runtime length

temp=280 #simulation temperature

hoomd.context.initialize() #Initialize the execution context
system = hoomd.init.read_gsd(sys.argv[1]) #Read initial system state from an GSD file.

nl = hoomd.md.nlist.cell() #Cell list based neighbor list
#### BOND DATA ####
harmonic = hoomd.md.bond.harmonic() #Harmonic bond potential and its parameters
harmonic.bond_coeff.set('1', k=20.000000, r0=3.800000)

nl.reset_exclusions(exclusions=['1-2', 'body']) #setting the exclusions from short range pair interactions
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl) #Ashbaugh-Hatch potential and its parameters
nb.pair_coeff.set('1', '1', epsilon = 0.200000, sigma = 6.180000, lam = 0.596471, r_cut = 24.720000)
nb.pair_coeff.set('1', '2', epsilon = 0.200000, sigma = 6.050000, lam = 0.258236, r_cut = 24.200000)
nb.pair_coeff.set('1', '3', epsilon = 0.200000, sigma = 5.680000, lam = 0.552354, r_cut = 22.720000)
nb.pair_coeff.set('1', '4', epsilon = 0.200000, sigma = 5.930000, lam = 0.552354, r_cut = 23.720000)
nb.pair_coeff.set('1', '5', epsilon = 0.200000, sigma = 6.100000, lam = 0.537648, r_cut = 24.400000)
nb.pair_coeff.set('1', '6', epsilon = 0.200000, sigma = 5.340000, lam = 0.545000, r_cut = 21.360000)
nb.pair_coeff.set('1', '7', epsilon = 0.200000, sigma = 5.610000, lam = 0.559707, r_cut = 22.440000)
nb.pair_coeff.set('1', '8', epsilon = 0.200000, sigma = 6.180000, lam = 0.618530, r_cut = 24.720000)
nb.pair_coeff.set('1', '9', epsilon = 0.200000, sigma = 6.370000, lam = 0.537648, r_cut = 25.480000)
nb.pair_coeff.set('1', '10', epsilon = 0.200000, sigma = 6.320000, lam = 0.706765, r_cut = 25.280000)
nb.pair_coeff.set('1', '11', epsilon = 0.200000, sigma = 6.020000, lam = 0.590589, r_cut = 24.080000)
nb.pair_coeff.set('1', '12', epsilon = 0.200000, sigma = 5.870000, lam = 0.637648, r_cut = 23.480000)
nb.pair_coeff.set('1', '13', epsilon = 0.200000, sigma = 6.130000, lam = 0.640589, r_cut = 24.520000)
nb.pair_coeff.set('1', '14', epsilon = 0.200000, sigma = 5.880000, lam = 0.405295, r_cut = 23.520000)
nb.pair_coeff.set('1', '15', epsilon = 0.200000, sigma = 6.270000, lam = 0.670000, r_cut = 25.080000)
nb.pair_coeff.set('2', '2', epsilon = 0.200000, sigma = 5.920000, lam = -0.080000, r_cut = 23.680000)
nb.pair_coeff.set('2', '3', epsilon = 0.200000, sigma = 5.550000, lam = 0.214118, r_cut = 22.200000)
nb.pair_coeff.set('2', '4', epsilon = 0.200000, sigma = 5.800000, lam = 0.214118, r_cut = 23.200000)
nb.pair_coeff.set('2', '5', epsilon = 0.200000, sigma = 5.970000, lam = 0.199412, r_cut = 23.880000)
nb.pair_coeff.set('2', '6', epsilon = 0.200000, sigma = 5.210000, lam = 0.206765, r_cut = 20.840000)
nb.pair_coeff.set('2', '7', epsilon = 0.200000, sigma = 5.480000, lam = 0.221471, r_cut = 21.920000)
nb.pair_coeff.set('2', '8', epsilon = 0.200000, sigma = 6.050000, lam = 0.280294, r_cut = 24.200000)
nb.pair_coeff.set('2', '9', epsilon = 0.200000, sigma = 6.240000, lam = 0.199412, r_cut = 24.960000)
nb.pair_coeff.set('2', '10', epsilon = 0.200000, sigma = 6.190000, lam = 0.368529, r_cut = 24.760000)
nb.pair_coeff.set('2', '11', epsilon = 0.200000, sigma = 5.890000, lam = 0.252353, r_cut = 23.560000)
nb.pair_coeff.set('2', '12', epsilon = 0.200000, sigma = 5.740000, lam = 0.299412, r_cut = 22.960000)
nb.pair_coeff.set('2', '13', epsilon = 0.200000, sigma = 6.000000, lam = 0.302354, r_cut = 24.000000)
nb.pair_coeff.set('2', '14', epsilon = 0.200000, sigma = 5.750000, lam = 0.067060, r_cut = 23.000000)
nb.pair_coeff.set('2', '15', epsilon = 0.200000, sigma = 6.140000, lam = 0.331765, r_cut = 24.560000)
nb.pair_coeff.set('3', '3', epsilon = 0.200000, sigma = 5.180000, lam = 0.508236, r_cut = 20.720000)
nb.pair_coeff.set('3', '4', epsilon = 0.200000, sigma = 5.430000, lam = 0.508236, r_cut = 21.720000)
nb.pair_coeff.set('3', '5', epsilon = 0.200000, sigma = 5.600000, lam = 0.493530, r_cut = 22.400000)
nb.pair_coeff.set('3', '6', epsilon = 0.200000, sigma = 4.840000, lam = 0.500883, r_cut = 19.360000)
nb.pair_coeff.set('3', '7', epsilon = 0.200000, sigma = 5.110000, lam = 0.515589, r_cut = 20.440000)
nb.pair_coeff.set('3', '8', epsilon = 0.200000, sigma = 5.680000, lam = 0.574413, r_cut = 22.720000)
nb.pair_coeff.set('3', '9', epsilon = 0.200000, sigma = 5.870000, lam = 0.493530, r_cut = 23.480000)
nb.pair_coeff.set('3', '10', epsilon = 0.200000, sigma = 5.820000, lam = 0.662648, r_cut = 23.280000)
nb.pair_coeff.set('3', '11', epsilon = 0.200000, sigma = 5.520000, lam = 0.546472, r_cut = 22.080000)
nb.pair_coeff.set('3', '12', epsilon = 0.200000, sigma = 5.370000, lam = 0.593530, r_cut = 21.480000)
nb.pair_coeff.set('3', '13', epsilon = 0.200000, sigma = 5.630000, lam = 0.596471, r_cut = 22.520000)
nb.pair_coeff.set('3', '14', epsilon = 0.200000, sigma = 5.380000, lam = 0.361178, r_cut = 21.520000)
nb.pair_coeff.set('3', '15', epsilon = 0.200000, sigma = 5.770000, lam = 0.625883, r_cut = 23.080000)
nb.pair_coeff.set('4', '4', epsilon = 0.200000, sigma = 5.680000, lam = 0.508236, r_cut = 22.720000)
nb.pair_coeff.set('4', '5', epsilon = 0.200000, sigma = 5.850000, lam = 0.493530, r_cut = 23.400000)
nb.pair_coeff.set('4', '6', epsilon = 0.200000, sigma = 5.090000, lam = 0.500883, r_cut = 20.360000)
nb.pair_coeff.set('4', '7', epsilon = 0.200000, sigma = 5.360000, lam = 0.515589, r_cut = 21.440000)
nb.pair_coeff.set('4', '8', epsilon = 0.200000, sigma = 5.930000, lam = 0.574413, r_cut = 23.720000)
nb.pair_coeff.set('4', '9', epsilon = 0.200000, sigma = 6.120000, lam = 0.493530, r_cut = 24.480000)
nb.pair_coeff.set('4', '10', epsilon = 0.200000, sigma = 6.070000, lam = 0.662648, r_cut = 24.280000)
nb.pair_coeff.set('4', '11', epsilon = 0.200000, sigma = 5.770000, lam = 0.546472, r_cut = 23.080000)
nb.pair_coeff.set('4', '12', epsilon = 0.200000, sigma = 5.620000, lam = 0.593530, r_cut = 22.480000)
nb.pair_coeff.set('4', '13', epsilon = 0.200000, sigma = 5.880000, lam = 0.596471, r_cut = 23.520000)
nb.pair_coeff.set('4', '14', epsilon = 0.200000, sigma = 5.630000, lam = 0.361178, r_cut = 22.520000)
nb.pair_coeff.set('4', '15', epsilon = 0.200000, sigma = 6.020000, lam = 0.625883, r_cut = 24.080000)
nb.pair_coeff.set('5', '5', epsilon = 0.200000, sigma = 6.020000, lam = 0.478824, r_cut = 24.080000)
nb.pair_coeff.set('5', '6', epsilon = 0.200000, sigma = 5.260000, lam = 0.486177, r_cut = 21.040000)
nb.pair_coeff.set('5', '7', epsilon = 0.200000, sigma = 5.530000, lam = 0.500883, r_cut = 22.120000)
nb.pair_coeff.set('5', '8', epsilon = 0.200000, sigma = 6.100000, lam = 0.559707, r_cut = 24.400000)
nb.pair_coeff.set('5', '9', epsilon = 0.200000, sigma = 6.290000, lam = 0.478824, r_cut = 25.160000)
nb.pair_coeff.set('5', '10', epsilon = 0.200000, sigma = 6.240000, lam = 0.647942, r_cut = 24.960000)
nb.pair_coeff.set('5', '11', epsilon = 0.200000, sigma = 5.940000, lam = 0.531766, r_cut = 23.760000)
nb.pair_coeff.set('5', '12', epsilon = 0.200000, sigma = 5.790000, lam = 0.578824, r_cut = 23.160000)
nb.pair_coeff.set('5', '13', epsilon = 0.200000, sigma = 6.050000, lam = 0.581766, r_cut = 24.200000)
nb.pair_coeff.set('5', '14', epsilon = 0.200000, sigma = 5.800000, lam = 0.346471, r_cut = 23.200000)
nb.pair_coeff.set('5', '15', epsilon = 0.200000, sigma = 6.190000, lam = 0.611177, r_cut = 24.760000)
nb.pair_coeff.set('6', '6', epsilon = 0.200000, sigma = 4.500000, lam = 0.493530, r_cut = 18.000000)
nb.pair_coeff.set('6', '7', epsilon = 0.200000, sigma = 4.770000, lam = 0.508236, r_cut = 19.080000)
nb.pair_coeff.set('6', '8', epsilon = 0.200000, sigma = 5.340000, lam = 0.567060, r_cut = 21.360000)
nb.pair_coeff.set('6', '9', epsilon = 0.200000, sigma = 5.530000, lam = 0.486177, r_cut = 22.120000)
nb.pair_coeff.set('6', '10', epsilon = 0.200000, sigma = 5.480000, lam = 0.655295, r_cut = 21.920000)
nb.pair_coeff.set('6', '11', epsilon = 0.200000, sigma = 5.180000, lam = 0.539119, r_cut = 20.720000)
nb.pair_coeff.set('6', '12', epsilon = 0.200000, sigma = 5.030000, lam = 0.586177, r_cut = 20.120000)
nb.pair_coeff.set('6', '13', epsilon = 0.200000, sigma = 5.290000, lam = 0.589119, r_cut = 21.160000)
nb.pair_coeff.set('6', '14', epsilon = 0.200000, sigma = 5.040000, lam = 0.353824, r_cut = 20.160000)
nb.pair_coeff.set('6', '15', epsilon = 0.200000, sigma = 5.430000, lam = 0.618530, r_cut = 21.720000)
nb.pair_coeff.set('7', '7', epsilon = 0.200000, sigma = 5.040000, lam = 0.522942, r_cut = 20.160000)
nb.pair_coeff.set('7', '8', epsilon = 0.200000, sigma = 5.610000, lam = 0.581765, r_cut = 22.440000)
nb.pair_coeff.set('7', '9', epsilon = 0.200000, sigma = 5.800000, lam = 0.500883, r_cut = 23.200000)
nb.pair_coeff.set('7', '10', epsilon = 0.200000, sigma = 5.750000, lam = 0.670000, r_cut = 23.000000)
nb.pair_coeff.set('7', '11', epsilon = 0.200000, sigma = 5.450000, lam = 0.553824, r_cut = 21.800000)
nb.pair_coeff.set('7', '12', epsilon = 0.200000, sigma = 5.300000, lam = 0.600883, r_cut = 21.200000)
nb.pair_coeff.set('7', '13', epsilon = 0.200000, sigma = 5.560000, lam = 0.603825, r_cut = 22.240000)
nb.pair_coeff.set('7', '14', epsilon = 0.200000, sigma = 5.310000, lam = 0.368531, r_cut = 21.240000)
nb.pair_coeff.set('7', '15', epsilon = 0.200000, sigma = 5.700000, lam = 0.633236, r_cut = 22.800000)
nb.pair_coeff.set('8', '8', epsilon = 0.200000, sigma = 6.180000, lam = 0.640589, r_cut = 24.720000)
nb.pair_coeff.set('8', '9', epsilon = 0.200000, sigma = 6.370000, lam = 0.559707, r_cut = 25.480000)
nb.pair_coeff.set('8', '10', epsilon = 0.200000, sigma = 6.320000, lam = 0.728824, r_cut = 25.280000)
nb.pair_coeff.set('8', '11', epsilon = 0.200000, sigma = 6.020000, lam = 0.612648, r_cut = 24.080000)
nb.pair_coeff.set('8', '12', epsilon = 0.200000, sigma = 5.870000, lam = 0.659707, r_cut = 23.480000)
nb.pair_coeff.set('8', '13', epsilon = 0.200000, sigma = 6.130000, lam = 0.662648, r_cut = 24.520000)
nb.pair_coeff.set('8', '14', epsilon = 0.200000, sigma = 5.880000, lam = 0.427354, r_cut = 23.520000)
nb.pair_coeff.set('8', '15', epsilon = 0.200000, sigma = 6.270000, lam = 0.692060, r_cut = 25.080000)
nb.pair_coeff.set('9', '9', epsilon = 0.200000, sigma = 6.560000, lam = 0.478824, r_cut = 26.240000)
nb.pair_coeff.set('9', '10', epsilon = 0.200000, sigma = 6.510000, lam = 0.647942, r_cut = 26.040000)
nb.pair_coeff.set('9', '11', epsilon = 0.200000, sigma = 6.210000, lam = 0.531766, r_cut = 24.840000)
nb.pair_coeff.set('9', '12', epsilon = 0.200000, sigma = 6.060000, lam = 0.578824, r_cut = 24.240000)
nb.pair_coeff.set('9', '13', epsilon = 0.200000, sigma = 6.320000, lam = 0.581766, r_cut = 25.280000)
nb.pair_coeff.set('9', '14', epsilon = 0.200000, sigma = 6.070000, lam = 0.346471, r_cut = 24.280000)
nb.pair_coeff.set('9', '15', epsilon = 0.200000, sigma = 6.460000, lam = 0.611177, r_cut = 25.840000)
nb.pair_coeff.set('10', '10', epsilon = 0.200000, sigma = 6.460000, lam = 0.817059, r_cut = 25.840000)
nb.pair_coeff.set('10', '11', epsilon = 0.200000, sigma = 6.160000, lam = 0.700883, r_cut = 24.640000)
nb.pair_coeff.set('10', '12', epsilon = 0.200000, sigma = 6.010000, lam = 0.747942, r_cut = 24.040000)
nb.pair_coeff.set('10', '13', epsilon = 0.200000, sigma = 6.270000, lam = 0.750883, r_cut = 25.080000)
nb.pair_coeff.set('10', '14', epsilon = 0.200000, sigma = 6.020000, lam = 0.515589, r_cut = 24.080000)
nb.pair_coeff.set('10', '15', epsilon = 0.200000, sigma = 6.410000, lam = 0.780295, r_cut = 25.640000)
nb.pair_coeff.set('11', '11', epsilon = 0.200000, sigma = 5.860000, lam = 0.584707, r_cut = 23.440000)
nb.pair_coeff.set('11', '12', epsilon = 0.200000, sigma = 5.710000, lam = 0.631766, r_cut = 22.840000)
nb.pair_coeff.set('11', '13', epsilon = 0.200000, sigma = 5.970000, lam = 0.634707, r_cut = 23.880000)
nb.pair_coeff.set('11', '14', epsilon = 0.200000, sigma = 5.720000, lam = 0.399413, r_cut = 22.880000)
nb.pair_coeff.set('11', '15', epsilon = 0.200000, sigma = 6.110000, lam = 0.664119, r_cut = 24.440000)
nb.pair_coeff.set('12', '12', epsilon = 0.200000, sigma = 5.560000, lam = 0.678824, r_cut = 22.240000)
nb.pair_coeff.set('12', '13', epsilon = 0.200000, sigma = 5.820000, lam = 0.681765, r_cut = 23.280000)
nb.pair_coeff.set('12', '14', epsilon = 0.200000, sigma = 5.570000, lam = 0.446472, r_cut = 22.280000)
nb.pair_coeff.set('12', '15', epsilon = 0.200000, sigma = 5.960000, lam = 0.711177, r_cut = 23.840000)
nb.pair_coeff.set('13', '13', epsilon = 0.200000, sigma = 6.080000, lam = 0.684707, r_cut = 24.320000)
nb.pair_coeff.set('13', '14', epsilon = 0.200000, sigma = 5.830000, lam = 0.449413, r_cut = 23.320000)
nb.pair_coeff.set('13', '15', epsilon = 0.200000, sigma = 6.220000, lam = 0.714119, r_cut = 24.880000)
nb.pair_coeff.set('14', '14', epsilon = 0.200000, sigma = 5.580000, lam = 0.214119, r_cut = 22.320000)
nb.pair_coeff.set('14', '15', epsilon = 0.200000, sigma = 5.970000, lam = 0.478825, r_cut = 23.880000)
nb.pair_coeff.set('15', '15', epsilon = 0.200000, sigma = 6.360000, lam = 0.743530, r_cut = 25.440000)

### ELECTROSTATICS ###

yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
yukawa.pair_coeff.set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'], ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'], epsilon=0.000000, kappa=0.000000, r_cut=0.000000)
yukawa.pair_coeff.set('2','2', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','9', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('2','14', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','2', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','9', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('9','14', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('14','2', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('14','9', epsilon=-4.150796, kappa=0.100000, r_cut=35.000000)
yukawa.pair_coeff.set('14','14', epsilon=4.150796, kappa=0.100000, r_cut=35.000000)

### MAKE PARTICLE GROUPS ###
all=hoomd.group.all()
nonrigid = hoomd.group.nonrigid()
rigid = hoomd.group.rigid()
ghost = hoomd.group.rigid_center()
integrate_group = hoomd.group.union(name='int_grp',a=nonrigid,b=ghost)
outp_group = hoomd.group.difference(name='outp',a=all, b=ghost)
### Slab extension and NVT Integration ###
num_particles=len(system.particles)
print(num_particles)
position = np.zeros((num_particles,3))
#manipulating the snapshot to unwrap the particles in the desired dimension
for i in range(num_particles):
   position[i,0] = system.particles[i].position[0]
   position[i,1] = system.particles[i].position[1]
   position[i,2] = system.particles[i].position[2] + system.box.Lz*system.particles[i].image[2]
snap = system.take_snapshot(all=True)
#updating the particles positions after unwrapping
for i in range(num_particles):
   snap.particles.position[i] = position[i]
# setting the new z box size in the snapshot
snap.box = hoomd.data.boxdim(Lx=system.box.Lx, Ly=system.box.Lx, Lz=system.box.Lx*8, xy=0, xz=0, yz=0)
system.restore_snapshot(snap)
hoomd.md.integrate.mode_standard(dt=0.2045814) #time step is 10fs
integrator = hoomd.md.integrate.langevin(group=integrate_group, kT=temp*0.0019872067,seed=555) #Langevin integrator
#setting the friction factor for the Langevin integrator [mass/1000fs]
integrator.set_gamma('1',gamma=6.413095/1000)
integrator.set_gamma('2',gamma=6.310447/1000)
integrator.set_gamma('3',gamma=4.256497/1000)
integrator.set_gamma('4',gamma=5.577242/1000)
integrator.set_gamma('5',gamma=6.261567/1000)
integrator.set_gamma('6',gamma=2.788621/1000)
integrator.set_gamma('7',gamma=3.474412/1000)
integrator.set_gamma('8',gamma=5.533250/1000)
integrator.set_gamma('9',gamma=7.635103/1000)
integrator.set_gamma('10',gamma=7.977265/1000)
integrator.set_gamma('11',gamma=4.842571/1000)
integrator.set_gamma('12',gamma=4.747255/1000)
integrator.set_gamma('13',gamma=6.701489/1000)
integrator.set_gamma('14',gamma=5.626122/1000)
integrator.set_gamma('15',gamma=7.195180/1000)
writegsd=hoomd.dump.gsd('restart.gsd',period=1000000, group=all,overwrite=True,truncate=True,dynamic=['property', 'momentum']) #writes the simulation snapshot for restarting purposes
hoomd.analyze.log(filename='thermo.log', quantities=['potential_energy','pair_ashbaugh_energy', 'pair_yukawa_energy','bond_harmonic_energy', 'temperature','lx','ly','lz'], period=1000000, overwrite=True, header_prefix='#') #Log a number of calculated quantities to a file.
hoomd.analyze.log(filename='pressure.log', quantities=['pressure_xx', 'pressure_yy', 'pressure_zz'], period=1000, overwrite=True, header_prefix='#')
hoomd.dump.gsd('prod.gsd',period=1000000, group=outp_group,overwrite=True,truncate=True,dynamic=['property', 'momentum']) #writes the last state of the simulation
hoomd.dump.gsd('prodtraj.gsd',period=500000, group=outp_group, overwrite=True,dynamic=['property', 'momentum']) #writes the trajectory of the simulation
hoomd.run_upto(nsteps) #Runs the simulation up to a given time step number.
writegsd.write_restart()