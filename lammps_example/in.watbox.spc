units		   real
atom_style	  full
boundary		p p p
dielectric	  1
special_bonds   lj/coul 0.0 0.0 1.0 

pair_style	  lj/charmm/coul/long/opt 9 10.00000 
bond_style      harmonic  
angle_style     harmonic  
dihedral_style  none
improper_style  none
kspace_style	pppm 1e-05

read_data	   data.watbox.spc

pair_modify	 mix geometric
neighbor		2.0 multi
neigh_modify	every 2 delay 4 check yes
thermo_style	multi
thermo_modify		line multi format float %14.6f
variable		input index in.watbox.spc
variable		sname index watbox.spc

print                          .
print ==========================================
print "500 steps CG Minimization"
print ==========================================
print                          .

dump            1 all atom 25 ${sname}.min.lammpstrj
dump_modify     1 image yes scale yes
thermo          10
min_style       sd
minimize        1.0e-4 1.0e-4 500 5000
min_style       cg
minimize        1.0e-4 1.0e-4 500 5000
#now minimize the entire system
minimize        1.0e-4 1.0e-4 500 5000
undump          1


print                          .
print =====================================
print "NVT dynamics to heat system"
print =====================================
print            .

reset_timestep  0
timestep        1.0
fix             shakeH all shake 0.0001 20 500 m 1.008 a 1
velocity        all create 0.0 12345678 dist uniform
thermo          100
thermo_style    multi
timestep        1.0
dump            1 all custom 1000 ${sname}.heat.lammpstrj id type xu yu zu vx vy vz
fix             4 all nvt temp 1.0 300.0 100.0
run             10000
unfix           4
undump          1

print                          .
print ================================================
print "NPT dynamics with an isotropic pressure of 1atm."
print ================================================
print                       .

timestep        2.0
fix             2 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 2000.0
thermo          100
thermo_style    multi
dump            1 all custom 5000 ${sname}.npt.lammps id type xu yu zu vx vy vz
run             100000 # run for 15 ns
unfix           2
undump          1

log             ${sname}.2pt.eng
timestep	2.0

#compute         atomPE all pe/atom
#compute         atomKE all ke/atom
#variable        atomEng atom c_atomPE+c_atomKE

print "================================================"
print "NVT dynamics for 20ps dumping velocities"
print "================================================"
thermo          2
thermo_style    custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press vol
thermo_modify	line multi
fix             1 all nvt temp 300.0 300.0 100.0
dump            1 all custom 2 ${sname}.2pt.lammps id type xu yu zu vx vy vz 
run             10000
undump          4
undump          1
unfix           1
#uncompute       atomPE
#uncompute       atomKE

