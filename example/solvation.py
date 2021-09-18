from ReadDataLAMMPS import *
from ReplicateSystem import *
from HydrateLAMMPS import *
from WriteDataLAMMPS import *

# reading lammps data file
print('Reading lammps data file(s) ...')
Solute = 'go_opls_long.data'
Solvent = 'water_tip4p_long.data'

SoluteData = readdatalammps(Solute,'full')
SolventData = readdatalammps(Solvent,'full')

print('Importing lammps data file(s) is finished.')

# replication of go for liquid crystal
#datafile = 'mydata.lammps'
#Replicate = 2,2,1
#Init_Atoms,Init_Bonds, Init_Angles, Lx, Ly, Lz = main_replicate(datafile, Replicate)

# solvation
print('Initiat solvation that make take while ...')
lx = 10 # size for the waterbox
ly = 10 # size for the waterbox
lz = 10 # size for the waterbox
d = 2.8 # distance between solute and solvent
solvated = hydrate(SoluteData,SolventData,lx,ly,lz,d)
print('Solvating is finished.')

# writing lammps data file
print('Exporting lammps data file.')
name = 'test.data'
write(solvated,name)
