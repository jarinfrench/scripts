units metal
atom_style atomic
boundary f f f

variable a equal 3.52 # This should be changed to match the material
variable Tmelt equal 2000 # Some high temperature that is guaranteed to create a melt, but not so high as to cause problems. This is material dependent
variable melt_est equal 1300 # estimate of melting temperature (used for initial velocity assignment)

# create the simulation cell
lattice bcc $a
region 1 block 0 10 0 5 0 5 units lattice # longer in the x dimension to accomodate two phases
region solid block 0 5 0 5 0 5 units lattice
region liquid block 5 10 0 5 0 5 units lattice
create_box 1 1
create_atoms 1 region liquid

mass 1 238.02891

# create the different groups
group liquid region liquid

velocity all create ${melt_est} 100 mom yes dist gaussian loop all

# set the timestep
timestep 0.002 # Model dependent
