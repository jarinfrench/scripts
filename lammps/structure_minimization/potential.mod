# Note: this script should be modified for different pair styles

# choose potential
kspace_style pppm 1.0e-5
pair_style buck/coul/long
pair_coeff * * 0.0 1.0 0.0
pair_coeff 1 2 755.1311 0.429 0.0 # Ce-O
pair_coeff 2 2 9533.421 0.234 224.88 # O-O

# setup neighbor style
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes

# setup minimization style
min_style cg

# setup output
thermo 100
thermo_style custom step temp pe lx ly lz vol
thermo_modify norm no flush yes
