# Note: this script can be modified for different pair styles
# see melting.in for more info

# choose potential
pair_style adp
pair_coeff * * U_Mo.alloy.corrected.adp.txt U

# setup neighbor style
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes

# setup output
thermo 100
thermo_style custom step temp pe
thermo_modify norm no
