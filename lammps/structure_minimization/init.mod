units metal
atom_style charge
boundary p p p

variable Tanneal equal 1200 # annealing temperature

read_data myfile # use the desired structure to minimize

mass 1 140.1160 # Ce
mass 2 15.9994 # O

set type 1 charge +4
set type 2 charge -2

timestep 0.002
