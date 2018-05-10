#! /usr/bin/env python

# This follows the example set in http://atsimpotentials.readthedocs.io/en/latest/potentials/eam_tabulation.html
import math
from atsim.potentials import EAMPotential
from atsim.potentials import Potential

def embed(rho):
    a0 = 3.52 # lattice constant, angstroms
    eSub = 4.45 # sublimation energy, eV
    B = 1.128781 # eV/angstrom^3 (1.8037e12 erg/cm^3) Bulk modulus
    aStar = 1/(a0 - 1) / (eSub / (9 * B * a0**3 / 4))**(0.5)
    return -eSub*(1 + (rij * aStar)) * math.exp(rij * aStar)

def density(rij):
    return rij**8*(math.exp(-3.58321*rij)+2**11*math.exp(-2*3.58321*rij))

def pair_NiNi(rij):
    return 1.39664*(math.exp(-2*1.22848*(rij-2.14146))-2*math.exp(-1.22848*(rij-2.14146)))

def main():
    # Requires the element interaction (two arguments, i.e. 'Ni', 'Ni'), and the
    # pair potential function.
    pairPotentials = [ Potential('Ni', 'Ni', pair_NiNi) ]
    # requires the element, the atomic number, and the atomic weight, along with
    # the embedding and density functions
    eamPotentials = [ EAMPotential("Ni", 28, 58.6934, embed, density)]


    # The number of density values to tabulate
    nrho = 50001
    drho = 0.001 # step between each density

    # The above two parameters determine the maximum density of the embedding
    # function as (nrho - 1)*drho.

    # The number of r values to tabulate
    nr = 4851
    dr = 0.001

    # Similarly here, this specifies the maximum distance between atoms as (nr - 1) * dr
    # The cutoff distance here is 4.85 Angstroms, so (4851 - 1) * dr = 4.850

    from atsim.potentials import writeFuncFL

    with open("Ni_Foiles_Hoyt.eam", 'wb') as outfile:
        writeFuncFL(
            nrho, drho,
            nr, dr,
            eamPotentials,
            pairPotentials,
            out= outfile,
            title='Foiles-Hoyt Ni')
if __name__ == "__main__":
    main()
