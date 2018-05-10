#! /usr/bin/env python

# This follows the example set in http://atsimpotentials.readthedocs.io/en/latest/potentials/eam_tabulation.html
import math
from atsim.potentials import EAMPotential
from atsim.potentials import Potential

def embed(rho):
    return ERROR

def density(rij):
    beta=3.1878 # angstroms ^ -1
    return rij**8*(math.exp(-beta*rij)+2**11*math.exp(-2*beta*rij))

def pair_AuAu(rij):
    Dm=0.6569 #eV
    Rm=2.5472 # angstroms
    alpha_m=1.2032 # angstroms ^-1
    return Dm*(1-math.exp(-alpha_m*(rij-Rm)))**2-Dm

def main():
    # requires the element, the atomic number, and the atomic weight, along with
    # the embedding and density functions
    eamPotentials = [ EAMPotential("Au", 79, 196.96657, embed, density)]
    # Requires the element interaction (two arguments, i.e. 'Ni', 'Ni'), and the
    # pair potential function.
    pairPotentials = [ Potential('Au', 'Au', pair_AuAu) ]


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

    with open("Au1_Foiles.eam", 'wb') as outfile:
        writeFuncFL(
            nrho, drho,
            nr, dr,
            eamPotentials,
            pairPotentials,
            out= outfile,
            title='Foiles-Hoyt Ni')
if __name__ == "__main__":
    main()
