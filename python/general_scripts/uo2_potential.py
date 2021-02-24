from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import math

# constants

# The uo2 potential consists of three distinct parts: the coulomb, buckingham, and morse terms
# MAKE SURE THE UNITS WORK OUT
def coulomb(r, charge1, charge2):
    # Note that 1/4*pi*e0 is in units of
    temp = [14.400774945963 * charge1 * charge2 / x for x in r]
    return temp

def buckingham(r, type1, type2):
    if type1 == type2:
        if type1 == 1:
            Aij = 294.7593 #eV
            rho_ij = 0.327022 #angstroms
            Cij = 0.0 #ev Angstroms^6
        elif type1 == 2:
            Aij = 1633.666
            rho_ij = 0.327022
            Cij = 3.95063
    else:
        Aij = 693.9297
        rho_ij = 0.327022
        Cij = 0.0

    temp = [Aij * np.exp(-x/rho_ij) - Cij / x**6 for x in r]
    return temp

def morse(r, type1, type2):
    if type1 == type2:
        if type1 == 1:
            Dij = 0.0 #eV
            beta_ij = 0.0 # angstroms^-1
            rij0 = 0.0 # angstroms
        elif type1 == 2:
            Dij = 0.0
            beta_ij = 0.0
            rij0 = 0.0
    else:
        Dij = 0.57745
        beta_ij = 1.65
        rij0 = 2.369

    temp = [Dij * (np.exp(-2*beta_ij * (x - rij0)) - 2*np.exp(-beta_ij * (x - rij0))) for x in r]
    return temp

x = [i/1000 for i in range(1, 10000)]

valsUU = [(xx + yy + zz) for xx, yy, zz in zip(coulomb(x, 2.4, 2.4), buckingham(x, 1, 1), morse(x, 1, 1))]
valsOO = [(xx + yy + zz) for xx, yy, zz in zip(coulomb(x, -1.2, -1.2), buckingham(x, 2, 2), morse(x, 2, 2))]
valsUO = [(xx + yy + zz) for xx, yy, zz in zip(coulomb(x, 2.4, -1.2), buckingham(x, 1, 2), morse(x, 1, 2))]

x = [y / 5.453 for y in x]


plt.figure(1)
plt.plot(x, coulomb(x, 2.4, 2.4), 'r-', label="Coulomb")
plt.plot(x, buckingham(x, 1, 1), 'g-', label="Buckingham")
plt.plot(x, morse(x, 1, 1), 'b-', label="Morse")
plt.plot(x,valsUU, 'c-', label="U-U Interaction")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(2)
plt.plot(x, coulomb(x, -1.2, -1.2), 'r-', label="Coulomb")
plt.plot(x, buckingham(x, 2, 2), 'g-', label="Buckingham")
plt.plot(x, morse(x, 2, 2), 'b-', label="Morse")
plt.plot(x,valsOO, 'c-', label="O-O Interaction")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(3)
plt.plot(x, coulomb(x, 2.4, -1.2), 'r-', label="Coulomb")
plt.plot(x, buckingham(x, 1, 2), 'g-', label="Buckingham")
plt.plot(x, morse(x, 1, 2), 'b-', label="Morse")
plt.plot(x,valsUO, 'c-', label="U-O Interaction")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(4)
plt.plot(x, valsUU, 'r-', label="U-U Interaction")
plt.plot(x, valsOO, 'g-', label="O-O Interaction")
plt.plot(x, valsUO, 'b-', label="U-O Interaction")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(5)
plt.plot(x, coulomb(x, 2.4, 2.4), 'r-', label="U-U Coulomb")
plt.plot(x, coulomb(x, -1.2, -1.2), 'g-', label="O-O Coulomb")
plt.plot(x, coulomb(x, 2.4, -1.2), 'b-', label="U-O Coulomb")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(6)
plt.plot(x, buckingham(x, 1, 1), 'r-', label="U-U Buckingham")
plt.plot(x, buckingham(x, 2, 2), 'g-', label="O-O Buckingham")
plt.plot(x, buckingham(x, 1, 2), 'b-', label="U-O Buckingham")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.figure(7)
plt.plot(x, morse(x, 1, 1), 'r-', label="U-U Morse")
plt.plot(x, morse(x, 2, 2), 'g-', label="O-O Morse")
plt.plot(x, morse(x, 1, 2), 'b-', label="U-O Morse")
plt.plot(x, np.zeros(len(x)), 'k--')
plt.xlabel(r"Distance ($a_0$)")
plt.ylabel("Potential (eV)")
plt.ylim([-1000,1000])
plt.legend()

plt.show()
