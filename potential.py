#!/usr/bin/env python3

"""
This program calculates the potential energy of silicon atoms in a unit cell
using a Tersoff potential given the atomic coordinates and unit cell geometry.
It calls the Fortran module "tersoff_pot.f90" which has been compiled using f2py.

The functional form of the Tersoff potential and the
empirical parameters used can be found at: J. Tersoff, Phys. Rev. B 38, 9902 (1988).
"""

import numpy as np
from tersoff_pot import potential as pot #fortran potential
import basis_change
import optimizer

#-------------------------------------------------------------------------------
def wrap_pot(x, l_param):
    """
    This function wraps the function that computes the total potential energy of the atoms in one unit
    cell interacting via a Tersoff Potential.

    Parameters
    ----------
    x : array-like (3N, 1)
        Array containing the positions of the N atoms in the unit cell, 1 < N < 13

    l_param : array-like (6, 1)

    Returns
    -------
    V : float
        The total potential energy of the atoms.
    """
    x0 = x[:3] #The position of the 1st atom in the unit cell.
    var = np.append(x[3:], l_param) #Positions of N-1 remaining atoms and lattice parameters

    #Calculate the number of atoms per unit cell, Nc
    Nc = int((len(var) - 3)/3)

    #Calculate potential
    V = pot.u_pot(var, x0, Nc)

    return V

def fpot(var, x0):
    """
    This function computes the total potential energy of the atoms in one unit
    cell interacting via a Tersoff Potential.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        The coordinates of the N-1 atoms in the unit cell in the basis
        of the lattice vectors and the 6 lattice parameters (see Documentation).
        We exclude the atom at the origin from this list given it is not avariable
        but a fixed parameter.

    x0 : array-like (3, 1)
        The position of the first atom at the origin of the unit cell.

    Returns
    -------
    V : float
        The total potential energy of the atoms.
    """
    #Calculate the number of atoms per unit cell, Nc
    Nc = int((len(var) - 3)/3)

    var2 = basis_change.l_c(var)

    V = pot.u_pot(var2, x0, Nc)

    return V

def cub_fpot(var, x0):
    """
    This function computes the total potential energy of the atoms in one cubic
    unit cell interacting via a Tersoff Potential.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        The coordinates of the N-1 atoms in the unit cell in the basis
        of the lattice vectors and the 6 lattice parameters (see Documentation).
        We exclude the atom at the origin from this list given it is not avariable
        but a fixed parameter.

    x0 : array-like (3, 1)
        The position of the first atom at the origin of the unit cell.

    Returns
    -------
    V : float
        The total potential energy of the atoms.
    """

    a = var[-1]
    var2 = np.append(var,(a,a,np.pi/2,np.pi/2,np.pi/2))
    V = fpot(var2,x0)

    return V

def orth_fpot(var, x0):
    """
    This function computes the total potential energy of the atoms in one cubic
    unit cell interacting via a Tersoff Potential.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        The coordinates of the N-1 atoms in the unit cell in the basis
        of the lattice vectors and the 6 lattice parameters (see Documentation).
        We exclude the atom at the origin from this list given it is not avariable
        but a fixed parameter.

    x0 : array-like (3, 1)
        The position of the first atom at the origin of the unit cell.

    Returns
    -------
    V : float
        The total potential energy of the atoms.
    """

    var2 = np.append(var,(np.pi/2,np.pi/2,np.pi/2))
    V = fpot(var2,x0)

    return V

def tet_fpot(var, x0):
    """
    This function computes the total potential energy of the atoms in one cubic
    unit cell interacting via a Tersoff Potential.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        The coordinates of the N-1 atoms in the unit cell in the basis
        of the lattice vectors and the 6 lattice parameters (see Documentation).
        We exclude the atom at the origin from this list given it is not avariable
        but a fixed parameter.

    x0 : array-like (3, 1)
        The position of the first atom at the origin of the unit cell.

    Returns
    -------
    V : float
        The total potential energy of the atoms.
    """

    a = var[-1]
    var2 = np.append(var,(a,np.pi/2,np.pi/2,np.pi/2))
    V = fpot(var2,x0)

    return V
