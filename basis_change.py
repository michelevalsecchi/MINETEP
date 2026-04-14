"""This module contains two functions that allow to change the coordinate system
from cartesian to lattice vector representation and vice versa"""

import numpy as np

def c_l(var):
    """
    This function performs a linear transformation from the cartesian basis to the
    Bravais lattice vector basis.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        Array of the (3N - 3) atomic coordinates of all the unit cell atoms except
        the one at the origin and the 6 lattice parameters (a, b, c, θ, φ, ψ) in
        the Cartesian basis.

    Returns
    -------
    var2 : array-like (3N + 3, 1)
        Array of the (3N - 3) atomic coordinates of all the unit cell atoms except
        the one at the origin and the 6 lattice parameters (a, b, c, θ, φ, ψ) in
        the Bravais lattice basis.

    """
    #Number of atoms
    N = int((len(var)-3)/3)

    #Extracting positions and lattice parameters in the cartesian basis
    cpos = var[:-6]
    lpar = var[-6:]
    cpos = cpos.reshape(N-1,3)

    #Defining the Bravais lattice vectors in the Cartesian basis
    a = lpar[0]*np.array([1,0,0])
    b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
    c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])

    #Performing the transformation
    trans = np.array([a,b,c])
    lpos = np.dot(cpos,np.linalg.inv(trans))
    var2 = np.append(lpos.flatten(),lpar)

    return var2

def l_c(var):
    """
    This function performs a linear transformation from the cartesian basis to the
    Bravais lattice vector basis.

    Parameters
    ----------
    var : array-like (3N + 3, 1)
        Array of the (3N - 3) atomic coordinates of all the unit cell atoms except
        the one at the origin and the 6 lattice parameters (a, b, c, θ, φ, ψ) in
        the Bravais lattice basis

    Returns
    -------
    var2 : array-like (3N + 3, 1)
        Array of the (3N - 3) atomic coordinates of all the unit cell atoms except
        the one at the origin and the 6 lattice parameters (a, b, c, θ, φ, ψ) in
        the Cartesian basis.

    """
    #Number of atoms
    N = int((len(var)-3)/3)

    #Extracting positions and lattice parameters in the lattice basis
    lpos = var[:-6]
    lpar = var[-6:]
    lpos = lpos.reshape(N-1,3)

    #Defining the Bravais lattice vectors in the Cartesian basis
    a = lpar[0]*np.array([1,0,0])
    b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
    c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])

    #Performing the transformation
    trans = np.array([a,b,c])
    cpos = np.dot(lpos,trans)
    var2 = np.append(cpos.flatten(), lpar)

    return var2
