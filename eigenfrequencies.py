#!/usr/bin/env python3
"""
This program computes the eigenfrequencies of the dynamical matrix
in order to determine the dynamical stability of computed structures.
It requires a structure to be passed to the code.
If there are any imaginary eigenfrequencies, a structure is dictated
unstable.
"""

import numpy as np
import itertools
from scipy.optimize import minimize
from scipy.linalg import eigvalsh as eig

import potential as tpot
from hessian import hessian

def load_file(filename):
	"""
	Load crystal structure given filename, provided it follows
	the format of the .struct files as detailed in the README
	---------------------------------------------------------
	Parameters:
		filename: string
		name of input file
	Returns:
		param: lattice parameters (a, b, c, θ, φ, ψ)
		atoms: positions of 3N atoms in the unit cell
	"""

	with open(filename, 'r') as file:
		lines= file.readlines()

	lines = [j.strip('\n') for j in lines]
	lines = [j.split(': ') for j in lines]
	del(lines[6])

	data = np.array([])

	for i in range(len(lines)):
		A= lines[i][1]
		if i not in (3,4,5):
			A=A[1: -1]
		A= A.split()
		A= np.array(A, dtype=float)
		data= np.append(data, A)

	# lattice params
	param = np.zeros(6)

	for j in range(3):
		param[j]= np.sqrt(np.sum(np.square(data[3*j:3*j+3])))
		param[j + 3] = data[j + 9]

	# atomic posns
	atoms = data[12:]

	return param, atoms

def supercell_gen(param,atoms):
    """
	This function generates the 2*2*2 supercell of the unit cell with specified
    parameters.

    Parameters
    ----------
    param: float, array-like
           lattice parameters a, b, c, θ, φ, ψ
    atoms: float, array-like
            positions of all atoms in the unit cell

    Returns
    -------
    new_posns: float, array-like
                     positions of all atoms in the supercell
    """

    x= np.array([1, 0]) # translation vectors to neighbouring cells
    c = [p for p in itertools.product(x, repeat=3)]
    c = np.array(c) # includes itself
    tempc = np.copy(c)
    c[0,:] = tempc[-1,:]
    c[-1,:] = tempc[0,:]
    unitCellVec=[param[0], param[1], param[2]]
    c= c*unitCellVec

    posns= [atoms[i:i+3] for i in range(0, len(atoms), 3)]
    posns= np.array(posns)
    new_posns= []

    for translation in c:
        new_pos= posns + translation
        new_posns.append(new_pos)

    new_posns= np.array(new_posns)
    new_posns= new_posns.flatten()
    param[0:3]= 2*param[0:3]
    return param, new_posns

def stability(param,atoms):
	"""
	This function computes the natural frequencies
	using the dynamical matrix. If there are any
	imaginary frequencies, the structure is deemed
	unstable.

	Parameters
	----------
	param: float, array-like
			 Lattice parameters a, b, c, θ, φ, ψ
	atoms: float, array-like
             Positions of all atoms

	Returns
	-------
	None
	"""


	eps = 1e-3#eps= np.finfo(np.float32).eps
	m= 4.6637066e-26 # mass of silicon in kg (=J m^-2 s^2)
	m= m*6.242e+18 # J to eV
	m= m*1.0e-20 # length to Angstrom

	H= hessian.hess(tpot.wrap_pot, atoms, param, eps)
	dynmat= H/m
	eigvals= eig(dynmat)
	w= np.sqrt(eigvals+0j)
	w_real= w.real
	w_im= w.imag
	print(w*1e-12/(2*np.pi)) #Show frequencies in THz
	print('Frequencies are in THz.')
	if any(w.imag > 1e-9):
		print('Imaginary frequencies. Structure is unstable')
	else:
		print('No imaginary frequencies. Structure is stable')
