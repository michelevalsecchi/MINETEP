"""This module contains a function to plot energy vs volume and enthalpy vs Pressure
    of an optimized crystal structure contained in a file with extension .struct;
    the function also calculates and prints the bulk modulus of that specific
    crystal structure"""

import numpy as np
import scipy.misc
import matplotlib.pyplot as plt
from tersoff_pot import potential as pot
from optimizer import var_unpack

def eh_plot(filename,vstep=0.1):
    """
    Calculates and plots the energy vs volume and enthalpy vs pressure
    of a specific crystal structure.
    ---------------------------------
    Parameters:
    filename: str;
        Name of the file from which the function reads the atomic positions and
        lattice parameters
    vstep: float;
        Step for range.
    """

    #Initializing arrays
    x0, var = var_unpack(filename)
    vmin = var[-6]*var[-5]*var[-4]*np.sin(var[-3])*np.sin(var[-2])*np.sin(var[-1])
    vrange = np.arange(vmin-round(0.5*vmin),vmin+round(0.5*vmin) + vstep,vstep)
    nparam = np.zeros((len(var),len(vrange)))
    E = np.zeros(len(vrange))
    P = np.array([])
    H = np.array([])

    #Create the matrix nparam whose columns are the modified atomic positions
    #and lattice parameters
    for i in range(len(vrange)):
        lrat = np.cbrt(vrange[i]/vmin)
        nparam[:-6,i] = var[:-6]*lrat
        nparam[-6:-3,i] = var[-6:-3]*lrat
        nparam[-3:,i] = var[-3:]

    #Create Energy and Enthalpy arrays
    for j in range(len(vrange)):
        E[j] = pot.u_pot(nparam[:,j],x0)
        if vrange[j] <= vmin:
            P = np.append(P,-(pot.u_pot(nparam[:,j+2],x0)-pot.u_pot(nparam[:,j],x0))/(2*vstep))
            H = np.append(H,E[j]+vrange[j]*P[j])

    #Estimate the bulk modulus of the relaxed structure
    eps = 1e-3
    uparam=np.zeros(len(var))
    dparam=np.zeros(len(var))
    vu = vmin + eps
    vd = vmin - eps
    lratu = np.cbrt(vu/vmin)
    lratd = np.cbrt(vd/vmin)
    uparam[:-6] = var[:-6]*lratu
    uparam[-6:-3] = var[-6:-3]*lratu
    uparam[-3:] = var[-3:]
    dparam[:-6] = var[:-6]*lratd
    dparam[-6:-3] = var[-6:-3]*lratd
    dparam[-3:] = var[-3:]
    B_0 = vmin*(pot.u_pot(uparam,x0)+pot.u_pot(dparam,x0)-2*pot.u_pot(var,x0))/(eps**2)


    J_eV = 1.602*1e-19 #Defining constant
    m_Å = 1e-10
    Å_Bohr = 0.5292
    eV_Rydberg = 13.6057
    N = int((len(var) - 3)/3)
    E = E/(N*eV_Rydberg) #Energy per atom in Ry
    H = H/N #Enthalpy per atom in eV
    vrange = vrange/(N*(Å_Bohr**3)) #Volume per atom in Bohr^3
    P = P*1e-9*J_eV/(m_Å**3) #Pressure in GPa
    B_0 = np.round(B_0*1e-9*J_eV/(m_Å**3),decimals=2) #Bulk modulus in GPa


    return vrange,E,P,H,B_0
