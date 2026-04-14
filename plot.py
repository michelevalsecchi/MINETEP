#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import basis_change
import optimizer

def plot(filename):
    """
    This function plots the crystal structure included in the file named
    filename.

    Parameters
    ----------
    filename : string
        The file containing the information about the crystal structure. The
        file has to be a .struct file
    """

    x0,var = optimizer.var_unpack(filename)
    args = np.append(x0, var[:-6])
    param = var[-6:] # a, b, c, theta, phi, psi
    a = param[0]*np.array([1,0,0])
    b = param[1]*np.array([np.cos(param[3]),np.sin(param[3]),0])
    c = param[2]*np.array([np.cos(param[4]),np.sin(param[4])*np.cos(param[5]),np.sin(param[4])*np.sin(param[5])])
    posns= [args[i:i+3] for i in range(0, len(args), 3)]
    posns= np.array(posns)
    xs= [item[0] for item in posns]
    ys= [item[1] for item in posns]
    zs= [item[2] for item in posns]
    fig= plt.figure()
    ax= fig.add_subplot(111, projection='3d')
    x,y,z= zip(x0,x0,x0)
    u,v,w= zip(a,b,c)
    ax.quiver(x,y,z,u,v,w,arrow_length_ratio=0.01)
    ax.scatter(xs, ys, zs, s=800)
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis')
    plt.xlim((-4,4))
    plt.ylim((-4,4))
    plt.axis('square')
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, help="""file to read the crystal
                        structure from. It has to be a .struct file""")
    args = parser.parse_args()
    plot(args.filename)
