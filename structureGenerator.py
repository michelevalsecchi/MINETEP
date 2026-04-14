import random
import itertools

from math import asin, cos, pi, sin

import numpy as np

from matplotlib import pyplot as plt

class crystal:

    def __init__(self, natoms, atomicradius=1.11, cclass='tri'):
        """
        Initialize crystal class
        --------------------------------
        Parameters:
            natoms : int
                number of atoms in the unit cell
            atomic radius : float
                radius of atom in Angstrom
                default = 1.11 (Silicon)
            cclass: string 
                crystal class
                default = 't' (triclinic)
        --------------------------------
        returns:
            none
        """

        self.natoms = natoms
        self.radius = atomicradius
        self.cclass = cclass

        self.params = np.array([]) # initialise parameters. for triclinic cell this is a,b,c,theta,phi,psi
        self.cell_vec = np.array([[1.0,0,0],[0,1.0,0],[0,0,1.0]]) # initialise cell vectors 
        self.volume = 0

        # initialise grid
        self.dist = np.zeros((10,10,10))
        self.grid = np.zeros((10,10,10,3))
        self.atoms = np.zeros((10,10,10))

    def triclinicParams(self):
        """
        Generate a triclinic crystal i.e. no initial symmetry is assumed.
        The parameters of the crystal are described in terms of the conventional cell
        lengths a,b,c and angles theta, psi and phi where theta is the angle between
        unit cell vectors a and b (equivalent to the conventional theta), psi is the
        angle between a and c, and phi is defined in a diagram in the readme file.

        The unit cell volume is fixed as 1.6*number of atoms*volume of atom.
        Further constraints are specified on the unit cell lengths, assuming at least
        one atom must be able to fit in each unit cell vector direction, and a maximum
        of 2*the number of atoms are able to fit in each unit cell vector direction.

        Constraints are set on the angles theta and psi based on the constraints on the
        unit cell lengths specified:
        V/absin(theta) > 2*2*atomic_radius
        V/bcsin(psi) > 2*2*atomic radius
        Such that:
        sin(psi)>4r/c
        sin(theta)>4r/b

        The final angle phi is found from the constraints specified.
        Atoms are placed in the described unit cell using periodic boundary conditions,
        enforcing that the first atom lies at the origin of the unit cell.
        -------------------------------------------------------------------------------
        Parameters:
            none
        -------------------------------------------------------------------------------
        Returns:
            none
        """
        
        # randomly choose a,b,c with the constraints that
        # 4*radius< a,b,c< 4*radius*natoms

        min_dim= self.radius*4
        max_dim= self.radius*6      

        sin_phi = 1.1 # initialise to forbidden value to enter while loop 

        while sin_phi > 1. or sin_phi < -1.:
            self.params= np.random.uniform(low=min_dim, high=max_dim, size=3) #a,b,c

            # randomly choose theta and psi subject to the constraints
            # sin(psi)>4r/c, sin(theta)>4r/b
            # also ensure sin != 0 (can't have zero angles)
            theta_min= asin(min_dim/self.params[1])
            psi_min= asin(min_dim/self.params[2])
            self.params = np.append(self.params, (theta_min + (pi-theta_min)*np.random.rand(1)))
            self.params = np.append(self.params, (psi_min + (pi-psi_min)*np.random.rand(1)))  

            #fix volume of unit cell
            self.volume =1.6*self.natoms*(4/3)*pi*self.radius**3   

            # generate third angle based on volume constraint 
            sin_phi = self.volume/(self.params[0]*self.params[1]*self.params[2]*sin(self.params[3])*sin(self.params[4]))
        self.params = np.append(self.params, asin(sin_phi))

        # check the generated angles lie between 0 and pi
        assert self.params[3]<=pi and self.params[3]>0, "Theta angle is not between zero and pi"
        assert self.params[4]<=pi and self.params[4]>0, "Psi angle is not between zero and pi"

        # check phi lies between 0 and 2pi
        assert self.params[5]<=2*pi and self.params[5]>0, "phi angle is not above zero and less than 2pi" 
        assert (np.dot(self.cell_vec[0], np.cross(self.cell_vec[1], self.cell_vec[2])))- self.volume< 1.0e-9, "Volumes do not match" 

        self.cell_vec[0]=np.array([self.params[0], 0, 0])
        self.cell_vec[1]=np.array(self.params[1]*np.array([cos(self.params[3]),sin(self.params[3]),0]))
        self.cell_vec[2]= self.params[2]*np.array([cos(self.params[4]), sin(self.params[4])*cos(self.params[5]), sin(self.params[4])*sin(self.params[5])]) 

    def generateCellParams(self):
        """
        Generate the cell parameters depending on the crystal class selected
        ---------------------------------
        Parameters:
            none
        ---------------------------------
        Returns:
            none
        """
        if self.cclass=='tri':
            self.triclinicParams()
            self.volume =1.6*self.natoms*(4/3)*pi*self.radius**3   
        elif self.cclass == 'o':
            min_dim= self.radius*2
            max_dim= self.radius*6      
            self.volume =2*self.natoms*(4/3)*pi*self.radius**3   
            self.params= np.random.uniform(low=min_dim, high=max_dim, size=2) #a,b
            self.params = np.append(self.params,
                                    self.volume/self.params[0]/self.params[1])
            self.cell_vec[0] = self.params[0]*np.array([1.0,0.0,0.0])
            self.cell_vec[1] = self.params[1]*np.array([0.0,1.0,0.0])
            self.cell_vec[2] = self.params[2]*np.array([0.0,0.0,1.0])

        elif self.cclass == 't':
            min_dim= self.radius*2
            max_dim= self.radius*6      
            self.volume =2*self.natoms*(4/3)*pi*self.radius**3   
            self.params= np.random.uniform(low=min_dim, high=max_dim, size=1) #a,b
            self.params = np.append(self.params,
                                    self.volume/self.params[0]/self.params[0])
            self.cell_vec[0] = self.params[1]*np.array([1.0,0.0,0.0])
            self.cell_vec[1] = self.params[0]*np.array([0,1.0,0.0])
            self.cell_vec[2] = self.params[0]*np.array([0.0,0.0,1.0])

        elif self.cclass == 'c':
            min_dim= self.radius*2
            max_dim= self.radius*6      
            self.volume =2*self.natoms*(4/3)*pi*self.radius**3   
            self.params= np.array([self.volume**(1/3)])
            self.cell_vec[0] = self.params[0]*np.array([1.0,0.0,0.0])
            self.cell_vec[1] = self.params[0]*np.array([0.0,1.0,0.0])
            self.cell_vec[2] = self.params[0]*np.array([0.0,0.0,1.0])

        #initialize probability distribution of atoms
        self.dist = self.dist +  1/self.volume

        #create grid
        self.grid = np.zeros((10,10,10,3))
        for i in range(self.grid.shape[0]):
            for j in range(self.grid.shape[1]):
                for k in range(self.grid.shape[2]):
                    self.grid[i,j,k,:] += self.cell_vec[0]*i/10
                    self.grid[i,j,k,:] += self.cell_vec[1]*j/10
                    self.grid[i,j,k,:] += self.cell_vec[2]*k/10


    def add_atoms(self):
        """
        Add atoms at random positions in unit cell
        --------------------------------
        Parameters:
            none
        --------------------------------
        returns:
            none
        """
        choice_array = np.arange(1000)

        x = [1,-1,0] #generate translation vectors to neighbouring unit cells
        c = [p for p in itertools.product(x, repeat=3)]
        c = np.array(c)

        if self.dist.sum() > 0:


            #pick random positon from uniform distribution over unit cell
            if np.sum(self.atoms) == 0:
                itemindex = [0,0,0]
            else:
                temp = self.dist.flatten()/self.dist.sum()
                choice = np.random.choice(choice_array, p=temp)
                choice_array = choice_array.reshape((10,10,10))
                itemindex = np.where(choice_array==choice)

            self.atoms[itemindex[0],itemindex[1],itemindex[2]] = 1.0

            for i, vec in enumerate(c):
                #adjust distribution to prevent overlap between atoms
                #atoms in neighbouring unit cells are also considered

                distances = np.zeros((10,10,10,3))
                translation = vec[0]*self.cell_vec[0]+vec[1]*self.cell_vec[1]
                translation += vec[2]*self.cell_vec[2]
                distances[:,:,:,0] = -self.grid[:,:,:,0]+self.grid[itemindex[0],itemindex[1],itemindex[2],0]
                distances[:,:,:,0] += translation[0]
                distances[:,:,:,1] = -self.grid[:,:,:,1]+self.grid[itemindex[0],itemindex[1],itemindex[2],1]
                distances[:,:,:,1] += translation[1]
                distances[:,:,:,2] = -self.grid[:,:,:,2]+self.grid[itemindex[0],itemindex[1],itemindex[2],2]
                distances[:,:,:,2] += translation[2]
                distances = distances*distances
                distances = distances.sum(axis=3)
                distances = np.sqrt(distances)
                self.dist[(distances<=1.8*self.radius)] = 0.0
        else:
            raise Exception('not enough space, generate a new structure')

    def saveToFile(self, filename):
        """
        Saves crystal structure to .struct file
        -----------------------------------------
        Parameters:
            filename : string
                name of output file
        ----------------------------------------
        returns:
            None
        """
        array = np.array(np.where(self.atoms==1))
        with open(filename, "w") as file:
            if self.cclass == 'c':
                file.write("a: " + str(self.cell_vec[0])+ "\n")
                file.write("b: " + str(self.cell_vec[1]) + "\n")
                file.write("c: " + str(self.cell_vec[2]) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")
            elif self.cclass == 't':
                file.write("a: " + str(self.cell_vec[0])+ "\n")
                file.write("b: " + str(self.cell_vec[1]) + "\n")
                file.write("c: " + str(self.cell_vec[2]) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")
            elif self.cclass == 'o':
                file.write("a: " + str(self.cell_vec[0])+ "\n")
                file.write("b: " + str(self.cell_vec[1]) + "\n")
                file.write("c: " + str(self.cell_vec[2]) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")
            elif self.cclass == 'tri':    
                file.write("a: " + str(self.cell_vec[0])+ "\n")
                file.write("b: " + str(self.cell_vec[1]) + "\n")
                file.write("c: " + str(self.cell_vec[2]) + "\n")
                file.write("theta: " + str(self.params[3]) + "\n")
                file.write("psi: " + str(self.params[4]) + "\n")
                file.write("phi: " + str(self.params[5]) + "\n")
                file.write("-----------------------------" + "\n")
            for i in range(int(np.sum(self.atoms))):
                atom = self.grid[array[0,i],array[1,i],array[2,i],:]
                file.write("Si: "+ str(atom) + "\n")

    def add_all_atoms(self):
        try:
            for i in range(self.natoms):
                self.add_atoms()
        except:
            self.generateCellParams()
            self.atoms = np.zeros((10, 10, 10))
            self.dist *= 0.0
            self.dist +=  1/self.volume
            self.add_all_atoms()
