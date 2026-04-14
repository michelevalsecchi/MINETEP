#!/usr/bin/env python3
import numpy as np
import argparse
import structureGenerator as sg
import eigenfrequencies as ef
import optimizer as opt
import E_V as ev
import matplotlib.pyplot as plt

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--path', type=str, help="""path to directory to save structure files to
                        or load input files from; default: current directory.
                        For optimisation or stability testing, or extraction of
                        thermodynamic data, specify path to
                        input file (including the filename).""", default=".",
                        nargs='*')

    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('--mode', type=str, help="""which mode do you want to
                        select? s for structure generation, o for structure
                        optimisation, t for calculating thermodynamic
                        properties, p for finding eigenfrequencies""")

    parser.add_argument('-n', type=int, help="""number of atoms; only
                            required for structure generation""", default=0)

    parser.add_argument('--structure', type=str, help="""specify Bravais
                        lattice to be generated; default triclinic,'t'
                        tetragonal, 'o' orthorhombic, 'c' cubic. Default is
                        'tri'. If the symmetry of the structure should be
                        maintained during optimisation, this argument is
                        required. Otherwise all lattice parameters will be
                        changed.""",default="tri")

    parser.add_argument('-k', type=int, help="""number of structures to be
                        generated; only required for structure generation, default=1""", default=1)

    parser.add_argument('-r', type=float, help="""atomic radius; only
                            required for structure generation; default=1.11""", default=1.11)

    parser.add_argument('-bsh', type=float, help="""type of minimization: 1 for
                        basinhopping, 0 for simple L-BFGS-B. Required for structure
                        optimization.""",default=1)

    parser.add_argument('-T', type=float, help="""temperature of the basinhopping
                        algorithm. Required if bsh is 1""", default=0.35)

    parser.add_argument('-it', type=float, help="""number of iterations of the
                        basinhopping algorithm. Required if bsh is 1""", default=100)

    parser.add_argument('-step', type=float, help="""stepsize of the basinhopping
                        algorithm. Required if bsh is 1""", default=0.15)

    parser.add_argument('-scell', type=str, help="""Specify whether to use the supercell or unit cell
                        in computing the eigenfrequencies. s for supercell and u for unit cell.""", default='u')

    args = parser.parse_args()
    mode = args.mode

    print("""
     ▄▄       ▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄
    ▐░░▌     ▐░░▌▐░░░░░░░░░░░▌▐░░▌      ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
    ▐░▌░▌   ▐░▐░▌ ▀▀▀▀█░█▀▀▀▀ ▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀  ▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌
    ▐░▌▐░▌ ▐░▌▐░▌     ▐░▌     ▐░▌▐░▌    ▐░▌▐░▌               ▐░▌     ▐░▌          ▐░▌       ▐░▌
    ▐░▌ ▐░▐░▌ ▐░▌     ▐░▌     ▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄      ▐░▌     ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌
    ▐░▌  ▐░▌  ▐░▌     ▐░▌     ▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
    ▐░▌   ▀   ▐░▌     ▐░▌     ▐░▌   ▐░▌ ▐░▌▐░█▀▀▀▀▀▀▀▀▀      ▐░▌     ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌       ▐░▌     ▐░▌     ▐░▌    ▐░▌▐░▌▐░▌               ▐░▌     ▐░▌          ▐░▌
    ▐░▌       ▐░▌ ▄▄▄▄█░█▄▄▄▄ ▐░▌     ▐░▐░▌▐░█▄▄▄▄▄▄▄▄▄      ▐░▌     ▐░█▄▄▄▄▄▄▄▄▄ ▐░▌
    ▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░▌      ▐░░▌▐░░░░░░░░░░░▌     ▐░▌     ▐░░░░░░░░░░░▌▐░▌
     ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀        ▀▀  ▀▀▀▀▀▀▀▀▀▀▀       ▀       ▀▀▀▀▀▀▀▀▀▀▀  ▀
              """)

    print("===========================================================================")

    print("""Welcome to the MINETEP (= 'MINETEP is not an electronic structure
          total energy package') software package for crystal structure
          generation, structure optimisation, structure stability testing and
          calculation of thermodynamic properties""")

    print("===========================================================================")

    parser.print_help()

    if mode == "s":

        nstruc = args.k
        natoms = args.n
        r = args.r
        structure = args.structure

        for i in range(nstruc):

            if structure == "c":

                struc = sg.crystal(natoms, r,'c')
                struc.generateCellParams()
                struc.add_all_atoms()
                struc.saveToFile(args.path[0]+"/"+"input"+"Si"+str(natoms)+"_"+str(i)
                                 +"_cubic"+".struct")

            elif structure == "t":

                struc = sg.crystal(natoms, r,'t')
                struc.generateCellParams()
                struc.add_all_atoms()
                struc.saveToFile(args.path[0]+"/"+"input"+"Si"+str(natoms)+"_"+str(i)
                                 +"_tetragonal"+".struct")

            elif structure == "o":

                struc = sg.crystal(natoms, r,'o')
                struc.generateCellParams()
                struc.add_all_atoms()
                struc.saveToFile(args.path[0]+"/"+"input"+"Si"+str(natoms)+"_"+str(i)
                                 +"_orthorhombic"+".struct")
            else:

                struc = sg.crystal(natoms, r, 'tri')
                struc.generateCellParams()
                struc.add_all_atoms()
                struc.saveToFile(args.path[0]+"/"+"input"+"Si"+str(natoms)+"_"+str(i)
                                 +"_triclinic"+".struct")

        print("\n")
        print("\n")
        print("##############################################################################")
        print("##############################################################################")
        print("##############################################################################")
        print("##############################################################################")
        print("\n")
        print("\n")
        print("Structure generation successful. Structure files are stored "+
              args.path[0]+"/")

    if mode== "p":

        param, atoms= ef.load_file(args.path[0])

        if args.scell=='s':

            param,atoms= ef.supercell_gen(param, atoms) # include supercell

        ef.stability(param, atoms) # prints whether or not a structure is stable

    if mode== "t":
        string = np.empty([])
        for i in range(len(args.path)):
            string = np.append(string,"Structure %d"%(i+1))
        string = tuple(string[1:])
        for i, j in enumerate(args.path):
            vrange,E,P,H,b = ev.eh_plot(j)
            plt.figure(1) #Energy vs volume plot
            plt.plot(vrange,E)
            plt.xlabel('Volume per atom (Bohr^3)')
            plt.ylabel('Energy per atom (Ry)')
            plt.grid(True)
            print('Bulk modulus of structure %d: %4f GPa'%((i+1),b))
        string = list(string)
        plt.legend(string,loc='best')
        plt.savefig("Energy vs Pressure")
        for i, j in enumerate(args.path):
            vrange,E,P,H,b = ev.eh_plot(j)
            plt.figure(2) #Enthalpy vs pressure plot
            plt.plot(P,H)
            plt.xlim(0,20)
            plt.xlabel('Pressure (GPa)')
            plt.ylabel('Enthalpy per atom (eV)')
            plt.grid(True)
        plt.legend(string,loc='best')
        plt.savefig('Enthalpy vs Pressure')


    if mode == "o": # optimise structure based on specified crystal system

        if args.structure == 'c':

            opt_var = opt.cub_opt(args.path[0],args.r,bool(args.bsh),args.T,args.it,args.step)
            x0 = np.array([0,0,0])
            atoms = np.append(x0,opt_var[:-6])
            lpar = opt_var[-6:]
            a = lpar[0]*np.array([1,0,0])
            b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
            c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])

            with open("opt"+"Si"+str(args.n)+"_"+'0'
                             +"_cubic"+".struct", "w+") as file:

                file.write("a: " + str(a)+ "\n")
                file.write("b: " + str(b) + "\n")
                file.write("c: " + str(c) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")

                for i in range(args.n):
                    file.write("Si: "+ str(atoms[3*i:3*i+3]) + "\n")


        if args.structure == 'o':

            opt_var = opt.orth_opt(args.path[0],args.r,bool(args.bsh),args.T,args.it,args.step)
            x0 = np.array([0,0,0])
            atoms = np.append(x0,opt_var[:-6])
            lpar = opt_var[-6:]
            a = lpar[0]*np.array([1,0,0])
            b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
            c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])


            with open("opt"+"Si"+str(args.n)+"_"+'0'
                             +"_orthorhombic"+".struct", "w+") as file:

                file.write("a: " + str(a)+ "\n")
                file.write("b: " + str(b) + "\n")
                file.write("c: " + str(c) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")

                for i in range(args.n):
                    file.write("Si: "+ str(atoms[3*i:3*i+3]) + "\n")


        if args.structure == 't':

            opt_var = opt.tet_opt(args.path[0],args.r,bool(args.bsh),args.T,args.it,args.step)
            x0 = np.array([0,0,0])
            atoms = np.append(x0,opt_var[:-6])
            lpar = opt_var[-6:]
            a = lpar[0]*np.array([1,0,0])
            b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
            c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])

            with open("opt"+"Si"+str(args.n)+"_"+'0'
                             +"_tetragonal"+".struct", "w+") as file:

                file.write("a: " + str(a)+ "\n")
                file.write("b: " + str(b) + "\n")
                file.write("c: " + str(c) + "\n")
                file.write("theta: " + str(np.pi/2)+ "\n")
                file.write("psi: " + str(np.pi/2) + "\n")
                file.write("phi: " + str(np.pi/2) + "\n")
                file.write("-----------------------------" + "\n")
                for i in range(args.n):
                    file.write("Si: "+ str(atoms[3*i:3*i+3]) + "\n")

        if args.structure == 'tri':

            opt_var = opt.tri_opt(args.path[0],args.r,bool(args.bsh),args.T,args.it,args.step)
            x0 = np.array([0,0,0])
            atoms = np.append(x0,opt_var[:-6])
            lpar = opt_var[-6:]
            a = lpar[0]*np.array([1,0,0])
            b = lpar[1]*np.array([np.cos(lpar[3]),np.sin(lpar[3]),0])
            c = lpar[2]*np.array([np.cos(lpar[4]),np.sin(lpar[4])*np.cos(lpar[5]),np.sin(lpar[4])*np.sin(lpar[5])])


            with open("opt"+"Si"+str(args.n)+"_"+'0'
                             +"_triclinic"+".struct", "w+") as file:

                file.write("a: " + str(a)+ "\n")
                file.write("b: " + str(b) + "\n")
                file.write("c: " + str(c) + "\n")
                file.write("theta: " + str(lpar[3])+ "\n")
                file.write("psi: " + str(lpar[4]) + "\n")
                file.write("phi: " + str(lpar[5]) + "\n")
                file.write("-----------------------------" + "\n")

                for i in range(int(atoms.size/3)):
                    file.write("Si: "+ str(atoms[3*i:3*i+3]) + "\n")

if __name__ == "__main__":
    main()
