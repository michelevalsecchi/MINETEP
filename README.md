# MINETEP

MINETEP is a software package for the generation and optimization of Silicon crystal structures, interacting via a Tersoff potential. It has functionalities for testing the stability of the structures and the calculation of thermodynamic properties. It has been developed by Michele Valsecchi, Mario Zauchner, and Christopher Keegan, all part of the CDT in Theory and Simulation of Materials (Imperial College London, 2018/2019)

# Installation
To install the fortran modules, make sure that you have the latest version of numpy and use
```
f2py3 -c -m hessian hessian.f90
f2py3 -c -m tersoff_pot tersoff_pot.f90
```
On some machines, the command might be ```f2py```. There is a chance that f2py has to be called directly via numpy. In this case, it is convenient to define an alias.
For example:
```bash
alias f2py3='python3 -m numpy.f2py'
```
For more information, please check the f2py documentation.

Ensure that you have the latest versions of matplotlib, scipy and numpy installed.
# Usage

The package is used as follows:
```
usage: minetep.py [-h] [--path PATH] [--mode MODE] [-n N]
                  [--structure STRUCTURE] [-k K] [-r R] [-bsh BSH] [-T T]
                  [-it IT] [-step STEP] [-scell SCELL]

optional arguments:
  -h, --help            show this help message and exit
  --path PATH           path to directory to save structure files to or load
                        input files from; default: current directory. For
                        optimisation or stability testing, or extraction of
                        thermodynamic data, specify path to input file
                        (including the filename).
  -n N                  number of atoms; only required for structure
                        generation
  --mode MODE		s for structure generation, o for optimization
  --structure STRUCTURE
                        specify Bravais lattice to be generated; default
                        triclinic,'t' tetragonal, 'o' orthorhombic, 'c' cubic.
                        Default is 'tri'. If the symmetry of the structure
                        should be maintained during optimisation, this
                        argument is required. Otherwise all lattice parameters
                        will be changed.
  -k K                  number of structures to be generated; only required
                        for structure generation, default=1
  -r R                  atomic radius ; only required for structure
                        generation; default=1.11
  -bsh BSH              type of minimisation: 1 for basinhopping, 0 for simple
                        L-BFGS-B. Required for structure optimisation.
  -T T                  temperature of the basinhopping algorithm. Required if
                        bsh is 1
  -it IT                number of iterations of the basinhopping algorithm.
                        Required if bsh is 1
  -step STEP            stepsize of the basinhopping algorithm. Required if
                        bsh is 1
  -scell SCELL          Specify whether to use the supercell or unit cell in
                        computing the eigenfrequencies. s for supercell and u
                        for unit cell.

required named arguments:
  --mode MODE           which mode do you want to select? s for structure
                        generation, o for structure optimisation, t for
                        calculating thermodynamic properties, p for finding
                        eigenfrequencies
```

## Structure generation

MINETEP can generate random crystal structures based on a number of constraints.
The main constraint of the generation heuristic is the total volume of the unit cell (for triclinic crystals), which is given by the following expression:

```math
V = \frac{32}{15}\pi r^3N
```
and for orthogonal crystal systems:
```math
V = \frac{8}{3}\pi r^3N
```
where $`N`$ is the number of atoms in the unit cell and $`r`$ is the atomic radius (both are input parameters of the program). The choice of maximum volume was made such that the volume is large enough to allow non close-packed structures, while ensuring that atoms do not overlap.

To generate a structure, set the '--mode' flag to 's'.

### For example:

```bash
./minetep.py --mode s --path input_files -n 4 -k 1000 -r 2.45 --structure c
```

This command generates 1000 random crystal (cubic) structures, containing four atoms in each unit cell, each with an atomic radius of 2.45 Å.
The structure files are stored in the ./input_files directory.

The output files have the following structure:


```
a: [8.87442469   0.          0.        ]
b: [-2.92308285  3.80038398  0.        ]
c: [0.34954952   8.06635251  0.67943803]
theta: 2.22644410204132
psi: 1.5276417754077363
phi: 0.08403277531683041
-----------------------------
Si: [0. 0. 0.]
Si: [1.55248656 1.9933089  0.13588761]


```
where a, b and c are the unit cell vectors. Below the cell parameters, the coordinates of the Si atoms within the unit cell are shown.
Note that the coordinates are not given with respect to unit cell vectors, but are represented in an orthogonal cartesian coordinate system. The first atom is always placed at the origin of the coordinate system.


Diagram detailing the coordinate system used:
![alt text](diagram.png?raw=true "Image detailing coordinate system")

## Structure optimisation

MINETEP allows the user to minimise the total energy of a crystal of silicon with an arbitrary number of atoms in the unit cell.
This tool works by specifying the path of the starting guess (--path), the number of atoms in the unit cell (-n), the crystal system (--structure) and the desired optimisation method (-bsh 1 or 0, default 1).

The program uses by default the basinhopping algorithm (Wales, David J. 2003, Energy Landscapes, Cambridge University Press, Cambridge, UK) coupled with L-BFGS-B to semi-randomly jump between local minima of the potential energy surface (PES), storing the lowest value found. This method has proven to be very successful in systems with very high energy barriers between local minima and complicated PESs.
The basinhopping algorithms (imported from scipy.optimize.basinhopping - check there for further information) relies on three main parameters:

  - The number of jumps (-it) between local minima before stopping the search
  - The stepsize (-step), which is the absolute value of coordinate displacement at every jump. After every minimisation, the algorithm displaces every atomic lattice coordinate and every lattice vector  magnitude of the last accepted jump by a random number between -step and +step (-2*step and +2*step for lattice angles), and minimises the function using the displaced coordinates as starting guess
  - The temperature (-T), which provides a criterion to accept or reject a jump. A jump is always accepted if it lowers the energy; if the energy is raised, the jump is accepted with probability $`exp(-(E_{new}-E_{old})/T)`$, much like a Boltzmann factor. This ensures an effective circumvention of energy barriers between local minima.

The stepsize should be of the order of the separation between the coordinates of local minima in the function to be minimised and the temperature of the order of the difference in energy of these local minima; from our experience flexible and effective parameters were 0.15 and 0.35, respectively, which are the default values. As a rule of thumb, as the number of atoms increases temperature should be increased (for example to 0.4/0.5) and stepsize should be slightly increased (for example to 0.2). The first runs should have bigger temperatures and stepsizes to allow exploration of the PES, while further refinements should use a smaller stepsizes and temperature (e.g. 0.1 and 0.25 respectively).
During every run the user can witness the 'movement' of the algorithm along the PES and see what jumps are accepted and what value the energy takes. This information has to be watched carefully as it guides later parameter settings.
The number of jumps per run can be chosen according to the time available. Generally the time taken per minimisation scales as the cube of the number of atoms (because of the analytic form of the Tersoff potential), and one may want to lower the total number of jumps as the number of atoms increases. After an unsuccessful run, the user should decide based on the information gathered during the run whether to re-start from another random structure or continuing from the last minimum found, loading the opt_<original_file>.struct file.

Every run writes the energy minima found in a file named 'energy_temp.txt' that can be opened with the command

```bash
./optimizer.py
```
to produce a histogram with the number of times that a particular minimum has been found. This may provide insight on important local minima of structures with many minima. After moving on to different structures the user should delete the file 'energy_temp.txt' from the working directory.

### For example:

The command:

```bash
./minetep.py --mode o --path 'crystal.struct' -n 4 --structure c -bsh 1 -T 0.3 -it 150 -step 0.1
```
will minimise the energy of a cubic crystal structure with 4 atoms in the unit cell through the basin-hopping algorithm with parameters T = 0.3, stepsize = 0.1 and 150 total jumps. The initial guess is read from the file 'crystal.struct'.

The command:

```bash
./minetep.py --mode o --path 'crystal.struct' -n 4 --structure tri -bsh 0
```
will minimise the energy of a triclinic crystal structure with 4 atoms in the unit cell using a single minimisation with the L-BFGS-B algorithm using the structure stored in crystal.struct as a starting guess

## Eigenfrequencies
The eigenfrequencies of the dynamical matrix are computed in order to determine the stability of crystal structures. It requires a structure to be passed to the code. If there are any imaginary eigenfrequencies, the structure is deemed unstable.

The Hessian is computed using the centred difference algorithm, from which the dynamical matrix is constructed and the eigenvalues are computed.

To generate the eigenfrequencies, set the '--mode' flag to 'p', and make sure you are passing the filename of the crystal structure to evaluate the stability of in the '--path' flag.

### For example:

```bash
./minetep.py --mode p --path input_file -scell u
```

This command will calculate the eigenfrequencies for a file named input_file in the current path in the unit cell only. If there are any imaginary eigenfrequencies, the structure will be deemed unstable.

## Thermodynamic properties
The energy vs. volume and enthalpy vs. pressure curves are computed and plotted for given crystal structures. The bulk modulus of the crystal structures is printed after the run. It requires one or more crystal structures to be passed to the code.

To generate the curves, set the '--mode' flag to 't', and make sure you are passing the filename(s) of the crystal structure(s) in the '--path' flag.

Note: if a structure passed in the arguments hasn't been optimised and/or it is unstable (see Testing - Structure optimisation) the bulk modulus printed out will in general be wrong.

### For example:

The command:

```bash
./minetep.py --mode t --path optfile1.struct optfile2.struct
```
will save a figure called 'Energy vs Volume.png' and another one called 'Enthalpy vs Pressure.png' that compare the two behaviours for the two structures and print their bulk moduli.

# Key files

The key files of the program are the following:
* The main file of the program is [minetep.py](minetep.py), which ties the indiviual modules together into a CLI program
* [structureGenerator.py](structureGenerator.py) contains the structure generation code
* [potential.py](potential.py) contains the Tersoff potential
* [optimizer.py](optimizer.py) contains the structure optimisation code
* [eigenfrequencies.py](eigenfrequencies.py) contains the code for computing the eigenfrequencies of the crystal system
* [E_V.py](E_V.py) contains the code for plotting Energy vs. Volume and Enthalpy vs Pressure
* [plot.py](plot.py) contains the code for plotting any crystal structure contained in files with extension .struct; in the terminal type [./plot.py 'filename.struct'] to show a 3D plot of the unit cell
* All other files are helper files that are used in one of the main files.
* To run one of the test files, move the relevant test files out of the 'test_file' directory and run it as usual.
* The fortran module [hessian.f90](hessian.f90) is a flexible module which can be used to compute the hessian matrix of a general function of n variables.
* The fortran module [tersoff_pot.f90](tersoff_pot.f90) computes the Tersoff potential of a supercell of silicon atoms.


# Testing

## Structure generator
The testing code for the structure generator is provided in [structureGenerator_test.py](testing/structureGenerator_test.py). To run it, simply move it up one directory level.
The code will then randomly generate thousands of crystal structures and print 'success' for every structure that fulfils the conditions imposed on the structures (i.e. if the atoms overlap or not).

## Structure optimisation
A code to test that the forces acting on every atom are zero after the minimisation is included in the file grad.py; to print the vector of the forces move the file to the working directory and type [./grad.py 'opt_filename.struct'] in the terminal where 'opt+filename.struct' is the file written after an optimisation. We successfully used this tool to test the stability of our structures. An additional test of reliability of the optimisation comes from the success in finding the experimentally observed structures; check out the folder test_files to look at some of the structures we found while running the code (included are among the others cubic diamond, simple cubic and hexagonal). The example files were generated using a variety of values for temperature, number of iterations and stepsizes, changed using the rationale explained above; we tested the structures obtained with grad.py to confirm that the gradient was zero and checked that the hessian was positive definite at that minimum through the --mode p (option u) included in this package.

## Hessian
The hessian module was tested using the [hessian_test.py](testing/hessian_test.py) file. To run the file again, move it one directory up and run it as usual.
The output will be the time taken to compute the hessian in the python implementation, vs. the time taken in the fortran implementation.
Following that an example of the hessian of the function $`\exp(-x^2-y^2)`$ is printed as computed by the fortran and python implementation.
The output will be the time taken to compute the hessian in the python implentation, vs. the time taken in the fortran implementation.
Following that an example of the hessian of the function $`exp(-x^2-y^2)`$ (evaluated at the point $`x=y=1`$) is printed as computed by the fortran and python implementations.

## Silicon diamond test file
A file containing the diamond structure in the correct format is provided in [8_cd.struct](test_files/tri_8/8_cd.struct).
# License
This project is available under the ```GNU General Public License v3.0```. See [license.txt](license.txt) for details.
