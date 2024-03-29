# Project code for IST Project Advanced Topics in Computational Physics 2022/2023: Project 2


#### Students: 

Catarina Corte-Real ist191035

Francisco Valério Raposo ist196531

Kevin Steiner ist1107611

#### Supervisors:

Prof. David Hilditch

Prof. Alex Vañó-Viñuales

### Explanation

This C++ code solves the wave equation $\Box \phi = 0$ and is parallelized with MPI. It can solve the current state of the field at a given step in time and generate an animation.  There are some example scripts in the section Examples which explain how to do this. Making and running the C++ code as is, will run the simulation and only write some basic information about the simulation parameters to the terminal and not save any data. This is done with

```shell
user:dir$ make parallel_main.exe
user:dir$ mpiexec -np 4 ./bin/parallel_main.exe 2000 0 0
```

Where dir is the base directory of the repository. The command will run on 4 cores and divide the space into 2000 points. Different initial conditions and time discretizations can be changed in the C++ code. The last two zeros are necessary arguments which suppress saving of the solution of the equation and the running of the point convergence test respectively.


### CODE RUNNING INSTRUCTIONS:

Build and run code in main folder with
```shell
user:dir$ make parallel_main.exe
user:dir$ mpiexec -np numberOfCores ./bin/parallel_main.exe NumberOfCells Saving PointConvTest
```

The last two arguments are bools which activate (1) or deactivate (0) the saving of files for an animation or the execution of the point convergence test respectively.

### Animation

```shell
user:dir$ bash animation.sh
```

* T = 2
* N = 3000
* dt = 1E-4

## GITHUB DIRECTORIES ORGANIZATION:

```
.
├── bin				# folder of executables
├── lib				# files necessary for compilation of the code
├── main			# C++ main files for unparallelized code (main.C) and parallelized code (parallel_main.C) 
├── pyt				# python scripts to make animations, and convergence and speedup graphs
│   └── Figures		        # location of figures made by python scripts
└── src				# documentation of C++ class WaveEquationSolver1D, for the unparallelized code

6 directories

```
