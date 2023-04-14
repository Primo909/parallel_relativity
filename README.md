# Project code for IST Project Advanced Topics in Computational Physics 2022/2023: Project 2


#### Students: 

Catarina Corte-Real ist191035

Francisco Valério Raposo ist196531

Kevin Steiner ist1107611

#### Supervisors:

Prof. David Hilditch

Prof. Alex Vañó-Viñuales

### Explanation

This C++ code solves the wave equation $ \box \phi = 0 $ and is parallelized with MPI. It can solve the current state of the field at a given steps in time and generate an animation, it can perform a point convergence and a norm convergence test of the code.  There are some example scripts in the section Examples which explain how to do this. Making and running the C++ code as is, will run the simulation and output the final state of the field into a file. This is done with

```shell
foo@bar:parallel_relativity$ make parallel_main.exe
foo@bar:parallel_relativity$ mpiexec -np 4 ./bin/parallel_main.exe 2000 0 0
```

The command will run on 4 cores and divide the space into 2000 points. Different initial conditions and time discretizations can be changed in the C++ code. The last two zeros are necessary arguments which suppress saving of the functions and the running of the point convergence test respectively.


### CODE RUNNING INSTRUCTIONS:

Build and run code in main folder with
```shell
foo@bar:parallel_relativity$ make parallel_main.exe
foo@bar:parallel_relativity$ mpiexec -np numberOfCores ./bin/parallel_main.exe NumberOfCells Saving PointConvTest
```

The last two arguments are books which activate (1) or deactivate (0) the saving of files for an animation or the execution of the point convergence test respectively.

### Examples
Run example scripts that generate an animation of the solution of the wave equation for the parameters

* T = 1.2
* N = 2000?
* dt = 1E-4?
```shell
user$: bash animation.sh
```

Or that perform a point convergence test for the parameters 
* T = ?
* N = ?
* dt = ?
```shell
user$: bash pointConv.sh
```




## GITHUB DIRECTORIES ORGANIZATION:

```
.
├── build                   # Compiled files (alternatively `dist`)
├── docs                    # Documentation files (alternatively `doc`)
├── src                     # Source files (alternatively `lib` or `app`)
├── test                    # Automated tests (alternatively `spec` or `tests`)
├── tools                   # Tools and utilities
├── LICENSE
└── README.md
```

* Figures/ 

* lib/: files necessary for compilation of the code

* main/: C++ main files for unparallelized code (main.C) and parallelized code (parallel_main.C) - what is file_configuration_test?? AND DELETE OTHER STUFF

* pyt/: python scripts to make animations, and convergence and speedup graphs

* src/: documentation of C++ class WaveEquationSolver1D, for the unparallelized code
