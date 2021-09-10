# Water-order
This software calculates the tetrahedral order parameter and the d5 parameters for droplets. Mainly aimed for molecular dynamics simulations of water droplets with various molecules/ions in it. Assumes that there is no periodic boundary condition.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (v1.2.4 already included in this repository).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS (Instructions for compiling below).

## How to run
The executable can be run with command line arguments (also without):
```
USAGE:

   Water_order  [-f <string>] [-s <string>] [-c <string>] [-t <OTO|d5>]
                    [-o <string>] [--rmax <float (Angstrom)>] [--bin-width
                    <float (Angstrom)>] [--start <positive integer>]
                    [--stop <positive integer>] [--] [--version] [-h]


Where:

   -f <string>,  --trajectory <string>
     Name of the trajectory file

   -s <string>,  --atom-info <string>
     Name of the file containing the atom names

   -c <string>,  --oxygen-name <string>
     Name of oxgyen atoms in the atom-info file

   -t <OTO|d5>,  --task <OTO|d5>
     Task requested to run: OTO = Oriental Tetrahedral Order, d5 = d5
     parameter; default is OTO

   -o <string>,  --output-file <string>
     Base name of output file; default is 'Water_order'

   --rmax <float (Angstrom)>
     Maximum distance from centre-of-mass considered for distribution (rmin
     = 0 always)

   --bin-width <float (Angstrom)>
     Size of each range considered for averaging

   --start <positive integer>
     Frame number to start calculation from; default is first frame
     (counting starts from 0)

   --stop <positive integer>
     Frame number to end calculation at; default is last frame (counting
     starts from 0)

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Water order analysis
```
This message is printed by running `Water_order -h`. For the -f argument, any trajectory file supported by chemfiles can be used (e.g. DCD, XTC, TRR, etc.). For the -s argument (i.e. atom information file), you should enter a pdb file. (Currently only VMD generated pdb files are supported). 

If the program is run without command line arguments, or some required command line arguments are missing, it will prompt for keyboard input automatically.

## Current issues
1) If you attempt to run Water_order on a large DCD file (from NAMD, CHARMM or LAMMPS) it will run out of memory. This is due to a problem in the VMD molfile plugin, which is used by chemfiles for reading DCD files. VMD molfile reads the whole trajectory into memory. (This is mentioned in chemfiles issues [here](https://github.com/chemfiles/chemfiles/issues/421) and [here](https://github.com/chemfiles/chemfiles/issues/370).)
2) Currently chemfiles does not support psf files. So, to obtain the masses of atoms, masses are inferred from atom names in the pdb file. This means only VMD generated pdb files are supported right now. This should be fixed within a few months. (Atom names OH2 = Water O; H1, H2 = Water H; CAL = Ca2+; FLU = F-)

## How to compile
1) First compile chemfiles library. The instructions for this can be found [here](http://chemfiles.org/chemfiles/latest/installation.html).
2) Compile water_order.cpp and link to chemfiles library (static or dynamic—both works)

On Windows (MSVC++)
```
cl /EHSc /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include -IC:\path\to\chemfiles\include Water_order.cpp /link C:\path\to\chemfiles.lib
```
On Linux (GCC)
```
g++ -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I/path/to/chemfiles/include/ Water_order.cpp -lchemfiles -L/path/to/chemfiles.a
```

This should also compile on other operating systems and/or other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the static chemfiles library, also use the /MD flag (you will get linker errors if you try otherwise). I am not sure what happens when linking to the dynamic chemfiles library, but the same runtime should be used to avoid unusual errors during runtime, even if the program compiles and links successfully.

----

Please cite if you use this software.

Feel free to open an issue if there are any problems.
