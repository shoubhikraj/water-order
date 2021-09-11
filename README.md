# Water-order
This software calculates the tetrahedral order parameter and the d5 parameters for droplets. Mainly aimed for molecular dynamics simulations of water droplets with various molecules/ions in it. Assumes that there is no periodic boundary condition.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (v1.2.4 already included in this repository).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS (Instructions for compiling are provided below).

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
This message is printed by running `Water_order -h`. For the -f argument, any trajectory file supported by chemfiles can be used (e.g. DCD, XTC, TRR, etc.). For the -s argument (i.e. atom information file), you should enter a psf file. (Currently only psf files are supported). 

If the program is run without command line arguments, or some required command line arguments are missing, it will prompt for keyboard input automatically.

Example usage:
```
Water_order -f trajectory.dcd -s topology.psf -c OH2 -t d5 --rmax 20 --bin-width 0.5
```
This calculates the radial distribution of d5 parameter of trajectory.dcd in the distance of 20 Angstroms from the centre of mass of droplet. The d5 values are averaged in ranges of 0.5 Angstroms. The topology file **must** be a psf file.

## Current issues
1) If you attempt to run Water_order on a large DCD file (from NAMD, CHARMM or LAMMPS) it will run out of memory. This is due to a problem in the VMD molfile plugin, which is used by chemfiles for reading DCD files. VMD molfile reads the whole trajectory into memory. (This is mentioned in chemfiles issues [here](https://github.com/chemfiles/chemfiles/issues/421) and [here](https://github.com/chemfiles/chemfiles/issues/370).)

## How to compile
1) First compile chemfiles library. The instructions for this can be found [here](http://chemfiles.org/chemfiles/latest/installation.html).
2) Compile water_order.cpp and link to chemfiles library (static or dynamic—both works)

On Windows (MSVC++)
```
cl /EHSc /Fe:Water_order /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include -IC:\path\to\chemfiles\include psfreader_stub.cpp Water_order.cpp /link C:\path\to\chemfiles.lib
```
On Linux (GCC)
```
g++ -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I/path/to/chemfiles/include/ psfreader_stub.cpp Water_order.cpp -lchemfiles -L/path/to/chemfiles.a
```

Water_order should also compile on other operating systems and/or other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the static chemfiles library, also use the /MD flag (you will get linker errors if you try otherwise). I am not sure what happens when linking to the dynamic chemfiles library, but the same runtime should be used to avoid unusual errors during runtime, even if the program compiles and links successfully.

For Windows, precompiled binaries are provided in the releases section: https://github.com/ShoubhikRaj/water-order/releases

----

Please cite if you use this software.

Feel free to open an issue if there are any problems.
