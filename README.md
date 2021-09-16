# Water-order
This software calculates the radial distribution of oriental tetrahedral order parameter (*q*<sub>T</sub>) and the fifth nearest neighbour (*d*<sub>5</sub>) parameters for droplets. Mainly aimed for molecular dynamics simulations of water droplets with various molecules/ions in it. Assumes that there is no periodic boundary condition.

It calculates the order parameter/d5 parameter for each oxygen, then takes the average of all oxygens in a certain distance range (`--bin-width`) from the centre-of-mass of the droplet, to get the radial distribution. For a trajectory, the distributions are also averaged over the whole trajectory.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (v1.2.4 already included in this repository).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS (Instructions for compiling are provided below).

# Detailed description
Tetrahedral order parameter(*q*<sub>T</sub>) is a parameter that is used to characterise the degree of order, particularly for water. The *q*<sub>T</sub> for an oxygen can be calculated by considering the 4 nearest oxygens of that oxygen, using this expression:
<p align="center">
   <img src="https://latex.codecogs.com/png.latex?q&space;=&space;1&space;-&space;\frac{3}{8}\sum^3_{j=1}\sum^4_{k=j&plus;1}\left(\mathrm{cos}\psi_{jk}&plus;\frac{1}{3}\right)^2" title="q = 1 - \frac{3}{8}\sum^3_{j=1}\sum^4_{k=j+1}\left(\mathrm{cos}\psi_{jk}+\frac{1}{3}\right)^2" />
</p>
Where <i>ψ</i><sub>jk</sub> indicates the angle formed by the oxygens j and k and the central oxygen. For ice, *q*<sub>T</sub> is close to 1 (0.8-0.9), whereas for 
liquid water *q*<sub>T</sub> is lower, around 0.5. In the gaseous state, *q*<sub>T</sub> should be close to 0. At intermediate temperatures, there are intermediate values
<br>
The fifth nearest neighbour parameter (*d*<sub>5</sub>) is simply the distance to the 5th nearest oxygen from the central oxygen. The value is lower for ice and is higher for liquid water. Both of these parameters are routinely used in calculations involving water, supercooled water and especially on supercooled water droplets.
<br>
This program is mainly aimed at calculating the radial distribution of *q*<sub>T</sub> and *d*<sub>5</sub> parameters for droplets/nanodroplets of water. The calculation procedure is as follows:
1) First the centre of mass is evaluated for one frame
2) The tetrahedral order parameter, or the fifth nearest neighbour parameter is calculated for each oxygen in the frame
3) The program considers spherical shells from the centre of mass (r=0) to the value entered through `--rmax`, by increasing the radius as `--bin-width`. (In other words, it considers radius ranges of size `--bin-width` upto `--rmax`.) Then the average of *q*<sub>T</sub> or *d*<sub>5</sub> values of all oxygens in that range is taken. This gives the radial distribution.
4) Steps 1-3 are iterated for each frame in the trajectory.
5) The radial distributions from all of the frames are summed and then divided by the total number of frames—to get the average radial distribution for the trajectory

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
```Shell
Water_order -f trajectory.dcd -s topology.psf -c OH2 -t d5 --rmax 20 --bin-width 0.5
```
This calculates the radial distribution of d5 parameter of trajectory.dcd in the distance of 20 Angstroms from the centre of mass of droplet. The d5 values are averaged in ranges of 0.5 Angstroms. The topology file **must** be a psf file.

## Current issues
1) If you attempt to run Water_order on a large DCD file (from NAMD, CHARMM or LAMMPS) it will run out of memory. This is due to a problem in the VMD molfile plugin, which is used by chemfiles for reading DCD files. VMD molfile reads the whole trajectory into memory. (This is mentioned in chemfiles issues [here](https://github.com/chemfiles/chemfiles/issues/421) and [here](https://github.com/chemfiles/chemfiles/issues/370).)

## How to compile
1) First compile chemfiles library. The instructions for this can be found [here](http://chemfiles.org/chemfiles/latest/installation.html).
2) Compile water_order.cpp and link to chemfiles library (static or dynamic—both works)

On Windows (MSVC++)
```Batchfile
cl /EHSc /Fe:Water_order /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include -IC:\path\to\chemfiles\include psfreader_stub.cpp Water_order.cpp /link C:\path\to\chemfiles.lib
```
On Linux (GCC)
```Shell
g++ -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I/path/to/chemfiles/include/ psfreader_stub.cpp Water_order.cpp -lchemfiles -L/path/to/chemfiles.a
```

Water_order should also compile on other operating systems and/or other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the static chemfiles library, also use the /MD flag (you will get linker errors if you try otherwise). I am not sure what happens when linking to the dynamic chemfiles library, but the same runtime should be used to avoid unusual errors during runtime, even if the program compiles and links successfully.

For Windows, precompiled binaries are provided in the releases section: https://github.com/ShoubhikRaj/water-order/releases

----

Please cite if you use this software.

Feel free to open an issue if there are any problems.
