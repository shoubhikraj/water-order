# Water-order : A fast C++ code to calculate order parameters for water droplets

**(Currently an alpha version, the code is being tested)**

This software calculates the radial distribution of several order parameters, including **(1)** oriental tetrahedral order parameter (*q*<sub>T</sub>) **(2)** fifth nearest neighbour (*d*<sub>5</sub>) parameter, **(3)** Translational tetrahedral order parameter (*S*<sub>k</sub>) and **(4)** Voronoi cell based density (*ρ*<sub>V</sub>) for water droplets. It is mainly aimed at analysis of molecular dynamics simulations of water droplets with various molecules/ions in it. The code assumes that there is no periodic boundary condition.

It calculates the requested order parameter for each oxygen, then takes the average of all oxygens in a certain distance range (`--bin-width`) from the centre-of-mass of the droplet, to get the radial distribution. For a trajectory, the distributions are also averaged over the whole trajectory.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (v1.2.4 already included in this repository).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS (Instructions for compiling are provided at the bottom of the page).

### Is this software suitable for you?
Water-order is mainly aimed at calculations on water droplets. For droplets, the main feature of interest is the radial distribution of the order parameters because it reveals the radial structure of the droplets (as they are mostly spherical). The radius is always calculated from the centre of mass of the droplet.

This is why water-order calculates the radial distribution of various parameters, by averaging over a certain range. It is also assumed that there is no periodic boundary condition (which is true for droplets). Water-order is not for periodic systems, neither does it output data for individual water molecules, although it is possible to make the program output them by modifying the "Water_order.cpp" source file.

## Detailed description
Tetrahedral order parameter(*q*<sub>T</sub>) is a parameter that is used to characterise the degree of order, particularly for water. The *q*<sub>T</sub> for an oxygen can be calculated by considering the 4 nearest oxygens of that oxygen, using this expression:

<p align="center">
   <img src="https://latex.codecogs.com/png.latex?q_T&space;=&space;1&space;-&space;\frac{3}{8}\sum^3_{j=1}\sum^4_{k=j&plus;1}\left(\mathrm{cos}\psi_{jk}&plus;\frac{1}{3}\right)^2" title="q_T = 1 - \frac{3}{8}\sum^3_{j=1}\sum^4_{k=j+1}\left(\mathrm{cos}\psi_{jk}+\frac{1}{3}\right)^2" />
</p>

Where <i>ψ</i><sub>jk</sub> indicates the angle formed by the oxygens j and k and the central oxygen. For ice, <i>q</i><sub>T</sub> is close to 1 (0.8-0.9), whereas for 
liquid water <i>q</i><sub>T</sub> is lower, around 0.5. In the gaseous state, *q*<sub>T</sub> should be close to 0. At intermediate temperatures, there are intermediate values
<br><br>
The fifth nearest neighbour parameter (*d*<sub>5</sub>) is simply the distance to the 5th nearest oxygen from the central oxygen. The value is lower for ice and is higher for liquid water. Both of these parameters are routinely used in calculations involving water, supercooled water and especially on supercooled water droplets.
<br><br>
The translational tetrahedral order parameter (*S*<sub>k</sub>) is similar to the oriental tetrahedral order parameter in that it considers the four nearest neighbours. However, it considers the distances instead of angles.

<p align="center">
   <img src="https://latex.codecogs.com/png.latex?S_k&space;=&space;1-\frac{1}{3}\sum^4_{k=1}\frac{(r_k-\bar{r})^2}{4\bar{r}^2}" title="S_k = 1-\frac{1}{3}\sum^4_{k=1}\frac{(r_k-\bar{r})^2}{4\bar{r}^2}" />
</p>

Where r<sub>k</sub> is the distance between central oxygen and the k-th nearest neighbour and <img src="https://latex.codecogs.com/png.latex?\bar{r}" title="\bar{r}" /> is the arithmetic mean of the distances to the 4 nearest neighbours.
<br><br>
The Voronoi cell based density constructs Voronoi cells for each heavy atom in the system (i.e. excluding hydrogen). The volume of the Voronoi cell associated with each water molecule (because it has one oxygen - which is a heavy atom), is calculated. The density of that molcule is then, the mass of one water molecule divided by the volume of the cell containing it.

This program is calculates the radial distribution of order parameters for droplets/nanodroplets of water. The calculation procedure is as follows:
1) First the centre of mass is evaluated for one frame
2) The order parameter is calculated for each water molecule (or rather, each oxygen) in the frame
3) The program considers spherical shells from the centre of mass (r=0) to the value entered through `--rmax`, by increasing the radius as `--bin-width`. (In other words, it considers radius ranges of size `--bin-width` upto `--rmax`.) Then the average of *q*<sub>T</sub> or *d*<sub>5</sub> values of all oxygens in that range is taken. This gives the radial distribution.
4) Steps 1-3 are iterated for each frame in the trajectory.
5) The radial distributions from all of the frames are summed and then divided by the total number of frames, to get the average radial distribution for the trajectory

The droplet of water can contain ions or other molecules, it only affects the centre of mass, it doesn't affect the rest of the calculation in any way.

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
This message is printed by running `Water_order -h`. For the -f argument, any trajectory file supported by chemfiles can be used (e.g. DCD, XTC, TRR, etc.). For the -s argument (i.e. atom information file), you should enter a PSF file (CHARMM/NAMD/XPLOR), or any valid file format containing mass and name information that can be read by chemfiles. A list of formats read by chemfiles can be found [here](http://chemfiles.org/chemfiles/latest/formats.html). (Right now it does not mention PSF because that functionality has been implemented very recently.)

If the program is run without command line arguments, or some required command line arguments are missing, it will prompt for keyboard input automatically.

Example usage:
```Shell
Water_order -f trajectory.dcd -s topology.psf -c OH2 -t d5 --rmax 20 --bin-width 0.5
```
This calculates the radial distribution of d5 parameter of trajectory.dcd in the distance of 20 Angstroms from the centre of mass of droplet. The d5 values are averaged in ranges of 0.5 Angstroms. The topology file can be a psf file, it can also be any file containing mass information that can be read by chemfiles.

## The output file
The program outputs the radial distribution as a CSV file, containing two columns. The first column is `r` and the second column is the `qT` or `d5` values. The data can be easily plotted with python, or xmgrace, or whatever visualisation software you prefer.

The base name of the output file is determined by the `-o` argument in the command line. When running OTO calculation, `_oto.csv` is appended to the base name. Similarly for d5 calculation, `_d5.csv` is appended to the base name. The default base name is `Water_order`, so the default output files are named `Water_order_oto.csv` and `Water_order_d5.csv`.

An example of plotting the output data with Python (using pandas and matplotlib):
```python
import pandas as pd
import matplotlib.pyplot as plt

# read the CSV file
data = pd.read_csv("Water_order_oto.csv")
# Then plot the x and y axes
plt.plot(data["r"],data["qT"])
# Label the axes
plt.xlabel("r")
plt.ylabel("$q_\mathrm{T}$(r)")
# Show the plot
plt.show()
```

## Current issues
1) If you attempt to run Water_order on a large DCD file (from NAMD, CHARMM or LAMMPS) it will run out of memory. This is due to a problem in the VMD molfile plugin, which is used by chemfiles for reading DCD files. VMD molfile reads the whole trajectory into memory. (This is mentioned in chemfiles issues [here](https://github.com/chemfiles/chemfiles/issues/421) and [here](https://github.com/chemfiles/chemfiles/issues/370).)
2) For histogramming, the sum is accumulated in a `double` type variable. If the trajectory is very large, this can cause overflow. The program checks for the overflow, and will automatically stop. If this happens, please open an issue, and I will try to implement a workaround. 

## How to compile
What is required:
1) C/C++ compiler that supports C++11 standard or higher.
2) CMake

### For Linux or Mac OS X (Or other Unix-like systems)

1) The source of the chemfiles is provided as a zip archive. This is the version which was last tested (You can use newer versions, but there is a chance that it will not be compatible). Chemfiles needs to be compiled for water_order to work. First, uncompress the zip archive (with `unzip` or any other tool), then build with cmake and GNU compilers (or whatever is the default C/C++ compiler in your system).
    ```shell
    unzip chemfiles-592c313.zip
    cmake -S./chemfiles-master -B./chemfiles-master/build -DCMAKE_INSTALL_PREFIX=./chemfiles-install/
    cmake --build . --config Release --target install
    ```
    Further instructions for compiling chemfiles library can be found [here](http://chemfiles.org/chemfiles/latest/installation.html).
2) Check that there is a `libchemfiles.a` file in the `chemfiles-install/lib` directory.
3) Finally, compile Water_order and link to chemfiles library.
```Shell
g++ -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I./chemfiles-install/include/ psfreader_stub.cpp Water_order.cpp -lchemfiles -L./chemfiles-install/lib
```
By default chemfiles is built as a static library, but you can also build a dynamic chemfiles library and link to it.

### For Windows

For Windows, precompiled binaries are provided in the releases section: https://github.com/ShoubhikRaj/water-order/releases.

1) First install Visual Studio 2019 Build Tools (a full install of Visual Studio 2019 can also be done, but only the build tools are required for the purposes of this). Open the Visual Studio Command Line for your platform. This can be found from the search bar in Start Menu. For 64-bit Windows, this is named "x64 Native Tools Command Prompt for VS 2019". (Note that earlier versions of Visual Studio should also work, but has not been tested)
2) Unzip the `chemfiles-592c313.zip` archive. On Windows, the file explorer can do this without any other software. If Winzip or Winrar are used, be careful, because they create an extra level of folder.
3) Build chemfiles with CMake.
    ```batchfile
    cd chemfiles-master
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=../../chemfiles-install -DCMAKE_BUILD_TYPE=Release
    cmake --build . --config Release --target install
4) Check that there is a `chemfiles.lib` in the `chemfiles-install\lib` directory.
5) Finally, compile Water_order and link to chemfiles library.
```batchfile
cl /EHSc /Fe:Water_order /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include /I.\chemfiles-install\include psfreader_stub.cpp Water_order.cpp /link .\chemfiles-install\lib\chemfiles.lib ws2_32.lib advapi32.lib
```
Again, chemfiles can be built as a dynamic library (DLL). In that case you still have to link to the `chemfiles.lib` import library. However, the executable will become dynamically linked to the `chemfiles.dll`, and that DLL has to be in the same folder as the executable, or in PATH.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the chemfiles library, also use the /MD flag (you will get linker errors or runtime errors if you try otherwise).

Water_order should also compile with other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

Note that the program can also be compiled without OpenMP, it will give the same results. Compiling with OpenMP means that the frames are processed parallely in each thread, which improves performance.


----

Please acknowledge if you use this software.

Feel free to open an issue if there are any problems.

## License and Contributing

The software is released under MIT license. Water-order also contains code from tclap-v1.2.4, released under MIT license and chemfiles-v0.11.0, released under BSD-3 clause license. Please read the the [LICENSE](https://github.com/ShoubhikRaj/water-order/blob/main/LICENSE) file for the full terms and conditions of the license.

Any contribution is greatly appreciated. Please open a pull request/issue if you wish to contribute to the code.
