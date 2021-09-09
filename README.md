# water-order
This software calculates the tetrahedral order parameter and the d5 parameters for droplets. Mainly aimed for molecular dynamics simulations of water droplets with various molecules/ions in it. Assumes that there is no periodic boundary condition.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (use v1.2.4).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS.

How to compile:
1) First compile chemfiles library. The instructions for this can be found [here](http://chemfiles.org/chemfiles/latest/installation.html).
2) Compile water_order.cpp and link to chemfiles library (static or dynamic—both works)

On Windows (MSVC++)
```
cl /EHSc /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include -IC:\path\to\chemfiles\include Water_order.cpp /link C:\path\to\chemfiles.lib
```
On Linux (GCC)
```
g++ -fexceptions -O3 -fopenmp -I./tclap-1.2.4/include -I/path/to/chemfiles/include/ Water_order.cpp -lchemfiles -L/path/to/chemfiles.lib
```

This should also compile on other operating systems and/or other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the static chemfiles library, also use the /MD flag (you will get linker errors if you try otherwise). I am not sure what happens when linking to the dynamic chemfiles library, but the same runtime should be used to avoid unusual errors during runtime, even if the program compiles and links successfully.

Please cite if you use this software.

Feel free to open an issue if there are any problems.
