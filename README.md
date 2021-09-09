# water-order
This software calculates the tetrahedral order parameter and the d5 parameters for droplets. Mainly aimed for molecular dynamics simulations of water droplets with various molecules/ions in it. Assumes that there is no periodic boundary condition.

Depends on the chemfiles library: https://github.com/chemfiles/chemfiles/ and TCLAP library http://tclap.sourceforge.net/ (use v1.2.4).

Compiles successfully on Windows with MSVC++ 2019, on Linux with GCC 10, and it should also compile on Mac OS.

How to compile:
1) First compile chemfiles library.
2) Compile water_order.cpp and link to chemfiles library (static or dynamic—both works)

On Windows (MSVC++)
```
cl /EHSc /O2 /fp:fast /MD water_order.cpp /link chemfiles.lib
```
