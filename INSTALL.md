ChemFiles currently builds with GNU, Clang and Intel C/C++ compilers on POSIX systems and Visual C++ on Windows systems. So it is recommended that water-order is built with the same compiler that you use for chemfiles. It is possible to use different compilers, but it is not recommended.

### For Linux or Mac OS X (Or other Unix-like systems)

GNU C/C++ compilers are available on most systems by default. Other compilers can usually be loaded by sourcing a script if required. By default chemfiles is built as a static library, but you can also build a dynamic chemfiles library and link to it.

### For Windows

For Windows, precompiled binaries are provided in the releases section: https://github.com/ShoubhikRaj/water-order/releases.

First install Visual Studio 2019 Build Tools (a full install of Visual Studio 2019 can also be done, but only the build tools are required for the purposes of this). Open the Visual Studio Command Line for your platform. This can be found from the search bar in Start Menu. For 64-bit Windows, this is named "x64 Native Tools Command Prompt for VS 2019". (Note that earlier versions of Visual Studio should also work, but has not been tested)

Again, chemfiles can be built as a dynamic library (DLL). In that case you still have to link to the `chemfiles.lib` import library. However, the executable will become dynamically linked to the `chemfiles.dll`, and that DLL has to be in the same folder as the executable, or in PATH.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the chemfiles library, also use the /MD flag (you will get linker errors or runtime errors if you try otherwise).

Water_order should also compile with other compilers, the only thing to make sure is that chemfiles is built with the same compiler, and that C++ exception handling is enabled.

