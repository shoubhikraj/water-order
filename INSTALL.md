ChemFiles currently builds with GNU, Clang and Intel C/C++ compilers on POSIX systems and Visual C++ on Windows systems. So it is recommended that water-order is built with the same compiler that you use for chemfiles. It is possible to use different compilers, but it is not recommended. Also make sure that C++ exception handling is enabled.

If you wish to modify the build process, look at the shell or batch install scripts that have been provided. The build process is not long or complicated, and there are detailed instructions for building chemfiles that is available on their webpages. Compiling water-order is just compiling the source file and linking with the chemfiles library.

### For Linux or Mac OS X (Or other Unix-like systems)

GNU C/C++ compilers are available on most systems by default. Other compilers can usually be loaded by sourcing a script if required. By default chemfiles is built as a static library, but you can also build a dynamic chemfiles library and link to it.

### For Windows

For Windows, precompiled binaries are provided in the releases section: https://github.com/ShoubhikRaj/water-order/releases.

First install Visual Studio 2019 Build Tools (a full install of Visual Studio 2019 can also be done, but only the build tools are required for the purposes of this). Open the Visual Studio Command Line for your platform. This can be found from the search bar in Start Menu. For 64-bit Windows, this is named "x64 Native Tools Command Prompt for VS 2019". (Note that earlier versions of Visual Studio should also work, but has not been tested)

Again, chemfiles can be built as a dynamic library (DLL). In that case you still have to link to the `chemfiles.lib` import library. However, the executable will become dynamically linked to the `chemfiles.dll`, and that DLL has to be in the same folder as the executable, or in PATH.

On Windows, there are multiple versions of C runtimes (static, dynamic, debug etc.). Make sure that chemfiles is build against the same runtime. For the release build of chemfiles, the dynamic multithreaded (/MD) runtime is linked to. When linking to the chemfiles library, also use the /MD flag (you will get linker errors or runtime errors if you try otherwise).

Note that as of now, chemfiles cannot be built with Intel Classic C/C++ or Intel LLVM C/C++ on Windows, so Visual C/C++ is the only compiler that is usable here. You can try the clang-cl provided by Visual Studio and see if it works.
