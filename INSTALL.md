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

