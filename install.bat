@echo *** Compiling ChemFiles library ***
@rem A script to install Water_order on windows
@cd chemfiles-3fc67b1
@mkdir build
@cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../chemfiles-install -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target install
@cd ..
@cd ..
@setlocal enableextensions
@setlocal enabledelayedexpansion
@if defined CXX (
@echo *** C++ compiler is %CXX% ***
) else (
@echo Setting C++ compiler to default: Visual C++ compiler cl.exe
@set CXX=cl
)
@if exist "chemfiles-install\lib\chemfiles.lib" (
@echo *** Now compiling water_order ***
%CXX% /EHSc /Fe:Water_order /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include /I.\chemfiles-install\include Water_order.cpp /I.\voro-71c84c8\src Water_order.cpp voro-71c84c8\src\voro++.cc /link .\chemfiles-install\lib\chemfiles.lib ws2_32.lib advapi32.lib
) else (
@echo !!! ChemFiles compilation did not finish successfully !!!
@echo !!! Please check the CMake logs !!! 
)
endlocal