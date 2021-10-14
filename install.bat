echo *** Compiling ChemFiles library ***
@cd chemfiles-3fc67b1
@mkdir build
@cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../chemfiles-install -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target install
@cd ..
@cd ..
if exist "chemfiles-install\lib\chemfiles.lib" (
echo *** Now compiling water_order ***
cl /EHSc /Fe:Water_order /O2 /MD /fp:fast /openmp /I.\tclap-1.2.4\include /I.\chemfiles-install\include Water_order.cpp /link .\chemfiles-install\lib\chemfiles.lib ws2_32.lib advapi32.lib
) else (
echo !!! ChemFiles compilation did not finish successfully !!!
echo !!! Please check the CMake logs !!! 
) 