echo *** Compiling ChemFiles ***
cd chemfiles-3fc67b1
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../chemfiles-install -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target install
cd ..
cd ..
if [test -e chemfiles-install/lib/libchemfiles.a]
then
	echo *** Compiling Water_order ***
	g++ -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I./chemfiles-install/include Water_order.cpp -lchemfiles -L./chemfiles-install/lib
else
	echo !!! Error in compiling chemfiles !!!
	echo !!! Please check the CMake logs !!!
fi