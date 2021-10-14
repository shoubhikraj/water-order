echo "--- Compiling ChemFiles ---"
cd chemfiles-3fc67b1
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../chemfiles-install -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release --target install
cd ..
cd ..
if [ -e chemfiles-install/lib/libchemfiles.a ]
then
	echo "--- Now compiling Water_order ---"
	if [ -n "${CXX+1}" ]
	then
		echo "--- C++ compiler is $CXX ---"
	else
		echo Setting C++ compiler to default: GNU C++ compiler g++
		CXX=g++
	fi
	$CXX -fexceptions -o Water_order -O3 -fopenmp -I./tclap-1.2.4/include -I./chemfiles-install/include -I./voro-71c84c8/src Water_order.cpp voro-71c84c8/src/voro++.cc -lchemfiles -L./chemfiles-install/lib
else
	echo !!! Error in compiling chemfiles !!!
	echo !!! Please check the CMake logs !!!
fi