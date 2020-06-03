cd muscle
./autogen.sh
./configure --prefix=$PWD CXXFLAGS='-fopenmp'
make install
cd ..
./autogen.sh
./configure
echo "Fix MUSCLE-3.7 linker" 
make LDADD=-lMUSCLE-3.7
make install
