cd muscle
./autogen.sh
./configure --prefix=$PWD CXXFLAGS='-fopenmp=libgomp'
make install
cd ..
./autogen.sh
./configure
make LDADD=-lMUSCLE-3.7
make install
