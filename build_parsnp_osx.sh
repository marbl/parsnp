export CC="/usr/local/Cellar/gcc49/4.9.2/bin/gcc"
export CXX="/usr/local/Cellar/gcc49/4.9.2/bin/g++"
cd muscle
./autogen.sh
./configure --prefix=`pwd`
make install
cd ..
./autogen.sh
./configure
echo "Fix MUSCLE-3.7 linker" 
make LDADD=-lMUSCLE-3.7
make install
