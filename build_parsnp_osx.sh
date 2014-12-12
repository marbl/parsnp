export CC="/usr/local/Cellar/gcc48/4.8.0/bin/gcc"
export CXX="/usr/local/Cellar/gcc48/4.8.0/bin/g++"
cd muscle
./autogen.sh
./configure --prefix=`pwd`
make install
cd ..
./autogen.sh
./configure
make install
