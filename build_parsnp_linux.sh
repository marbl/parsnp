cd muscle
./autogen.sh
./configure --prefix=`pwd`
make install
cd ..
./autogen.sh
export LD_LIBRARY_PATH="/home/travis/build/marbl/parsnp/muscle/lib":$LD_LIBRARY_PATH
./configure --with-libmuscle=/home/travis/build/marbl/parsnp/muscle/lib 
make install
