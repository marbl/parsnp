cd muscle
./autogen.sh
./configure --prefix=`pwd`
make install
cd ..
./autogen.sh
./configure
make install
