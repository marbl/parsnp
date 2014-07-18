aclocal
automake --add-missing
autoconf
cd muscle
./autogen.sh
./configure --prefix=`pwd`
