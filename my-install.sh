export LD_LIBRARY_PATH=/home/Users/vh11/miniconda3/lib
./autogen.sh
./configure
make LDADD=-lMUSCLE-3.7 
make -j install