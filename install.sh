export LD_LIBRARY_PATH=/home/Users/vh11/miniconda3/lib
./autogen.sh
./configure
make LDADD=-lMUSCLE-3.7 
make install
./parsnp -r examples/mers_virus/ref/England1.fna -d examples/mers_virus/genomes/*.fna --skip-phylogeny