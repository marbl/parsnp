Parsnp is a command-line-tool for efficient microbial core genome alignment and SNP detection. Parsnp was designed to work in tandem with Gingr, a flexible platform for visualizing genome alignments and phylogenetic trees; both Parsnp and Gingr form part of the Harvest suite :

- [Harvest project page](http://harvest.readthedocs.org)
  -  url: http://harvest.readthedocs.org

Parsnp is primarily distributed as a binary for Linux or OS X (see Parsnp link above ). However, the source and build scripts are provided for incompatible platforms or for development purposes.


# Installation

## From Conda 
```
conda install parsnp --channel bioconda
```
will install Parsnp and all dependencies.

## From source

To build Parsnp from source, users must have Automake 1.15 installed. Parsnp also requires RaxML, Phipack, and numpy. Some additional features require Mash, FastANI and FastTree. All of these packages are available via Conda (many on the Bioconda channel).

### Build instructions
First, you must build the Muscle library
```
cd muscle
./autogen.sh
./configure --prefix=$PWD CXXFLAGS=’-fopenmp’
make install
``

Now we can build Parsnp
```
cd ..
./autogen.sh
./configure
make LDADD=-lMUSCLE-3.7
make install
```

If you wish to be able to move your Parsnp installation around after building, build the parsp binary as follows (after building the Muscle library)
```
./autogen.sh
export ORIGIN=\$ORIGIN
./configure LDFLAGS='-Wl,-rpath,$$ORIGIN/../muscle/lib'
make LDADD=-lMUSCLE-3.7 
make install
```

Note that the `parsnp` executable in `bin/` is not the same as the one in the root level. The former is an alias for Parsnp.py while the latter is the core algorithm of Parsnp that we build above.

# Running Parsnp
Parsnp can be run multiple ways, but the most common is with a set of genomes and a reference. 
```
parsnp -r <reference_genome> -d <genomes 
```
For example, 
```
parsnp -r ref/EMC_2012.gbk -d mers49/*.fna
```

More examples can be found in the [readthedocs tutorial](https://harvest.readthedocs.io/en/latest/content/parsnp/tutorial.html)

CITATION provides details on how to cite Parsnp.

LICENSE provides licensing information.
