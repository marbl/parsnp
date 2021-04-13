Parsnp is a command-line-tool for efficient microbial core genome alignment and SNP detection. Parsnp was designed to work in tandem with Gingr, a flexible platform for visualizing genome alignments and phylogenetic trees; both Parsnp and Gingr form part of the Harvest suite :

- [Harvest project page](http://harvest.readthedocs.org)
  -  url: http://harvest.readthedocs.org



# Installation
## From conda
Parsnp is available on the [Bioconda](https://bioconda.github.io/user/install.html#set-up-channels) channel. This is the recommended method of installation and can be installed via
```
conda install parsnp
```

## From source

To build Parsnp from source, users must have automake 1.15, autoconf, and libtool installed. Parsnp also requires RaxML, Phipack, Harvest-tools, and numpy. Some additional features require Mash, FastANI and FastTree. All of these packages are available via Conda (many on the Bioconda channel).

### Build instructions
First, you must build the Muscle library
```
cd muscle
./autogen.sh
./configure --prefix=$PWD CXXFLAGS='-fopenmp'
make install
```

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

## OSX Users (Catalina)
Recent OSX have a Gatekeeper, that's designed to ensure that only softwre from known developers runs on  tour Mac. Please refer to this link to enable the binaries shipped with Parsnp to run: https://support.apple.com/en-us/HT202491

# Running Parsnp
Parsnp can be run multiple ways, but the most common is with a set of genomes and a reference. 
```
parsnp -g <reference_genbank> -d <genomes 
```
```
parsnp -r <reference_fasta> -d <genomes 
```
For example, 
```
./parsnp -g examples/mers_virus/ref/England1.gbk -d examples/mers_virus/genomes/*.fna -c
```
More examples can be found in the [readthedocs tutorial](https://harvest.readthedocs.io/en/latest/content/parsnp/tutorial.html)


## Misc

CITATION provides details on how to cite Parsnp.

LICENSE provides licensing information.
