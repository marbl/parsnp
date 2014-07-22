# Parsnp 

Rapid bacterial core genome alignment and SNP detection

## Required for building:

* 64-bit Linux/*nix or OSX (>= v10.7)
* autoconf && automake
* gcc (>= v4.2.*)
* Python (>= 2.6.*)

## Build

For your convenience, precompiled binaries provides at:

    http://github.com/marbl/harvest

otherwise, to install from source:

    git clone https://github.com/marbl/parsnp.git parsnp_src
    cd parsnp_src
    
Before you start, if running OSX Mavericks, OpenMP is not supported via Clang, so you will not be able to build the source. You will need to install OpenMP and build gcc with OpenMP support. This can be accomplished a couple of ways:

    1.  Install Macports, then:
    
       - sudo port install gcc49
       - sudo port gcc-select mp-gcc49
       
    2.  Install Homebrew, then:
    
       -  brew install gcc49
       
    3.  Build & install gcc from source with OpenMP
    
       - Download & install gcc 4.9
       
          - https://gcc.gnu.org/install/
          
    4.  Download & install gcc prebuilt binaries with OpenMP support
    
       - http://hpc.sourceforge.net/
    
    5. If issues persist, we recommend using the precompiled binary until OpenMP is natively supported by Clang/OSX (likely to be so in Yosemite)
    
Once OpenMP support is added, the first (required!) step is to build libMUSCLE:

    cd muscle
    ./autogen.sh
    ./configure --prefix=`pwd`
    make install

Then, build Parsnp:

    cd ..
    ./autogen.sh
    ./configure
    make install

Once both installed (to cwd install by default):

    export PARSNPDIR=/path/to/parsnp/install

## Quick start

    python Parsnp.py –p <threads> –d <directory of genomes> –r <ref genome>

    Output Files:

    1) Newick formatted core genome SNP tree: <outdir>/parsnp.tree
    2) XMFA format: <outdir>/parsnp.xmfa
    3) Gingr input: <outdir>/parsnp.ggr
    4) VCF variants: <outdir>/parsnp.vcf

## External software dependencies:

* Muscle 3.8 (included as lib)
* PhiPack (included in distribution)
* FastTree2 (included in distribution)

## Docs

See http://harvest.readthedocs.org for full documentation
