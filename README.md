# Parsnp 

Rapid bacterial core genome alignment and SNP detection

## Required for building:

* Linux/*nix or OSX (10.7 or newer)
* gcc
* OpenMP
* Python 2.5+ 

## Build

For your convenience, precompiled binaries provides at:

    http://github.com/marbl/harvest

otherwise, to install from source:

    git clone https://github.com/marbl/parsnp.git parsnp_src
    cd parsnp_src
    
First (required!), build libMUSCLE:

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
