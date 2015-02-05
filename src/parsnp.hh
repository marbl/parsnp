/////////////////////////////////////////
// parsnp.hh
// main module for parsnp Aligner
// takes a single configuration file as input, output is XMFA
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.

#ifndef ALIGNER_H
#define ALIGNER_H
#define _CRT_SECURE_NO_DEPRECATE 1

// {{{ includes
#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <vector>
#include <map>
#include "LCR.hh"
#include "LCB.hh"
#include "TMum.hh"
// }}}

using namespace std;
typedef string::size_type ssize;
// {{{ class Aligner
class Aligner
{
    // {{{ member variables
    long aligned;
    long m0;
    long filtered;
    long filtered_clusters;
    float l;
    vector<float> gcCount;
    vector<float> atCount;
    bool calc_mumi;
    bool recomb_filter;
    bool filter;
    bool harsh;
    string prefix;
    float factor;
    float diag_diff;
    string outdir;
    //float minmum,minanchor,factor;
    string minmum,minanchor;
    vector<char *> clustalp;
    vector<long> startpos;
    vector<Cluster> clusters;
    vector<TMum> mums;
    vector<TMum> rndmums;
    vector<Cluster> rndclusters;
    vector<TMum> anchors;
    vector<TRegion> regions;
    vector<TRegion> r110;
    // }}}
    
    // {{{ public
public:
    int c, d, q, p, doAlign,cores;
    ssize n;
    bool gridRun;
    bool extendmums;
    bool doUnalign;
    map<string, int> hdr2idx;
    vector< map<int, string> > pos2hdr;
    vector<float>  coverage;
    vector<vector<bool> > mumlayout;
    vector<string> genomes;
    vector<string> files;
    vector<string> fasta;
    vector<string> headers;
    bool shustring;
    float anchorTime,coarsenTime,randomTime,clustersTime,iclustersTime;
    int random;

    Aligner();
    ~Aligner();

     /**
     * parsnp Aligner class constructor.
     * @see testMe()
     * @param genomes string vector containing seq data.
     * @param files genome file paths.
     * @param c min cluster length.
     * @param c min cluster length.
     * @return The test results
     */
    Aligner( vector<string>& ,vector<string>&, int, int , int, int, string, string, bool,  vector<char *>&, vector<string>&
	     , float, bool , vector<float>&, vector<float>& ,bool,int,bool,int,bool,map<string, int>&, vector< map <int,string> > &, vector<string>&, bool calc_mumi, float diag_diff, string prefix, string outdir, bool recomb_filter, bool doUnalign);
    /**
     * parsnp Aligner class destructor
     */

    //test
    void alignClusters( void );
    bool doWork(void );
    /**
     * A pure virtual member.
     * @see testMe()
     * @param c1 the first argument.
     * @param c2 the second argument.
     * @return The test results
     */
    TRegion determineRegion(Cluster,  bool );

    void filterRandom1( int );
    void filterRandomClustersSimple1( void );
    void draw(void);
    void expandCluster(Cluster );

    string  reversec( string seq );
    void  trim( TMum &);
    void setMumi( TRegion, vector<TMum>&, bool, bool );
    void setMums1( TRegion, vector<TMum>&, bool, bool );
    bool setUnalignableRegions( void );
    bool setInitialClusters(void);
    bool setInitialClusters(string);
    void setFinalClusters(void);
    void setFinalClusters(string);
    bool setInterClusterRegions(void);
    static void shuffleSeq(string &seq,int diff);
    static void getSubseq(string seq,int start,int end);
    static char shuffleChar(char b,bool enabled,bool &ft, int i,int diff);
    void formatPpmumInput( vector<long>, vector<long>,  char [] );
    void writeOutput(string psnp, vector<float>& coveragerow);
    
    
    // }}}
};
// }}}

#endif // ALIGNER_H
