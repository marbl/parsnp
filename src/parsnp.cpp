/////////////////////////////////////////
// parsnp.cpp
// main module for parsnp Aligner
// takes a single configuration file as input, output is XMFA
// author: Todd J Treangen
// email: treangen@gmail.com
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.


#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include "parsnp.hh"
#include "ext/iniFile.h"
#include "Converter.h"
#include "MuscleInterface.h"

#include <cstring>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <map>
#include <sstream>

#define OMP_NUM_THREADS 36

// }}}
// {{{ extern C
extern "C"
{
#include "csgmum/csg.c"
#include "csgmum/mum.c"
}
// }}}

typedef unsigned char  uchar;


// {{{ Aligner::Aligner()
Aligner::Aligner()
{
}
// }}}


// {{{ Aligner::Aligner( vector<string> genomes, vector<string> files, int d, int q,
Aligner::Aligner( vector<string>& genomes , vector<string>& files, int c, int d, int q, int p, string anchors, \
                 string mums, bool filter,  vector<char *>& clustalparams,\
                 vector<string>& fasta, float factor ,bool harsh,vector<float>& gcCount,vector<float>& atCount, bool shustring, int doAlign, bool gridRun, int cores, bool extendmums,map<string, int>& header_to_index, vector < map<int,string> > & pos_to_header, vector<string>& headers, bool calc_mumi, float diag_diff, string prefix, string outdir, bool recomb_filter)
: d(d), q(q),p(p)
{
    
    
    // {{{ variables
    this->diag_diff = diag_diff;
    this->calc_mumi = calc_mumi;
    this->hdr2idx = header_to_index;
    this->recomb_filter = recomb_filter;
    this->doAlign = doAlign;
    this->gridRun = gridRun;
    this->cores = cores;
    this->outdir = outdir;
    this->prefix = prefix;
    this->extendmums = extendmums;
    this->l = 0;
    this->filtered = 0;
    this->filtered_clusters = 0;
    this->m0 = 0;
    this->aligned = 0;
    this->random     = 0;
    this->n          = genomes.size();
    this->shustring = shustring;
    for( vector<string>::iterator it = genomes.begin(); it != genomes.end(); it++)
        this->genomes.push_back(*it);
    genomes.clear();
    for( vector<string>::iterator it = files.begin(); it != files.end(); it++)
        this->files.push_back(*it);
    files.clear();
    for( vector<char *>::iterator it = clustalparams.begin(); it != clustalparams.end(); it++)
        this->clustalp.push_back(*it);
    clustalparams.clear();
    for( vector<string>::iterator it = fasta.begin(); it != fasta.end(); it++)
        this->fasta.push_back(*it);
    fasta.clear();
    for( vector<string>::iterator it = headers.begin(); it != headers.end(); it++)
        this->headers.push_back(*it);
    headers.clear();
    for( vector<float>::iterator it = gcCount.begin(); it != gcCount.end(); it++)
        this->gcCount.push_back(*it);
    gcCount.clear();
    for( vector<float>::iterator it = atCount.begin(); it != atCount.end(); it++)
        this->atCount.push_back(*it);
    atCount.clear();
    for( vector<map<int,string> >::iterator it = pos_to_header.begin(); it != pos_to_header.end(); it++)
        this->pos2hdr.push_back(*it);
    pos_to_header.clear();
    this->filter     = filter;
    this->factor     = factor;
    this->harsh      = harsh;
    this->anchorTime = 0.0;
    this->coarsenTime = 0.0;
    this->randomTime = 0.0;
    this->clustersTime = 0.0;
    this->iclustersTime = 0.0;
    this->c   = c;
    this->d   = d;
    this->q   = q;
    this->p   = p;
    this->minanchor = anchors;
    this->minmum = mums;
    
    // }}}
}
// }}}
// {{{ Constructor Aligner::~Aligner()
Aligner::~Aligner()
{
}

/////////////////////////////////////////
// Aligner::doWork
// Recursive MUM search
// Initial region is entire genome length unless partition length set. These initial MUMs are used as anchors for LCBS
// Recursively searches in smaller regions until the min region size is hit
// Final MUMs are used to anchor alignments within LCBs
/////////////////////////////////////////
bool Aligner::doWork(void)
{
    
    // {{{ variables
    int minsize,ppmumcount;
    vector<TMum> mums;
    vector<TMum>::iterator it;
    vector<TMum>::iterator mt;
    vector<TRegion>::iterator rt;
    vector<Cluster>::iterator ct;
    
    TRegion currRegion, aRegion;
    TRegion lRegion, rRegion;
    
    // consider alternative using regions
    vector<TRegion>::iterator front;
 
    ppmumcount = 0;

    while ( ! this->regions.empty() )
    {
        currRegion = this->regions.front();
        this->regions.erase(this->regions.begin());
        mums.clear();
        
        this->setMums1(currRegion, mums, false, false);
        ppmumcount++;
        int mumcount = int(mums.size());
        front = this->regions.begin();
        
        bool adjacentLeft = false;
        bool adjacentRight = false;
        bool pushedBackLeft = false;
        bool pushedBackRight = false;
        
        // {{{ for (int i = 0; i < mumcount; i++)
        
        for ( int i = 0; i < mumcount; i++)
        {
            
            adjacentLeft = true;
            adjacentRight = true;
            pushedBackLeft = false;
            pushedBackRight = false;
            
            // left
            if ( i == 0 )
            {
                lRegion = this->determineRegion( mums.at(i),  1 );
                adjacentLeft = true;
            }
            
            else if ( lRegion == rRegion  )
            {
                adjacentLeft = true;
            }
            
            // right
            rRegion = this->determineRegion( mums.at(i),  0 );
            
            if ( lRegion.slength > this->q  )
            {
                this->regions.push_back(lRegion);
                pushedBackLeft = true;
            }
            if ( rRegion.slength > this->q )
            {
                this->regions.push_back(rRegion);
                pushedBackRight = true;
            }
            
            if ( i+1 == mumcount )
                adjacentRight = true;
            else
            {
                // left of next mum
                lRegion = this->determineRegion( mums.at(i+1),  1 );
                
                if ( lRegion == rRegion  )
                {
                    adjacentRight = true;
                }
            }
            if ( ! (adjacentLeft) && !( adjacentRight)&& 0 )
            {
                // mum not adjacent to either side
                // if size is deemed to be small, consider random and delete
                minsize = int(0.5*log10(double( mums.at(i).slength))/log10(2.0));
                if ( mums.at(i).length <= minsize )
                {
                    
                    //update mum layout
                    for ( ssize k = 0; k < n; k++)
                    {
                        //punt last change here, <= to <
                        for ( int m = mums.at(i).start[k]; m < mums.at(i).end[k]; m++)
                            this->mumlayout[k][m] = 0;
                    }
                    
                    
                    i--;
                    //update mumcount
                    mumcount -=1;

                    //remove mum region
                    if ( pushedBackLeft )
                        this->regions.pop_back();
                    if ( pushedBackRight)
                        this->regions.pop_back();
                }
            }
            
            else
                this->mums.push_back(mums.at(i));
            
        }
        // }}}
        
        if ( !(this->regions.empty()) )
            sort ( this->regions.begin(), this->regions.end() );
        
        ulong rsize = this->regions.size();
        if ( rsize )
        {
            for ( ulong mregion = 0; mregion < rsize-1; mregion++)
            {
                if( this->regions.at(mregion) == this->regions.at(mregion+1))
                {
                    this->regions.erase(this->regions.begin()+mregion);
                    mregion-=1;
                    rsize -=1;
                }
            }
        }
        
    }
    
    // }}}
    
    bool empty = true;
    if (this->mums.empty() || this->mums.size() == 0 )
        empty = false;
    return ( empty );

}
// }}}

/////////////////////////////////////////
// Aligner::filterRandom1
// Simple MUM filter
// MUMs <= minsize and not collinear to adjacent MUMs are removed
// helps filter out spurious matches 
/////////////////////////////////////////
// {{{ void Aligner::filterRandom(int rvalue)
void Aligner::filterRandom1(int rvalue )
{
    vector<long> start, end;
    vector<TMum> mums;
    TMum mt, nt, prev_mum;
    string sOutput;
    bool adjacentR = true;
    
    sort( this->mums.begin(), this->mums.end());
    
    ulong numums = this->mums.size();
    for ( ulong msize = 0; msize < numums-1; msize++)
    {
        
        //for all clusters, select those that have size < minsize
        // then check to see for those IF they have atleast 1 adjacent mum to them that is <= d
        // if not, delete
        adjacentR = true;
        if(msize>0)
            prev_mum =  this->mums.at(msize-1);
        mt = this->mums.at(msize);

        if ( mt.length <=  rvalue )
        {
            // nt = mt+1;
            nt = this->mums.at(msize+1);
            if ( nt.start.size() < this->n )
            {
                adjacentR = false;
                //break;
            }
            if ( msize == this->mums.size() )
            {
                adjacentR = false;
                //break;
            }
            if( adjacentR )
            {
                for ( int k = 0; k < this->n; k++)
                {
                    
                    if ( abs(nt.start.at(k)) - abs(mt.end.at(k)) < 0  || abs(nt.start.at(k)) - abs(mt.end.at(k)) >  5000 )
                    {
                        adjacentR = false;
                        break;
                    }
                    //for ( int m = mt->end.at(k)+1; m <  mt->end.at(k)+1000; m++)
                    for ( int m = mt.end.at(k)+1; m <  nt.start.at(k); m++)
                    {
                        if(this->mumlayout[k][m])
                        {
                            adjacentR = false;
                            break;
                        }
                    }
                    if ( msize != 0 )
                    {
                        if ( abs(mt.start.at(k)) - abs(prev_mum.end.at(k)) < 0  || abs(mt.start.at(k)) - abs(prev_mum.end.at(k)) >  5000 )
                        {
                            adjacentR = false;
                            break;
                        }
                        //for ( int m = mt->end.at(k)+1; m <  mt->end.at(k)+1000; m++)
                        for ( int m = prev_mum.end.at(k)+1; m <  mt.start.at(k); m++)
                        {
                            if(this->mumlayout[k][m])
                            {
                                adjacentR = false;
                                break;
                            }
                        }
                    }
                    if ( ! adjacentR )
                        break;
                }
            }
            if ( ! adjacentR )
            {
                // mum is not adjacent to any other mum and has small size, delete
                this->filtered+=1;
                this->rndmums.push_back(mt);
                
                for ( ssize k = 0; k < this->n; k++)
                {
                    for ( int m = mt.start.at(k); m <  mt.end.at(k); m++)
                    {
                        this->mumlayout[k][m] = 0;
                    }
                    
                }
                this->mums.erase(this->mums.begin()+msize);

                msize-=1;
                numums-=1;

            }
        }
    }
}
// }}}

/////////////////////////////////////////
// Aligner::filterRandomClustersSimple1
// Simple LCB filter
// Similar to MUM filter, but operates on LCBs
/////////////////////////////////////////
void Aligner::filterRandomClustersSimple1( void )
{
    vector<long> start, end;
    
    Cluster ct;
    string sOutput;
    bool adjacentR = true;
    
    sort( this->clusters.begin(), this->clusters.end());
    TMum mt;
    ulong nsize = 0;
    ulong numclusters = this->clusters.size();
    
    vector<long> clustersToDelete;
    for ( ulong csize = 0; csize < numclusters-1; csize++)
    {
        if ( this->clusters.at(csize).length <= this->c )
        {
            this->filtered_clusters+=1;
            
            this->rndclusters.push_back(this->clusters.at(csize));
            ulong mumid = 0;
            ulong numums =  this->clusters.at(csize).mums.size();
            for ( ulong msize = 0; msize < numums; msize++)
            {
                
                mt = this->clusters.at(csize).mums.at(msize);
                this->filtered+=1;
                
                this->rndmums.push_back(mt);
                
                for ( ssize k = 0; k < this->n; k++)
                {
                    for ( int m = mt.start.at(k); m <  mt.end.at(k); m++)
                    {
                        this->mumlayout[k][m] = 0;
                    }
                    
                }
                mumid = this->clusters.at(csize).mums.at(msize).getid();
                ulong nummums = this->mums.size();
                for ( ulong mmsize = 0; mmsize < nummums; mmsize++)
                {
                    if( this->mums.at(mmsize).getid() == mumid )
                    {
                        this->mums.erase(this->mums.begin()+mmsize);
                        mmsize-=1;
                        nummums-=1;
                        if(mmsize < 0 )
                            mmsize = 0;
                        break;
                    }
                }
            }
            this->clusters.erase(this->clusters.begin()+csize);
            csize-=1;
            if(csize < 0 )
                csize = 0;
            numclusters-=1;
            
        }
        
    }

}
/////////////////////////////////////////
// Aligner::writeOutput
// description: runs libMUSCLE on all LCBs and produces main output
// input: prefix
// output: XMFA file, alignment coverages per genome
/////////////////////////////////////////
// {{{ void Aligner::writeOutput(string psnp, vector<float>& coveragerow)
void Aligner::writeOutput(string psnp,vector<float>& coveragerow)
{
    
    if (this->doAlign)
        cerr << "Step 7: Writing output files & aligning LCBs..." << endl;
    //string psnp = "psnp";
    string prefix = this->outdir;
    prefix.append("/");
    
    
    string test =this->outdir;
    string lcbprefix = this->outdir;
    string lcbdir = this->outdir;
    test.append("/parsnpAligner.log");
    lcbdir.append("/blocks/");
    lcbprefix.append("/blocks/b");
    try
    {
        ofstream testfile ( test.c_str() );
        if (! testfile.good() )
            throw "ParSNP:: output directory does not exist, creating...";
    }
    catch ( char const * str )
    {
        string command = "mkdir ";
        command.append(this->outdir);
        int result = system(command.c_str());
        if ( result )
        {
            cerr << "ParSNP:: error creating output directory, exiting.." << endl;
            exit(1);
        }
    }
    if (this->recomb_filter)
    {
        string command = "mkdir ";
        command.append(lcbdir);
        int result = system(command.c_str());
    }
    string xmfasfile = prefix+ psnp + ".xmfa";
    ofstream xmfafile ( xmfasfile.c_str());
    
    string logfile = prefix+ psnp + ".log";
    ofstream log ( logfile.c_str());
    
    //determine the largest file
    long largest = 0;
    long cur = 0;
    for ( ssize i = 0; i < this->n; i++)
    {
        cur = this->genomes[i].size();
        
        if ( cur > largest )
            largest = cur;
    }
    
    string s1;
    long l1, l2;
    uint cluster_number = 0;
    vector<string> concatalign(this->n,"");
    int chunk = 1;
    int tt = 0;
    int pfilenum = 0;
    int filenum = 0;
    char b[9];
    int nnum = this->n;
    int doalign = this->doAlign;
    bool extendmums = this->extendmums;
    vector<Cluster> allclusters = this->clusters;
    vector<string> allgenomes = this->genomes;
    vector<string> allfasta = this->fasta;
    int total_clusters = 0;
    map <string,int> hdr2idx = this->hdr2idx;
    vector <map <int,string> > pos2hdr = this->pos2hdr;
    for (int zz = 0; zz < allclusters.size(); zz++)
    {
        Cluster ct = allclusters.at(zz);
        if(ct.type==1)
            total_clusters++;
    }
    xmfafile << "#FormatVersion MultiSNiP" << endl;
    xmfafile << "#SequenceCount " << nnum << endl;
    for (int zz = 0; zz < nnum; zz++)
    {
        if ( 1)
        {
            xmfafile << "##SequenceIndex " << zz+1 << endl;
            xmfafile << "##SequenceFile " << allfasta.at(zz) << endl;
            xmfafile << "##SequenceHeader " << this->headers.at(zz) << endl;
            xmfafile << "##SequenceLength " << this->genomes.at(zz).size() << "bp" << endl;
        }
    }
    xmfafile << "#IntervalCount " << total_clusters << endl;
    
    string mumfile2 = "allmums.out";
    ofstream mumfile(mumfile2.c_str());
    std::vector< vector<string> > tempalign2(allclusters.size(),std::vector<string>(nnum));//[allclusters.size()][nnum];
    int numt = this->cores;
    int z = 0;
#pragma omp parallel num_threads(numt)
    {
#pragma omp for schedule(dynamic)
        for ( z = 0; z < allclusters.size(); z++)
        {
            Cluster ct = allclusters.at(z);
            if (ct.type == 1 && ct.mums.size() > 0 && doalign != 0)
            {
                for (ssize i = 0; i < nnum; i++)
                {
                    tempalign2[z][i] = "";
                }
                if (recomb_filter)
                {
                    char bb[9];
                    sprintf(bb,"%d",z+1);// C-style string formed without null
                    
                    string lcbdir = lcbprefix + bb;
                    string lcbfile = lcbdir +"/"+"seq.fna";
                    if (1)
                    {
                        ofstream testfile ( lcbfile.c_str() );
                        if (! testfile.good() )
                        {
                            string command = "mkdir "+lcbdir;
                            int result = system(command.c_str());
                            if ( result )
                            {
                                cout << "ParSNP:: error creating LCB directory: " << lcbdir << endl;
                                exit(1);
                            }
                        }
                        testfile.close();
                    }
                }
            }
            else
                continue;
        }
    }
    
#pragma omp parallel num_threads(numt) shared(tempalign2) private(pfilenum,b)
    {
#pragma omp for schedule(dynamic)//(static, 1)
        for ( tt = 0; tt < allclusters.size(); tt++)
        {
            string s1;
            long l1,l2;
            Cluster ct = allclusters.at(tt);
            if (ct.type == 1 && ct.mums.size() > 0 && doalign != 0 )
            {
                //pass
            }
            else
            {
                continue;
            }
            
            vector<string> tmpalign(nnum,"");
            
            for ( vector<TMum>::iterator dt = ct.mums.begin(); dt != ct.mums.end(); dt++)
            {
                ::std::vector<string> regalign(nnum,"");
                for (ssize i = 0; i < nnum; i++)
                {
                    regalign[i] = "";
                }
                if (ct.mums.size() == 1 && ct.type)
                {
                    vector<TMum>::iterator et = ct.mums.end()-1;
                    //only one MUM, alignment is done!
                    for ( ssize i = 0; i < nnum; i++)
                    {
                        if ( 1)
                        {
                              string s1;
                            if ( !ct.mums.at(0).isforward.at(i) )
                            {
                                
                                string n1 = reversec(allgenomes[i].substr(dt->start.at(i),dt->length));
                                transform(n1.begin(), n1.end(), n1.begin(), ::tolower);
                                tempalign2[tt][i].append(n1);
                            }
                            else
                            {
                                string n1 = allgenomes[i].substr(dt->start.at(i),dt->length);
                                transform(n1.begin(),n1.end(), n1.begin(), ::tolower);
                                tempalign2[tt][i].append(n1);
                            }
                        }
                        
                        //}
                    }
                    //cout << endl;
                }
                //	  }
                //}
                else if ( dt+1 != ct.mums.end() )
                {

                    for ( ssize i = 0; i < nnum; i++)
                    {
                        if (1)
                        {
                            if ( !ct.mums.at(0).isforward.at(i) )
                            {
                                if (i == 0)
                                {
                                    cout << "Error!! MUM in - orientation for ref genome**" << endl;
                                }
                                if ( ct.type)
                                {

                                    string n1 = reversec(allgenomes[i].substr(dt->start.at(i),dt->length));
                                    transform(n1.begin(), n1.end(), n1.begin(), ::tolower);
                                    tempalign2[tt][i].append(n1);
                                    
                                    if (dt->start.at(i)-(dt+1)->end.at(i) >= 1)
                                    {
                                        s1 = reversec(allgenomes[i].substr((dt+1)->end.at(i),(dt->start.at(i)-(dt+1)->end.at(i))));
                                        if (s1.size() >= 1)
                                        {
                                            transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
                                            regalign[i] = s1;
                                        }
                                        else
                                        {
                                            regalign[i] = "";
                                        }
                                    }
                                    else
                                    {
                                        regalign[i] = "";
                                    }
                                    //}
                                }
                            }
                            else
                            {
                                if (ct.type)
                                {
                                    string n1 = allgenomes[i].substr(dt->start.at(i),dt->length);
                                    transform(n1.begin(), n1.end(), n1.begin(), ::tolower);
                                    tempalign2[tt][i].append(n1);
                                    string s1 = allgenomes[i].substr(dt->end.at(i),((dt+1)->start.at(i) - dt->end.at(i)));
                                    if (s1.size() >= 1)
                                    {
                                        transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
                                        regalign[i]  = s1;
                                    }
                                    else
                                    {
                                        regalign[i] = "";
                                    }
                                }
                            }
                        }
                    }
                    
                } //END if dt+1 != ct->mums.end
                
                uint max_length_region = 0;
                uint min_length_region = 1000000;
                vector<bool> skipped(nnum,0);
                uint total_skipped = 0;
                if (ct.type && dt+1 != ct.mums.end() && doalign)
                {
                    vector<string> seq2aln;
                    
                    for (ssize k = 0; k < nnum; k++)
                    {
                        
                        seq2aln.push_back("");
                        if ( 1)
                        {
                            if (regalign[k].size() > max_length_region)
                            {
                                max_length_region = regalign[k].size();
                            }
                            if (regalign[k].size() < min_length_region)
                            {
                                min_length_region = regalign[k].size();
                            }
                        }
                    }
                    if (max_length_region>1 and min_length_region >0)
                    {
                        for ( ssize w = 0; w < nnum; w++)
                        {
                            if ( 1)
                            {
                                int width = 80;
                                ssize t = 0;
                                string s1 = regalign[w];
                                
                                if (s1.size() == 0)
                                {
                                    seq2aln[w] = "N";
                                  }
                                else
                                {
                                    seq2aln[w] = s1;
                                    for( t = 0; t+width < s1.size();)
                                    {
                                        t+= width;
                                    }
                                    
                                }
                            }
                        }
                        if (total_skipped >= nnum-1)
                        {
                            for ( ssize j = 0; j < nnum; j++)
                            {
                                
                                if (1)
                                {
                                    if (skipped.at(j))
                                    {
                                        string ss(max_length_region,'-');
                                        tempalign2[tt][j].append(ss);
                                    }
                                    else
                                    {
                                        tempalign2[tt][j].append(regalign[j]);
                                    }
                                }
                                // }
                            }
                        }
                        else
                        {
                            
                            vector<string> alignment_result;
                            bool ok = true;
                            if (seq2aln.size() > 1)
                            {
                                for (int eze = 0; eze < seq2aln.size(); eze++)
                                {
                                    if (seq2aln[eze].size() < 1 || seq2aln[eze] == "")
                                    {
                                        ok = false;
                                        break;
                                    }
                                }
                                if (ok)
                                {
                                    
                                    bool align_success = false;
                                    MuscleInterface gmi = MuscleInterface();
                                    align_success = gmi.CallMuscleFast(alignment_result,seq2aln);
                                    
                                }
                            }

                            for (int izi =0; izi < alignment_result.size(); izi++)
                            {
                                tempalign2[tt][izi].append(alignment_result[izi]);
                            }
                        }
                    }
                    else if (max_length_region>0 and min_length_region >=0)
                    {
                        for ( ssize p = 0; p < nnum; p++)
                        {
                            if ( 1)
                            {
                                if (regalign[p].size() > 0)
                                {
                                    transform(regalign[p].begin(),regalign[p].end(), regalign[p].begin(), ::toupper);
                                    tempalign2[tt][p].append(regalign[p]);
                                }
                                if (regalign[p].size() < max_length_region)
                                {
                                    int gap_len = max_length_region - regalign[p].size();
                                    for (int ize = 0; ize < gap_len; ize++)
                                    {
                                        tempalign2[tt][p].append("-");
                                    }
                                }
                            }
                            //}
                        }
                    }
                }	// END if (ct->type)
                if (dt+1 == ct.mums.end() && ct.mums.size() > 1 && ct.type)
                {
                    //add last MUM!
                    //vector<TMum>::iterator et = ct->mums.end()-1;
                    //only one MUM, alignment is done!
                    for ( ssize i = 0; i < nnum; i++)
                    {
                        if ( 1)
                        {

                            if ( !ct.mums.at(0).isforward.at(i) )
                            {
                                string n1 = reversec(allgenomes[i].substr(dt->start.at(i),dt->length));
                                transform(n1.begin(),n1.end(), n1.begin(), ::tolower);
                                tempalign2[tt][i].append(n1);
                            }
                            else
                            {
                                string n1 = allgenomes[i].substr(dt->start.at(i),dt->length);
                                transform(n1.begin(),n1.end(), n1.begin(), ::tolower);
                                tempalign2[tt][i].append(n1);
                            }
                        }
                        //}
                    }
                }
            }// END for ( vector<TMum>::iterator dt = ct->mums.begin(); dt != ct->mums.end(); dt++)
            
        }// end for ( vector<Cluster>::iterator ct = this->clusters.begin(); ct != this->clusters.end(); ct++)
    }

    for ( ssize z = 0; z < allclusters.size(); z++)
    {
	    Cluster ct = allclusters.at(z);
	    if (ct.type == 1 && ct.mums.size() > 0 && doalign != 0 && tempalign2[z][0].size() > (this->c * 1))
        {
            sprintf(b,"%d",(int)z+1);// C-style string formed without null
            ofstream clcbfile;
            if (recomb_filter)
            {
                string lcbdir = lcbprefix + b;
                string lcbfile = lcbdir +"/"+"seq.fna";
                clcbfile.open( lcbfile.c_str());
            }
            for ( ssize i = 0; i < nnum; i++)
            {
                if ( 1)
                {
                    int width = 80;
                    
                    ssize k = 0;
                    string s1s = tempalign2[z][i];
                    concatalign.at(i) += s1s;
                    while ( s1s.find('\n') != string::npos )
                    {
                        s1s.erase(  s1s.find('\n'), s1s.find('\n')+1 );
                    }
                    
                    if ( ct.mums.at(0).isforward.at(i) )
                    {
                        xmfafile << ">" << i+1 << ":" << ct.start.at(i) <<  "-" << ct.end.at(i)-1 << " ";
                        if (recomb_filter)
                        {
                            clcbfile << ">" << i+1 << ":" << ct.start.at(i) <<  "-" << ct.end.at(i)-1 << " ";
                        }
                    }
                    else
                    {
                        xmfafile << ">" << i+1 << ":" << ct.mums.back().start.at(i) <<  "-" << ct.mums.front().end.at(i)-1 << " ";
                        if (recomb_filter)
                        {
                            clcbfile << ">" << i+1 << ":" << ct.mums.back().start.at(i) <<  "-" << ct.mums.front().end.at(i)-1 << " ";
                        }
                    }
                    bool hit1 = false;
                    bool hit2 = false;
                    string hdr1 = "";
                    string lasthdr1 = "";
                    int seqstart = 0;
                    int laststart = 0;
                    for(map<int,string>::iterator ctt = pos2hdr[i].begin(); ctt != pos2hdr[i].end(); ctt++)
                    {
                        if (hit1 && (ct.start.at(i) < (*ctt).first))
                        {
                            hit2 = true;
                            hdr1 = lasthdr1;//(*ctt).second;
                            seqstart = laststart;//(*ctt).first;
                            break;
                        }
                        else if (ct.start.at(i) >= (*ctt).first && !hit2)
                        {
                            hit1 = true;
                            laststart = (*ctt).first;
                            lasthdr1 = (*ctt).second;
                            continue;
                        }
                        else if (hit1 & hit2)
                        {
                            hdr1 = lasthdr1;
                            seqstart = laststart;
                            break;
                        }
                        
                    }
                    if (hit1 && !hit2)
                    {
                        hdr1 =  lasthdr1;
                        seqstart = laststart;
                    }
                    if ( !ct.mums.at(0).isforward.at(i) )
                    {
                        xmfafile << "- cluster" << b << " "  << hdr1 << ":p" << (ct.start.at(i)-seqstart)+1 << endl;
                        if(recomb_filter)
                        {
                            clcbfile << "- cluster" << b << " "  << hdr1 << ":p" << (ct.start.at(i)-seqstart)+1 << endl;
                        }
                    }
                    else
                    {
                        xmfafile << "+ cluster" << b << " "  << hdr1 << ":p" << (ct.start.at(i)-seqstart)+1 << endl;
                        if(recomb_filter)
                        {
                            clcbfile << "+ cluster" << b << " "  << hdr1 << ":p" << (ct.start.at(i)-seqstart)+1 << endl;
                        }
                    }
                    for( k = 0; k+width < s1s.size();)
                    {
                        if(recomb_filter)
                        {
                            clcbfile << s1s.substr(k,width) << endl;
                        }
                        xmfafile << s1s.substr(k,width) << endl;
                        k+= width;
                    }
                    if (recomb_filter)
                    {
                        clcbfile << s1s.substr(k,s1s.size()-k) << endl;
                    }
                    xmfafile << s1s.substr(k,s1s.size()-k) << endl;
                }
            }
            xmfafile << "=" << endl;
        }
	    else
	    	continue;
    }
    
    
    float percent;
    //cout << this->aligned <<   " " << this->genomes.at(0).size() << endl;
    
    log << "Number of sequences analyzed:" << setiosflags(ios::fixed) << setprecision(1) << setw(10) << this->n << endl << endl;
    for ( ssize i = 0; i < this->n; i ++)
    {
        if ( 1)
        {
            log <<  "Sequence "<< i+1 << " : "  << this->files.at(i) << endl;
            log <<  this->fasta.at(i) << endl;
            log <<  "Length:"<< setw(10) << this->genomes.at(i).size() << " bps" <<endl;
            log <<  " GC:"<< setw(10) << setiosflags(ios::fixed) << setprecision(1) << this->gcCount.at(i) << endl;
            log <<  " AT:"<< setw(10) << setiosflags(ios::fixed) << setprecision(1) << this->atCount.at(i) << endl;
        }
    }
    
    log <<  setw(2) << setiosflags(ios::left) << "d value:   " << setw(2) << this->d << endl;
    log <<  setw(2) << "q value:   " << setw(2) << this->q << endl << endl;
    log << setw(2) << "Mum anchor size:   " << setw(2) << this->l << endl;
    log << setw(2) <<"Number of MUM anchors found:   "<< setw(2) << this->m0 << endl;
    if((this->mums.size()+this->filtered)>=this->m0)
        log << setw(2) <<  "Number of MUMs found:   "   << setw(2) << (this->mums.size()+this->filtered)-(this->m0) << endl;
    else
        log << setw(2) <<  "Number of MUMs found:   "   << setw(2) << 0 << endl;
    log << setw(2) <<  "Total MUMs found((Anchors+MUMs)-filtered):   "   << setw(2) <<this->mums.size() << endl << endl;
    log << setw(2) << "Random MUM length:   "<< setw(2) << this->random << endl;
    log << setw(2) << "Minimum Cluster length:   "<< setw(2) << this->c << endl;
    log << setw(2) << "Number of MUMs filtered:   "<< setw(2) << this->filtered << endl;
    log << setw(2) << "Number of Clusters filtered:   "<< setw(2) << this->filtered_clusters << endl << endl;
    
    
    percent = (float)(this->aligned)/(float)(this->genomes.at(0).size());
    //count number of actual clusters
    long ccount = 0;
    for(vector<Cluster>::iterator ct = this->clusters.begin(); ct != this->clusters.end(); ct++)
    {
        if (ct->type && ct->mums.size()>0)
            ccount++;
        
        
    }
    log << setw(2) << "Number of clusters created:   "<< setw(2) << ccount << endl;
    if( this->clusters.size() == 0)
        log << setw(2) << "Number of clusters created:   " << setw(2) << "NONE" <<  endl;
    log << setw(2) << "Average number of MUMs per cluster:   "<< setw(2) << this->mums.size()/ccount << endl;
    
    
    long avg = 0;
    vector<long> coverage;
    long totcoverage =0;
    long totsize = 0;
    
    
    for ( ssize i = 0; i < n; i ++ )
    {
        coverage.push_back(0);
    }
    
    
    for ( ssize i = 0; i < n; i ++ )
    {
        if ( 1)
        {
            for ( vector<Cluster>::iterator ct = this->clusters.begin(); ct != this->clusters.end(); ct++)
            {
                if (!ct->type || ct->mums.size()<=0)
                    continue;

                if(ct->mums.at(0).isforward.at(i))
                {
                    coverage.at(i) += abs((ct->mums.back().start.at(i)+ct->mums.back().length)-ct->mums.front().start.at(i));
                    if (i == 0)
                    {
                        avg += abs((ct->mums.back().start.at(i)+ct->mums.back().length)-ct->mums.front().start.at(i));
                    }
                }
                else
                {
                    coverage.at(i) += abs((ct->mums.front().start.at(i)+ct->mums.front().length)-ct->mums.back().start.at(i));
                    if (i == 0)
                    {
                        avg += abs((ct->mums.front().start.at(i)+ct->mums.front().length)-ct->mums.back().start.at(i));
                    }
                }
                
            }
        }
    }
    
    log << setw(2) << "Average cluster length:   "<< avg/ccount <<  " bps" << endl;
    for ( ssize i = 0; i < n; i ++ )
    {
        if (1)
        {
            percent = (float)coverage.at(i)/((float)this->gcCount.at(i)+(float)this->atCount.at(i));
            log << setw(2) << "Cluster coverage in sequence " << i+1 <<  ":   " <<  setiosflags(ios::fixed) << setprecision(1) << 100.00 * percent << "%" <<  endl;
            totcoverage += coverage.at(i);
            totsize += this->genomes.at(i).size();
        }
    }
    percent = (float)totcoverage/(float)totsize;
    log << setw(2) << "Total coverage among all sequences:   " << setiosflags(ios::fixed) << setprecision(1) << 100.00 * percent <<  "%" <<  endl << endl;
    log << setw(2) << " MUM anchor search elapsed time:   " << this->anchorTime << "s " << endl;
    log << setw(2) << " MUM coarsening elapsed time:   " <<  this->coarsenTime << "s " << endl;
    if( this->filter )
        log << setw(2) << " MUM filtering elapsed time:   " << this->randomTime << "s " << endl;
    log << setw(2) << " MUM clustering elapsed time:   " << this->clustersTime << "s " << endl;
    log << setw(2) << " Inter-clustering elapsed time:   " << this->iclustersTime << "s " << endl;
    log << setw(2) << " Total running time:   " << this->anchorTime + this->coarsenTime+this->randomTime+this->clustersTime + this->iclustersTime << "s " << endl;
    this->coverage.push_back(100.00 * percent);
    coveragerow.push_back(100.00 * percent);
    log.close();
}
// }}}
// {{{ TRegion Aligner::determineRegion( Clustet c1, bool LEFT)

/////////////////////////////////////////
// Aligner::determineRegion
// Identify collinear regions to the left & right of initial LCB
/////////////////////////////////////////
TRegion  Aligner::determineRegion( Cluster c1,  bool LEFT )
{

    int mumfound = 0;
    long currpos  = 0;
    vector<long> start;
    vector<long> end;
    vector<string> bps;
    
    // calculate leftmost positions
    if ( LEFT )
    {
        for( ssize i = 0; i < this->n; i++)
        {
            
            if (1)
            {
                mumfound = 0;
                currpos = c1.start[i];
                while ( ! mumfound )
                {
                    currpos -= 1;
                    // no mums found, at beginning of genome
                    if ( currpos < 0 )
                    {
                        mumfound = 1;
                        currpos = 0;
                    }
                    else
                        mumfound = this->mumlayout[i][currpos];
                }
                start.push_back(currpos+1);
                end.push_back(c1.start[i]-1);
            }
            else
            {
                start.push_back(c1.start[i]+1);
                end.push_back(c1.start[i]-1);
            }
        }
        
        TRegion r(start, end );
        start.clear();
        end.clear();
        
        return r;
        
    }
    else
    {
        // calculate rightmost positions
        for( ssize i = 0; i < this->n; i++)
        {
            if ( 1)
            {
                mumfound = 0;
                currpos = c1.end[i];
                while ( ! mumfound )
                {
                    currpos += 1;
                    // no mums found, at end of genome
                    if ( currpos >= long(genomes.at(i).size()) )
                    {
                        mumfound = 1;
                    }
                    else
                        mumfound = this->mumlayout[i][currpos];
                }
                start.push_back(c1.end[i]+1);
                end.push_back(currpos-1);
            }
            else
            {
                start.push_back(c1.end[i]+1);
                end.push_back(genomes.at(i).size()-1);
            }
            if ( currpos > c1.end[i]+1 )
            {
                //pass
            }
        }
        
        TRegion r(start, end );
        start.clear();
        end.clear();
        return r;
    }
    
    //return r;
    // if no mums found between this cluster and any adjacent
    //return mumfound;
}

// }}}

string Aligner::reversec( string seq )
{
    //test
    string genome;
    for( string::iterator start = seq.begin(); start != seq.end(); start++ )
    {
        switch(toupper(*start))
        {
                
            case 'A':
                
                genome.append("T");
                break;
            case 'G':
                
                genome.append("C");
                break;
            case 'C':
                
                genome.append("G");
                break;
            case 'T':
                
                genome.append("A");
                break;
                
            case 'N':
                genome.append("N");
                break;
            case '$':
                genome.append("N");
                break;
                
            case 'n':
                genome.append("N");//
                break;
            case 'x':
                genome.append("N");//
                break;
            case 'X':
                genome.append("N");//
                break;
            case 'Y':
                genome.append("N");//
                break;
            case 'R':
                genome.append("N");//
                break;
            case 'W':
                genome.append("N");//
                break;
            case 'S':
                genome.append("N");//
                break;
            case 'K':
                genome.append("N");//
                break;
            case 'V':
                genome.append("N");//
                break;
            case 'H':
                genome.append("N");//
                break;
            case 'D':
                genome.append("N");//
                break;
            case 'B':
                genome.append("N");//
                break;
            case 'M':
                genome.append("N");//
                break;
            case 'U':
                genome.append("T");//
                break;
            case '\r':
                break;
            case '\n':
                break;
            case '\t':
                break;
            case ' ':
                break;
            case '>':
                break;
            case '.':
                break;
                
            default:
                genome.append("N");
                continue;
        }
               
    }
    // }}}
    
    
    std::reverse(genome.begin(),genome.end());
    return genome;
}

/////////////////////////////////////////
// Aligner::trim
// Trim MUM overlaps
/////////////////////////////////////////
void Aligner::trim( TMum& mum )
{
    
    for (ssize j = 0; j < this->n; j++)
    {
        if ( 1)
        {
            if ( mum.isforward.at(j) || 1 )
            {
                //punt last change <= to <
                for ( long m = mum.start.at(j); m < mum.end.at(j); m++)
                {
                    //remove overlaps, if existing mum overlaps with current mum, trim
                    if (this->mumlayout[j][m] != 0)
                    {
                        mum.trimleft();
                        
                    }
                    else
                    {
                        //cant trim mum in middle, try from right
                        break;
                    }
                }
            }
            else
            {
                // punt last change none to -1
                for ( long m = mum.end.at(j)-1; m >= mum.start.at(j); m--)
                {
                    //remove overlaps, if existing mum overlaps with current mum, trim
                    if (this->mumlayout[j][m] != 0)
                    {
                        mum.trimright();
                    }
                    else
                    {
                        //cant trim mum any further
                        break;
                    }
                }
            }
            if ( mum.isforward.at(j) || 1 )
            {
                for ( long m = mum.end.at(j)-1; m >= mum.start.at(j); m--)
                {
                    //remove overlaps, if existing mum overlaps with current mum, trim
                    if (this->mumlayout[j][m] != 0)
                    {
                        mum.trimright();
                    }
                    else
                    {
                        //cant trim mum any further
                        break;
                    }
                }
            }
            else
            {
                for ( long m = mum.start.at(j); m < mum.end.at(j); m++)
                {
                    //remove overlaps, if existing mum overlaps with current mum, trim
                    if (this->mumlayout[j][m] != 0)
                    {
                        mum.trimleft();                       
                    }
                    else
                    {
                        //cant trim mum any further
                        break;
                        
                    }
                }
            }
        }
    }
    
}
/////////////////////////////////////////
// Aligner::setMums1
// description: For a given region r1, build CSG and identify MUMs
// input: region r1, bool anchors, bool first
// output: a list o MUMs
/////////////////////////////////////////
void Aligner::setMums1(TRegion r1, vector<TMum>& mums, bool anchors = true, bool first = false)
{
    vector<long> mumdata;
    ssize p = 0;
    ssize j = 0;
    float limit = 0;
    int minsize =0;
    vector<string> bps;
    string tbps;
    ulong num_mums = 0;
    SP *SPF,*SPR;
    UM *Master, *MasterRC, *Pair, *PairRC;
    bool extendmums = this->extendmums;
    string sOutput;
    bool firstrun = true;

    // the minsize applies to the genome region with the longest length
    // so that minsize has statistical importance w.r.t llength
    if (anchors)
    {
        Converter(this->minanchor,sOutput,80);
        limit = Calculator(sOutput,sOutput.length(),r1.slength);
        minsize = int(ceil(limit));
        this->l = minsize;
    }
    else
    {
        Converter(this->minmum,sOutput,80);
        limit = Calculator(sOutput,sOutput.length(),r1.slength);
        minsize = int(ceil(limit));
    }

    IRegion* rs = new IRegion[this->n];

    //record current position
    ssize partpos = 0;
    if ( this->p > r1.length.at(0) )
        p = r1.length.at(0);
    else
        p = this->p;
    
    while ( partpos < (ssize)r1.length.at(0) )
    {        
        if ( partpos + p > (ssize)r1.length.at(0) )
        {

            p = r1.length.at(0)-partpos;

            if ( p < 50 )
            {
                p = 50+p;
                partpos = partpos -50;
                //}
            }
        }
        
        rs[0].sequence = (char *)calloc( p+10, sizeof(char));
        rs[0].rc = (char *)calloc( p+10, sizeof(char));
        strcpy(rs[0].sequence,(char *)this->genomes.at(0).substr(r1.start.at(0) + partpos,p).append(1,(char)5).c_str());
        strcpy(rs[0].rc, (char *)reversec(this->genomes.at(0).substr(r1.start.at(0) + partpos,p)).append(1,(char)5).c_str());
        
        rs[0].len_region  = p;
        rs[0].ini_region =  r1.start.at(0)+partpos;
        partpos += p;
        
        if ( firstrun )
        {
            for ( ssize a = 1; a < this->n; a++)
            {
                rs[a].sequence = (char *)calloc( r1.length.at(a)+10, sizeof(char));
                rs[a].rc = (char *)calloc( r1.length.at(a)+10, sizeof(char));
                strcpy(rs[a].sequence,(char *)this->genomes.at(a).substr(r1.start.at(a),r1.length.at(a)).append(1,(char)5).c_str());
                
                strcpy(rs[a].rc,(char *)reversec(this->genomes.at(a).substr(r1.start.at(a),r1.length.at(a))).append(1,(char)5).c_str());
                rs[a].len_region  = r1.length.at(a);
                rs[a].ini_region =  r1.start.at(a);
            }
        }
        
        CSG * csg = 0;
        time_t start,end;
        double dif = 0;
        time ( &start);
        
        if(anchors)
            cerr <<  endl << "        Step 2a: constructing compressed suffix graph...\n";
        csg = new_CSG(csg,int(factor)*rs[0].len_region,rs[0].sequence,rs[0].len_region,0);
        build_CSG(csg, rs[0].sequence, rs[0].len_region, 0);
        find_leaves(csg);
        
        time ( &end);
        dif = difftime(end,start);
        if(anchors)
            printf("                 compressed suffix graph construction elapsed time: %.0lf seconds\n\n",dif);

        Master   = new UM[rs[0].len_region];
        MasterRC   = new UM[rs[0].len_region];
        Pair     =  new UM[rs[0].len_region];
        PairRC   =  new UM[rs[0].len_region];
        SPF = new SP[this->n-1];
        SPR = new SP[this->n-1];
        
        
        time ( &start);
        if(anchors)
            cerr << "        Step 2b: performing initial search for exact matches in the sequences...\n";
        
        for (int i = 0; i < int(rs[0].len_region); i++)
        {
            Pair[i].UP = Pair[i].EP = Master[i].UP = 0;
            PairRC[i].UP = PairRC[i].EP = 0;
            Master[i].EP = rs[0].len_region;
            MasterRC[i].EP = rs[0].len_region;
        }
        ulong *tempMSP;
        tempMSP = new ulong[rs[0].len_region];
        for (ssize i = 0; i < this->n-1; i++)
        {

            SPF[i].MSP = new ulong[rs[0].len_region];
            SPR[i].MSP = new ulong[2];
            SPF[i].forward = new char[rs[0].len_region];
            SPR[i].forward = new char[rs[0].len_region];
            for ( j = 0; j < rs[0].len_region;j++)
            {
                SPF[i].forward[j] = (char)1;
                SPR[i].forward[j] = (char)0;
                tempMSP[j] = 0;
            }

            Find_UM(csg, rs[i+1].sequence, SPF[i].MSP, Pair);
            Find_UM(csg, rs[i+1].rc, tempMSP, PairRC);
            Intersect_UM(csg, Master, Pair, rs[0].len_region, SPF[i].MSP);
            Intersect_UM(csg, MasterRC, PairRC, rs[0].len_region, tempMSP);
            Merge_Master( Master, MasterRC, rs[0].len_region, rs[i+1].len_region, SPF, SPR, tempMSP,i);
        }

        vector<Mum> list_mums;
        vector<Mum> list_shus;
        num_mums = 0;
        int M_EP = 0;
        int M_UP = 0;
        int next_UP = 0,next_EP=0;
        int M_LON = 0;
        int curr_EP = 1, curr_UP = 0;
        //int S_EP = 0;
        int k_SP = 0;
        int smallest_shu =  rs[0].len_region;

        for ( ssize k = 0; k < rs[0].len_region; k++)
        {
            
            
            //SP can be >= when in UM part of MUM
            //SP--UP--------EP
            // SP-UP--------EP
            //  SPUP--------EP
            // ...
            //need to compare with next UP, not previous(M_UP)
            if ( k < rs[0].len_region -1 )
            {
                next_UP = Master[k+1].UP;
                next_EP = Master[k+1].EP;
            }
            else
            {
                next_EP =  k + 1+rs[0].len_region;
                next_UP = k + 1+rs[0].len_region;
            }

            curr_EP = Master[k].EP;
            curr_UP = Master[k].UP;
            k_SP = k;
            if ((Master[k].EP > M_EP) && (Master[k].UP < Master[k].EP))
            {
                
                M_EP = Master[k].EP;
                
                M_LON = M_EP-k;
                if (M_LON >= minsize)
                {
                    Mum  tmum;
                    tmum.forward = (char *)malloc(sizeof(char) * this->n);
                    tmum.DSP =  (ulong *)malloc(sizeof(ulong) * this->n);
                    tmum.forward[0] = 1;
                    
                    if ( tmum.forward[0] )
                        tmum.DSP[0] =  k + 1 + rs[0].ini_region;
                    else
                        tmum.DSP[0] =  (k + 1 + rs[0].ini_region);
    
                    
                    for (ssize j = 1; j < this->n; j++)
                    {
                        tmum.forward[j] = SPF[j-1].forward[k];
                        
                        if ( tmum.forward[j] )
                            tmum.DSP[j] =  SPF[j-1].MSP[k] + 1 + rs[j].ini_region;
                        else
                        {
                            tmum.DSP[j] =  SPF[j-1].MSP[k] + 1 + rs[j].ini_region;
                        }
                    }
                    tmum.LON = M_LON;
                    num_mums++;
                    
                    list_mums.push_back(tmum);
                }
            }
            M_UP = Master[k].UP;
            M_EP = Master[k].EP;
        }
        
       
        delete[] Master;
        delete [] MasterRC;
        delete[] Pair;
        for ( ssize i = 0; i < this->n-1; i++)
        {
            delete SPF[i].MSP;
            delete SPR[i].MSP;
            delete SPF[i].forward;
            delete SPR[i].forward;
        }
        delete [] tempMSP;
        delete [] SPF;
        delete [] SPR;
        

        if (1)
        {
            vector<TMum> unext_mums;
            vector<TMum> unext_mums2;
            for (unsigned int i = 0; i < num_mums; i++)
            {
                vector<long> startpos;
                bool badmum = false;
                for (j = 0; j < this->n; j++)
                {
                    if ( (list_mums[i].DSP[j]-r1.start.at(j) > (unsigned int)r1.length[j] ) )
                    {
                        badmum = true;
                    }
                    else if(0)
                    {
                        //existing mum with this startpos, dont add this mum
                        //since this mum is either not unique or not maximal
                        badmum = true;
                    }
                    else
                    {
                        //check again to see if this mum is not maximal                      
                        int start = 0;
                        bool maximal = false;
                        
                        for ( ssize k = n-1; k < n; k++)
                        {
                            start = this->mumlayout[j][list_mums[i].DSP[j]-1];
                            for ( unsigned int m = list_mums[i].DSP[j]-1; m <= list_mums[i].DSP[j]-1+list_mums[i].LON; m++)
                            {
                                if (this->mumlayout[j][m] != start || this->mumlayout[j][m] == 0 || 1)
                                {
                                    maximal = true;
                                    break;
                                }
                            }
                            if ( maximal )
                            {
                                //for the current sequence this mum is maximal, continue
                                maximal = false;
                            }
                            else
                            {
                                //else this mum is non-maximal, consider it a bad mum and exit loop
                                badmum = true;
                                
                                break;
                            }
                        }
                        
                    }
                    startpos.push_back(list_mums[i].DSP[j]-1);
                }
                if ( badmum )
                    continue;
                else
                {
                    
                    vector<int> isforward;
                    for ( int t = 0 ; t < this->n; t++)
                    {
                        isforward.push_back(list_mums[i].forward[t]);
                    }
                    bool ok = 1;
                    int extend_window = 3*list_mums[i].LON;
                    TMum mum(startpos,list_mums[i].LON,isforward,this->genomes,ok );

                    if ( !ok || mum.length < 5)
                        continue;

                    trim(mum);

                    if ( mum.length < 2 || mum.start.size() <= 1)
                    {
                        continue;
                    }

                    string prevs1 = genomes[0].substr(mum.start.at(0),mum.length);
                    if (!mum.isforward.at(0))
                    {

                        continue;
                    }
                    string s1;
                    long l1, l2;
                    badmum = false;
                    for ( int k = 0; k < this->n; k++)
                    {
                        l1 = mum.start.at(k);
                        if ( !mum.isforward.at(k) && 0  )
                        {
                            l1-=1;
                        }
                        l2 = mum.length;
                        s1 = genomes[k].substr(l1, l2);
                        
                        
                        if ( !mum.isforward.at(k) )
                        {
                            s1 = reversec(s1);
                            if ( prevs1 != s1 )
                            {
                                badmum = true;
                            }
                        }
                        if ( badmum )
                            break;
                        
                        
                    }
                    if ( badmum && 1)
                        continue;

                    for ( ssize k = 0; k < n; k++)
                    {
                        for ( int m = mum.start[k]; m < mum.end[k]; m++)
                            this->mumlayout[k][m] =1;//+= 1;
                    }
                    
                    mum.slength = r1.slength;
                    
                    for( int l = 0; l < 1; l++)
                    {
                        mums.push_back(mum);
                    }

                }
            }
        }

        firstrun = false;
        time ( &end);
        dif = difftime(end,start);
        if(anchors)
            printf("                 MUM anchor search elapsed time: %.0lf seconds\n\n",dif);
        free_CSG(csg);
        
    }
    //} //end omp parallel for on partpos
    
    for ( ssize a = 0; a < n; a++)
    {
        delete rs[a].sequence;
        delete rs[a].rc;
    }
    
    delete[] rs;
    
}
/////////////////////////////////////////
// Aligner::setMumi
// description: For set of genomes, build CSG and calculate MUMi
// input: entire genomic region r1, bool anchors, bool first
// output: MUMi values per genome
/////////////////////////////////////////
void Aligner::setMumi(TRegion r1, vector<TMum>& mums, bool anchors = true, bool first = false)
{
    vector<long> mumdata;
    ssize p = 0;
    ssize j = 0;
    float limit = 0;
    int minsize =0;
    vector<string> bps;
    string tbps;
    ulong num_mums = 0;
    bool extendmums = this->extendmums;
    string sOutput;
    bool firstrun = true;
    
    // the minsize applies to the genome region with the longest length
    // so that minsize has statistical importance w.r.t llength
    if (anchors)
    {
        Converter(this->minanchor,sOutput,80);
        limit = Calculator(sOutput,sOutput.length(),r1.slength);
        minsize = int(ceil(limit));
        this->l = minsize;
    }
    else
    {
        //minsize = int(this->minmum*log10(double(r1.slength))/log10(2.0));
        Converter(this->minmum,sOutput,80);
        limit = Calculator(sOutput,sOutput.length(),r1.slength);
        minsize = int(ceil(limit));
    }
    IRegion* rs = new IRegion[this->n];
    //record current position
    ssize partpos = 0;
    
    //check if genome is larger than partition
    if ( this->p > r1.length.at(0) )
        p = r1.length.at(0);
    else
        p = this->p;

    while ( partpos < (ssize)r1.length.at(0) && partpos < p)
    {       
        if ( partpos + p > (ssize)r1.length.at(0) )
        {
            // if this is true, we have some partition at
            // the end of the sequence which is less than p,
            // assign the following value:
            p = r1.length.at(0)-partpos;
            if ( p < 50 )
            {
                p = 50+p;
                partpos = partpos -50;
            }
        }
        
        rs[0].sequence = (char *)calloc( p+10, sizeof(char));
        rs[0].rc = (char *)calloc( p+10, sizeof(char));
        
        strcpy(rs[0].sequence,(char *)this->genomes.at(0).substr(r1.start.at(0) + partpos,p).append(1,(char)5).c_str());
        strcpy(rs[0].rc, (char *)reversec(this->genomes.at(0).substr(r1.start.at(0) + partpos,p)).append(1,(char)5).c_str());
        rs[0].len_region  = p;
        rs[0].ini_region =  r1.start.at(0)+partpos;
        partpos += p;
        //window overlap in attempt to avoid losing mums at the partition
        if ( firstrun )
        {
            for ( ssize a = 1; a < this->n; a++)
            {
                rs[a].sequence = (char *)calloc( r1.length.at(a)+10, sizeof(char));
                rs[a].rc = (char *)calloc( r1.length.at(a)+10, sizeof(char));
                strcpy(rs[a].sequence,(char *)this->genomes.at(a).substr(r1.start.at(a),r1.length.at(a)).append(1,(char)5).c_str());
                strcpy(rs[a].rc,(char *)reversec(this->genomes.at(a).substr(r1.start.at(a),r1.length.at(a))).append(1,(char)5).c_str());
                rs[a].len_region  = r1.length.at(a);
                rs[a].ini_region =  r1.start.at(a);
            }
        }
              
        time_t start,end;
        double dif = 0;
        time ( &start);
        
        if(anchors)
            cerr <<  endl << "        Step 2a: constructing compressed suffix graph...\n";
        
        
        time ( &end);
        dif = difftime(end,start);
        if(anchors)
            printf("                 compressed suffix graph construction elapsed time: %.0lf seconds\n\n",dif);
        
        time ( &start);
        if(anchors)
            cerr << "        Step 2b: Calculting pairwise MUMi distances...\n";
        
        
        FILE* mumifile;
        string mmf = this->outdir;
        cout << mmf << endl;
        mmf.append("/all.mumi");
        cout << mmf << endl;
        mumifile = fopen(mmf.c_str(),"w");
        int numt = this->cores;
        int iix = 0;

        CSG * csg = 0;// = new CSG;
        csg = new_CSG(csg,int(factor)*rs[0].len_region,rs[0].sequence,rs[0].len_region,0);
        build_CSG(csg, rs[0].sequence, rs[0].len_region, 0);
        find_leaves(csg);
#pragma omp parallel shared(rs,csg) num_threads(numt)
        {
#pragma omp for schedule(dynamic) //static, 1)
            for (iix = 0; iix < this->n-1; iix++)
            {
                map<int, int> amums;
                
                
                
                int total_M_LON = 0;
                UM* Master   = new UM[rs[0].len_region];
                UM* MasterRC   = new UM[rs[0].len_region];
                UM* Pair     =  new UM[rs[0].len_region];
                UM* PairRC   =  new UM[rs[0].len_region];
                SP* SPF = new SP[1];
                SP* SPR = new SP[1];
                SPF[0].MSP = new ulong[rs[0].len_region];
                SPR[0].MSP = new ulong[2];
                SPF[0].forward = new char[rs[0].len_region];
                SPR[0].forward = new char[rs[0].len_region];
                ulong* tempMSP = new ulong[rs[0].len_region];
                
                for ( int j = 0; j < rs[0].len_region;j++)
                {
                    SPF[0].forward[j] = (char)1;
                    SPR[0].forward[j] = (char)0;
                    tempMSP[j] = 0;
                    amums[j] = 0;
                }
                int tid = omp_get_thread_num();
                
                for (int ix = 0; ix < int(rs[0].len_region); ix++)
                {
                    Pair[ix].UP = Pair[ix].EP = Master[ix].UP = 0;
                    PairRC[ix].UP = PairRC[ix].EP = 0;
                    Master[ix].EP = rs[0].len_region;
                    MasterRC[ix].EP = rs[0].len_region;
                }
                
                Find_UM(csg, rs[iix+1].sequence, SPF[0].MSP,Pair);
                Find_UM(csg, rs[iix+1].rc, tempMSP, PairRC);                
                Intersect_UM(csg, Master, Pair, rs[0].len_region, SPF[0].MSP);
                Intersect_UM(csg, MasterRC, PairRC, rs[0].len_region, tempMSP);
                Merge_Master( Master, MasterRC, rs[0].len_region, rs[iix+1].len_region, SPF, SPR, tempMSP, 0);
                int M_EP1 = 0;
                int M_UP1 = 0;
                int next_UP1 = 0,next_EP1=0;
                int M_LON1 = 0;
                int curr_EP1 = 1, curr_UP1 = 0;
                int k_SP1 = 0;
                int totsize = 0;
                for ( ssize k = 0; k < rs[0].len_region; k++)
                {
                    if ( k < rs[0].len_region -1 )
                    {
                        next_UP1 = Master[k+1].UP;
                        next_EP1 = Master[k+1].EP;
                    }
                    else
                    {
                        next_EP1 =  k + 1+rs[0].len_region;
                        next_UP1 = k + 1+rs[0].len_region;
                    }
                    
                    curr_EP1 = Master[k].EP;
                    curr_UP1 = Master[k].UP;
                    k_SP1 = k;

                    if ((Master[k].EP > M_EP1) && (Master[k].UP < Master[k].EP) && (Master[k].EP-k < rs[0].len_region))
                    {
                        totsize +=1;
                        M_EP1 = Master[k].EP;
                        
                        M_LON1 = M_EP1-k;
                        int ii = 0;
                        if (M_LON1 >= 15)
                        {
                            while (ii < M_LON1)
                            {
                                amums[k+ii] = 1;
                                ii++;
                            }
                        }
                        if (M_LON1 >= 15)
                        {
                            total_M_LON += M_LON1;
                        }
                    }
                }
                std::map<int,int>::iterator ip;
                total_M_LON = 0;
                for (ip = amums.begin(); ip != amums.end(); ++ip)
                {
                    total_M_LON += ip->second;
                }
                
                float avg_region_len = float(rs[0].len_region+rs[iix+1].len_region)/2.0;
                int minlen = std::min(rs[0].len_region,rs[iix+1].len_region);
                minlen = rs[0].len_region;
                if (float(r1.length.at(0))/float(rs[iix+1].len_region) > 1.3 || float(r1.length.at(0))/float(rs[iix+1].len_region) < 0.7)
                    total_M_LON = 0;
                if (total_M_LON > minlen)
                    total_M_LON = minlen;
#pragma omp critical
                fprintf(mumifile,"%d:%f\n",iix+1,1.0-(float(total_M_LON)/float(minlen)));
                delete [] Master;
                delete [] MasterRC;
                delete [] Pair;
                delete [] PairRC;
                delete [] tempMSP;
                delete [] SPF[0].forward;
                delete []  SPF[0].MSP;
                delete [] SPR[0].forward;
                delete []  SPR[0].MSP;
                delete SPF;
                delete SPR;
                
            }
        }
        
        //mumifile.close();
        fclose(mumifile);
        free_CSG(csg);
        time ( &end);
        dif = difftime(end,start);
        if(anchors)
            printf("                 MUMi pairwise distance calculation finished: %.0lf seconds\n\n",dif);
            
    }
    
    for ( ssize a = 0; a < n; a++)
    {
        delete rs[a].sequence;
        delete rs[a].rc;
    }
    
    delete[] rs;
    return;
    
}
// {{{ bool Aligner::setInitialClusters(void)
/////////////////////////////////////////
// Aligner::setInitialClusters
// description: Initialize LCBs based on MUMs found
/////////////////////////////////////////
bool Aligner::setInitialClusters(void)
{
    
    vector<long> start, end;
    vector<TMum> mums;
    ssize i = 0;
    while (i < this->n)
    {
        start.push_back(0);
        end.push_back(this->genomes.at(i).length());
        i++;
    }
    
    TRegion sRegion( start, end  );
    if (!this->calc_mumi)
        this->setMums1( sRegion, mums, true, true);
    else
    {
        this->setMumi( sRegion, mums, true, true);
        return -1;
    }
    this->mums = mums;
    this->anchors = mums;
    for ( i = 0; i < mums.size(); i++)
    {
         Cluster cluster(this->mums.at(i));
         this->clusters.push_back(cluster);
    }
    
    TRegion lRegion;
    TRegion rRegion;
    this->m0 = int(mums.size());
    for ( i = 0; i < this->clusters.size(); i ++ )
    {
        lRegion = this->determineRegion(this->clusters.at(i), 1 );
        if ( lRegion.slength > this->q && i == 0 )
        {
            this->regions.push_back(lRegion);
            this->r110.push_back(lRegion);
        }
        else if ( lRegion.slength > this->q && !(lRegion == rRegion))
        {
            this->regions.push_back(lRegion);
            this->r110.push_back(lRegion);
        }
        rRegion = this->determineRegion(this->clusters.at(i),  0 );
        if ( rRegion.slength > this->q && !(rRegion == lRegion) )
        {
            this->regions.push_back(rRegion);
            this->r110.push_back(rRegion);
        }
    }
    return this->m0;
}

// }}}
/////////////////////////////////////////
// Aligner::setInitialClusters
// description: Initialize LCBs based matches contained in anchoFileName
/////////////////////////////////////////
// {{{ bool Alinger::setInitialClusters( string anchorFileName )
bool Aligner::setInitialClusters( string anchorFileName )
{
    
    vector<TMum> mums;
    vector<string> rawmums;
    vector<long> startpos;
    long mumlength;
    char header[80];
    int length;
    char * buffer ;
    int i = 0;
    
    //parse anchorFile
    ifstream is(anchorFileName.c_str());
    is.seekg (0, ios::end);
    length = is.tellg();
    is.seekg (0, ios::beg);
    
    // get length of file:
    // allocate memory:
    buffer = new char [length-2];
    
    is.getline(header,80);
    is.getline(header,80);
    // read data as a block:
    is.read (buffer,length);
    
    string data(buffer);
    int pos = 0;
    int oldpos = 0;
    int mpos = 0;
    int oldmpos = 0;
    long start = 0;
    pos = data.find("\n",0);
    while( pos != int(string::npos) )
    {
        string tempmum( data.substr(oldpos , pos-oldpos ));
        oldpos = pos+1;
        pos = data.find("\n",pos+1);
        rawmums.push_back(tempmum);
        
    }

    for( ssize i = 0; i < rawmums.size(); i++)
    {
        string tempmum;
        if ( i % (this->n + 1) )// i % 3
        {
            continue;
        }
        tempmum = rawmums.at(i);
        mpos = tempmum.find("  ",0);
        oldmpos = 0;
        
        while( mpos != int(string::npos) )
        {
            string mstart = tempmum.substr( oldmpos, mpos-oldmpos );
            start = atol(mstart.c_str());
            if( ! start )
            {
                start = 0;
            }
            startpos.push_back( start );
            oldmpos = mpos+1;
            mpos = tempmum.find("  ",mpos+1);
            
        }
        
        mpos = tempmum.find(" ",oldmpos+1);
        string mlength = tempmum.substr( oldmpos, mpos-oldmpos );
        mumlength = atol(mlength.c_str());
        TMum mum(startpos,mumlength);
        mum.slength = 10000;
        mums.push_back(mum);
        startpos.clear();
        mumlength = 0;

    }
    
    is.close();
    
    for(int i =0; i < int(mums.size());i++)
    {
        for ( ssize k = 0; k < this->n; k++)
        {
            for ( int m = mums.at(i).start[k]; m < mums.at(i).end[k]; m++)
                this->mumlayout[k][m] =1;//+= 1;
        }
        this->anchors.push_back(mums.at(i));
        this->mums.push_back(mums.at(i));
        Cluster cluster(mums.at(i));
        this->clusters.push_back(cluster);
    }
    
    this->m0 = int(mums.size());
    
    TRegion lRegion;
    TRegion rRegion;
    
    for ( i = 0; i < int(this->clusters.size()); i ++ )
    {
        lRegion = this->determineRegion(this->clusters.at(i), 1 );
        if ( lRegion.slength > this->q && i == 0 )
        {
            this->regions.push_back(lRegion);
            this->r110.push_back(lRegion);
        }
        else if ( lRegion.slength > this->q && !(lRegion == rRegion))
        {
            this->regions.push_back(lRegion);
            this->r110.push_back(lRegion);
        }
        rRegion = this->determineRegion(this->clusters.at(i),  0 );
        if ( rRegion.slength > this->q && !(rRegion == lRegion) )
        {
            this->regions.push_back(rRegion);
            this->r110.push_back(rRegion);
        }
    }
    
    return this->m0;
}

// }}}
/////////////////////////////////////////
// Aligner::setUnalignableRegions
// description: Output all non-collinear regions and subsets unable to be aligned
/////////////////////////////////////////
bool Aligner::setUnalignableRegions( void )
{
    string unalignfile = this->outdir+"/parsnp.unalign";
    ofstream unalnfile ( unalignfile.c_str());

    bool stop = 0;
    
    vector<long> lastpos(this->mumlayout.size(),0);
    while ( !stop )
    {
        
        for ( ssize k = 0; k < this->n; k++)
        {
            long startpos = -1;
            long endpos = -1;
            
            
            for( long m = lastpos.at(k); m < this->mumlayout[k].size(); m++)
            {
                if ( this->mumlayout[k][m] == 0 )
                {
                    
                    this->mumlayout[k][m] = 1;
                    if (startpos < 0)
                    {
                        startpos = m;
                        endpos = m;
                    }
                    else
                        endpos +=1;
                    
                }
                else
                {
                    if (startpos >= 0)
                    {
                        lastpos.at(k) = endpos+1;
                        break;
                    }
                    else
                        continue;
                }
            }
            if (startpos != endpos)
            {                
                unalnfile << ">" << k+1  << ":" << startpos << "-" << endpos << " + " << this->fasta.at(k).substr(0,this->fasta.at(k).size()) << endl;
                string  s1 =  this->genomes[k].substr(startpos, endpos-startpos);
                
                long pos = 0;
                while ( pos+80 < s1.size() )
                {
                    
                    unalnfile << s1.substr(pos,80) << endl;
                    pos = pos + 80;
                }
                if ( pos + 1 < s1.size() )
                    unalnfile << s1.substr(pos,s1.size()) << endl;
                unalnfile << "=" << endl;
                
            }
            else if (startpos == -1 && k == this->n-1)
                stop = 1;
            
        }
    }
    
    unalnfile.close();
    return 0;
}
/////////////////////////////////////////
// Aligner::setInterClusterRegions
// need to sort clusters wrt each genome
// only add collinear intercluster regions.. regions that are left out
// due to the d value
/////////////////////////////////////////
bool Aligner::setInterClusterRegions( void )
{
    vector<Cluster> interclusters;
    sort( this->clusters.begin(), this->clusters.end() );
    
    
    for( vector<Cluster>::iterator ct  = this->clusters.begin(); ct+1 < this->clusters.end(); ct++)
    {
        bool add = true;
        vector<long> start;
        vector<long> end;
        long  length = 1;
        //initialize "start mum" with end of current cluster
        int flag = 0;
        for( ssize a = 0; a < this->n; a++)
        {
            
            vector<Cluster>::iterator nt = ct+1;
            
            // if current cluster is not at end, and the start in each sequence is greater than the end of the previous cluster
            if (( nt != this->clusters.end() && nt->start.at(a) - ct->end.at(a) <=0 )&&(1))
            {
                add = false;
                break;
            }
            long stop = this->genomes.at(a).length();
            if (!ct->mums.at(0).isforward.at(a) && 0 )
                start.push_back(ct->end.at(a)-ct->mums.back().length);
            else
                start.push_back(ct->end.at(a));
            for ( long m = ct->end.at(a)+1; m <= stop; m++)
            {
                flag = this->mumlayout[a][m];
                if (flag)
                {
                    end.push_back(m-1);
                    break;
                }
                else
                    continue;
                
            }
            //no mums found, reached end of genome
            if (! flag)
                end.push_back(stop-1);
            
        }
        
        if ( !add )
            continue;
        TMum bmum(start,length);
        TMum emum(end,length);
        Cluster acluster(bmum,0);
        acluster.addMum(emum);
        
        for( ssize a = 0; a < this->n; a++)
        {
            if ((acluster.end.at(a) - acluster.start.at(a)) < 5)
            {
                add = false;
                break;
            }
        }
        if (add)
            interclusters.push_back(acluster);

        
    }

    this->clusters.insert(this->clusters.begin(), interclusters.begin(),interclusters.end());
    return 1;
}
// {{{ void Aligner::setFinalClusters(string mumFileName)
/////////////////////////////////////////
// Aligner::setFinalClusters
// set final LCBs from match file mumFileName
/////////////////////////////////////////
void Aligner::setFinalClusters(string mumFileName)
{
    vector<TMum> mums;
    vector<string> rawmums;
    vector<long> startpos;
    long mumlength;
    char header[80];
    int length;
    char * buffer ;

    ifstream is(mumFileName.c_str());
    is.seekg (0, ios::end);
    length = is.tellg();
    is.seekg (0, ios::beg);

    // allocate memory:
    buffer = new char [length-2];
    
    is.getline(header,80);
    is.getline(header,80);
    
    // read data as a block:
    is.read (buffer,length);
    
    string data(buffer);
    int pos = 0;
    int oldpos = 0;
    int mpos = 0;
    int oldmpos = 0;
    long start = 0;
    pos = data.find("\n",0);
    
    while( pos != int(string::npos) )
    {
        string tempmum( data.substr(oldpos , pos-oldpos ));
        oldpos = pos+1;
        pos = data.find("\n",pos+1);
        rawmums.push_back(tempmum);
        
    }

    for( int i = 0; i < int(rawmums.size()); i++)
    {
        string tempmum;
        if ( i % (this->n + 1) )// i % 3
        {
            continue;
        }
        tempmum = rawmums.at(i);
        mpos = tempmum.find("  ",0);
        oldmpos = 0;
        
        while( mpos != int(string::npos) )
        {            
            string mstart = tempmum.substr( oldmpos, mpos-oldmpos );
            start = atol(mstart.c_str());
            startpos.push_back( start );
            oldmpos = mpos+1;
            mpos = tempmum.find("  ",mpos+1);
            
        }

        mpos = tempmum.find(" ",oldmpos+1);
        string mlength = tempmum.substr( oldmpos, mpos-oldmpos );
        mumlength = atol(mlength.c_str());
        TMum mum(startpos,mumlength);
        mum.slength = 10000;
        mums.push_back(mum);
        startpos.clear();
        mumlength = 0;

    }
    
    is.close();
    for(int i =0; i < int(mums.size());i++)
    {
        for ( ssize k = 0; k < this->n; k++)
        {
            for ( int m = mums.at(i).start[k]; m < mums.at(i).end[k]; m++)
                this->mumlayout[k][m] = 1;
        }
        this->anchors.push_back(mums.at(i));
        this->mums.push_back(mums.at(i));
        Cluster cluster(mums.at(i));
        this->clusters.push_back(cluster);
    }
    
    this->m0 = int(mums.size());
    this->setFinalClusters();
}


/////////////////////////////////////////
// Aligner::setFinalClusters
// set final LCBs from found MUMs
/////////////////////////////////////////
// {{{ void Aligner::setFinalClusters(void)
void Aligner::setFinalClusters(void)
{
    vector<long> start, end;
    vector<TMum> mums;
    vector<TMum>::iterator mt, nt;
    bool addmum  = true;
    long mumcount = 0;
    
    this->clusters.clear();
    // converting old to new clusters
    sort( this->mums.begin(), this->mums.end() );
    mt = this->mums.begin();
    Cluster cluster(*(mt));
    for ( nt  = this->mums.begin()+1; nt < this->mums.end(); nt++)
    {
        
        if (nt->length < this->random)
        {
            addmum = true;
            continue;
        }
        if ( ! addmum  )
        {
            Cluster newcluster( *(nt-1) );
            cluster = newcluster;
        }
        
        mumcount = 0;
        addmum = true;
        float max_length_region = 0;
        float min_length_region= this->d+10;
        float avg_length_region = 0;
        for ( ssize k = 0; k < this->n; k ++ )
        {
            if ( 1)
            {
                if (nt->isforward.at(k))
                    avg_length_region += nt->start.at(k) - cluster.end.at(k);
                else
                    avg_length_region += nt->start.at(k) - cluster.end.at(k);
                
                if ( nt->isforward.at(k) && (nt->start.at(k) - cluster.end.at(k)) >  max_length_region)
                {
                    max_length_region = nt->start.at(k) - cluster.end.at(k);
                }
                else if ( !nt->isforward.at(k) && (cluster.mums.back().start.at(k) - nt->end.at(k))  >  max_length_region)
                {
                    max_length_region = nt->start.at(k) - cluster.end.at(k);
                }
                
                if ( nt->isforward.at(k) && (nt->start.at(k) - cluster.end.at(k)) <  min_length_region)
                {
                    min_length_region = nt->start.at(k) - cluster.end.at(k);
                }
                else if ( !nt->isforward.at(k) && (cluster.mums.back().start.at(k) - nt->end.at(k))  <  min_length_region)
                {
                    min_length_region = cluster.mums.back().start.at(k) - nt->end.at(k);
                }
                if (1 && ( nt->isforward.at(k) != cluster.mums.back().isforward.at(k)) || (nt->isforward.at(k) != cluster.mums.front().isforward.at(k)) )
                {
                    //the overlap is with 2 inverted MUMs, JOIN!
                    addmum = false;
                    
                }
                else if (1 && (nt->isforward.at(k)) && (nt->start.at(k) - cluster.end.at(k) <0))
                {
                    addmum = false;
                }
                else if (1 && (!nt->isforward.at(k)) && ((nt->start.at(k) - cluster.end.at(k)) >= 0 ))
                {
                    addmum = false;
                }
                else if (1 && nt->isforward.at(k) && (nt->start.at(k) - cluster.end.at(k)) >  this->d )
                {
                    addmum = false;
                }
                else if (1 && !nt->isforward.at(k) && (cluster.mums.back().start.at(k) - nt->end.at(k)) >  this->d )
                {
                    addmum = false;
                }
                if ( addmum )
                {
                    // check if there exists a mum between the two mums attempting to join into a cluster
                    if (nt->isforward.at(k))
                    {
                        for ( int m = nt->start.at(k)-1; m > cluster.mums.back().start.at(k)+cluster.mums.back().length; m--)
                        {
                            // if a mum is found, dont add mum to cluster
                            if ( this->mumlayout[k][m] && 0)
                            {
                                addmum = false;
                                break;
                            }
                        }
                    }
                    else
                    {
                        for ( int m = cluster.mums.back().start.at(k)-1; m > nt->end.at(k); m--)
                        {
                            if ( this->mumlayout[k][m] && 0)
                            {
                                addmum = false;
                                break;
                            }
                            
                        }
                    }
                }
                
                if ( ! addmum )
                    break;
            }
        }
        if ( addmum )
        {

            if (min_length_region == 0)
                min_length_region = 1;
            if (max_length_region == 0)
                max_length_region = 1;
            avg_length_region /= float(this->n);
            if (this->diag_diff > 1.0)
            {
                if (max_length_region-min_length_region < this->diag_diff)
                {
                    cluster.addMum( *(nt) );
                    mumcount++;
                }
                
            }
            else if (min_length_region/max_length_region >= 1.0-diag_diff)
            {
                //overlapping inverted mums
                cluster.addMum( *(nt) );
                mumcount++;
            }
            else
            {
                addmum = false;
                this->clusters.push_back(cluster);
            }
        }
        else
        {
            this->clusters.push_back(cluster);
        }
        
    }
    if ( ! addmum  )
    {
        Cluster newcluster( *(nt-1) );
        cluster = newcluster;
        this->clusters.push_back(cluster);
    }
    else
        this->clusters.push_back(cluster);
}
// {{{ void Aligner::shuffleSeq(string &seq, ind diff)
void Aligner::shuffleSeq(string &seq, int diff)
{
    string bases = "AGCT";
    long len = seq.size();
    time_t seconds;
    time(&seconds);
    srand((unsigned int) seconds);
    for(int i = 0; i < len; i ++)
    {
        //1%,5%,10%,20%,50%,100%
        if( rand()%300 <= diff )
        {
            seq.at(i) = bases.at(rand()%3);
        }
    }
}
// }}}
// {{{ void Aligner::getSubseq(string seq, int start, int end)
void Aligner::getSubseq(string seq,int start,int end)
{
    string bases = "AGCT";
    string ssfile = "./output/subseq.fna";
    
    ofstream sfile ( ssfile.c_str());
    sfile << ">subseq " << start <<":" << end << endl;
    for(int i = start; i < end; i++)
    {
        sfile << seq.at(i);
    }
}
// }}}
// {{{ char Aligner::shuffleChar(char b, bool enabled, bool &ft, int i, int diff)
char Aligner::shuffleChar(char b, bool enabled, bool &ft,int i,int diff)
{
    string bases = "AGCT";
    char rvalue = b;
    time_t seconds;
    time(&seconds);
    if(ft)
    {
        srand((unsigned int) seconds);
        ft = false;
    }
    if ( b != 'A' && b != 'G' && b != 'C' && b != 'T')
        return rvalue;
    if( enabled && i )
    {
        if( rand()%300 <= diff )
        {
            rvalue = char(bases.at(rand()%3));
        }
        
    }
    
    return rvalue;
}
// }}}

/////////////////////////////////////////
// Aligner::main
// 1) Instantiate parsnp from configuration file
// 2) Set initial region
// 3) Find anchors
// 4) Determine collinear regions, recursive MUM search
// 5) Set initial LCBs
// 5) Filter spurious matches/LCB
// 6) Set final LCBs
// 7) Align final LCBs
// 8) write output
/////////////////////////////////////////
// {{{ int main(int argc, char* argv[])
int main ( int argc, char* argv[] )
{
    
    // {{{ variables
    //task_scheduler_init init;
    vector<string> genomes,pwgenomes,files,pwfiles,fasta,headers;
    vector<char *> clustalparams;
    vector<bool> mumrow;
    vector< vector<float> > coverages;
    vector< vector<uint> > contig_intervals;
    unsigned int nparams;
    string genome,  ref,query;
    vector<float> gcCount;
    vector<float> atCount;
    char buffer[320], header[2520];
    int i = 0;
    int qfiles = 1 ;
    int c=0, d=0,q=0,p=0,doAlign=0,cores=2;
    bool gridRun = false;
    bool extendmums = false;
    bool calc_mumi = false;
    bool recomb_filter = false;
    float diag_diff = 1.0;
    int random=0;
    float factor;
    string mums,anchors,anchorfile,mumfile,prefix,outdir;
    bool reverseRef=false,reverse=false,anchorsOnly=false;
    bool reverseQuery=false;
    bool shustring = false;
    time_t start,end, tstart,tend;
    double dif;
    long long acount=0,aloc=0,nloc=0,ncount=0,gcount=0,gloc=0,ccount=0,cloc=0,tcount=0,tloc=0;
    map< string, int> header_to_index;
    vector< map < int, string> > pos_to_header;

    if ( argc < 2)
    {
        cout << "ERROR: No ini file specified!" << endl;
        exit(1);
    }
    
    time (&tstart);
    CIniFile iniFile( argv[1] );
    iniFile.ReadFile();
    string cparams = "ClusterParams";
    string cval = "c";
    c = iniFile.GetValueI( cparams, cval);
    d  = iniFile.GetValueI( "ClusterParams", "d");
    diag_diff  = iniFile.GetValueF( "ClusterParams", "diagdiff");
    if (diag_diff < 0.0 || diag_diff > 10000000)
    {
        diag_diff = 1.0;
    }
    q  = iniFile.GetValueI( "ClusterParams", "q");
    p = iniFile.GetValueI(  "ClusterParams","p");
    doAlign = iniFile.GetValueI( "ClusterParams","doalign");
    cores = iniFile.GetValueI( "ClusterParams","cores");
    gridRun = iniFile.GetValueB( "ClusterParams","gridRun");
    recomb_filter = iniFile.GetValueB( "ClusterParams","recombfilter");
    anchors = iniFile.GetValue( "MumParams", "anchors");
    anchorfile = iniFile.GetValue( "MumParams", "anchorfile");
    anchorsOnly = iniFile.GetValueB( "MumParams","anchorsonly");
    calc_mumi = iniFile.GetValueB( "MumParams","calcmumi");
    extendmums = iniFile.GetValueB( "MumParams","extendmums");
    mums = iniFile.GetValue( "MumParams","mums");
    mumfile = iniFile.GetValue( "MumParams","mumfile");
    random = iniFile.GetValueI( "MumParams","filter");
    factor = iniFile.GetValueF( "MumParams","factor");
    //factor = 1.5;
    prefix = iniFile.GetValue( "OutputParams","prefix","parsnp");
    outdir = iniFile.GetValue( "OutputParams","outdir","output");
    
    reverseRef = iniFile.GetValueB( "Reference", "reverse");
    reverseQuery = iniFile.GetValueB( "Query","reverse");
    
    qfiles  = iniFile.NumValues("Query")/2;
    
    bool harsh=false;
    vector<string> allfiles;
    
    string mumatom;
    
    string fname;
    int loc,len;
    time (&start);
    
    pos_to_header.resize(qfiles+1);
    while ( i <= qfiles )
    {
        vector< uint > intervals;
        
        contig_intervals.push_back(intervals);
        reverse = false;
        if ( i == 0 )
        {
            pos_to_header.at(i)[1] = "s1";
            ref = iniFile.GetValue( "Reference", "file");
            reverse = reverseRef;
            ifstream file(ref.c_str());
            loc = ref.rfind("/")+1;
            len = ((ref.size())-loc);
            fname.clear();
            if (loc == int(string::npos))
                fname.append(ref);
            else
                fname.append(ref.c_str(),loc,len);
            allfiles.push_back(fname);
            if ( ! file )
            {
                cout << " Cannot open reference file ! " << endl;
                exit(1);
            }
            else
            {
                file.close();
            }
            files.push_back(ref);
        }
        else
        {
            //Aligner::shuffleSeq(genomes.at(0),4);
            pos_to_header.at(i)[1] = "s1";//header;
            sprintf( buffer,  "file%d",i);
            query = iniFile.GetValue("Query", buffer );
            sprintf( buffer,  "reverse%d",i);
            reverseQuery = iniFile.GetValueB("Query", buffer );
            reverse = reverseQuery;
            ifstream file( query.c_str() );
            loc = query.rfind("/")+1;
            len = (query.size())-loc;
            fname.clear();
            if (loc == int(string::npos))
                fname.append(query);
            else
                fname.append(query.c_str(),loc,len);
            allfiles.push_back(fname);
            if ( ! file )
            {
                cout << " Cannot open query file: " << query <<  endl;
                exit(1);
            }
            else
            {
                file.close();
            }
            files.push_back(query);
        }
        ifstream file( files.at(i).c_str());
        char ch;
        
        file.seekg(0,ios_base::end);
        file.seekg(0, ios_base::beg);
        file.getline( header, 2500 );
        header_to_index[fname] = i;
        headers.push_back(header);
        fasta.push_back( fname);//header );
        unsigned int seqcount = 1;
        acount = 0;
        aloc   = 0;
        gcount = 0;
        gloc   = 0;
        ccount = 0;
        cloc   = 0;
        tcount = 0;
        tloc   = 0;
        ncount = 0;
        nloc = 0;
        
        
        long long counter = 0;
        bool stopimport = false;
        std::stringstream sstm2;
        while ( file.get(ch) and !stopimport)
        {
            counter +=1;
            
            switch(toupper(ch))
            {
                    
                case 'A':
                    acount +=1;
                    aloc   += counter;
                    if ( reverse )
                        genome.append("T");//84
                    else
                        genome.append("A");
                    
                    break;
                    
                case 'G':
                    gcount +=1;
                    gloc   += counter;
                    if ( reverse )
                        genome.append("C");//67
                    else
                        genome.append("G");
                    break;
                    
                case 'C':
                    ccount +=1;
                    cloc   += counter;
                    if ( reverse )
                        genome.append("G");
                    else
                        genome.append("C");//67
                    break;
                    
                case 'T':
                    tcount+=1;
                    tloc   += counter;
                    if ( reverse )
                        genome.append("A");
                    else
                        genome.append("T");//84
                    break;
                case 'X':
                    genome.append("N");
                    break;
                case 'Y':
                    genome.append("N");//
                    break;
                case 'S':
                    genome.append("N");//
                    break;
                case 'W':
                    genome.append("N");//
                    break;
                case 'K':
                    genome.append("N");//
                    break;
                case 'H':
                    genome.append("N");//
                    break;
                case 'U':
                    genome.append("T");
                    break;
                case 'R':
                    genome.append("N");
                    break;
                case 'M':
                    genome.append("N");
                    break;
                case 'V':
                    genome.append("N");
                    break;
                case 'D':
                    genome.append("N");
                    break;
                case 'B':
                    genome.append("N");
                    break;
                case '-':
                    genome.append("N");
                    break;
                case '\n':
                    break;
                case 'n':
                    ncount += 1;
                    cloc += counter;
                    genome.append("N");
                    break;
                case 'N':
                    ncount += 1;
                    cloc += counter;
                    genome.append("N");
                    break;
                case '>':
                    char tmpbuf[2520];
                    file.getline( tmpbuf, 2500 );
                    sstm2.str("");
                    sstm2.clear();
                    if (ncount+ccount+tcount+acount+gcount < 1000)
                        continue;
                    seqcount += 1;
                    sstm2 << "s" << seqcount;
                    pos_to_header.at(i)[ncount+ccount+tcount+acount+gcount] = sstm2.str();//header;
                    contig_intervals[i].push_back(ncount+ccount+tcount+acount+gcount);
                    
                    if (i > 0)
                    {
                        genome.append(d+10,'N');
                        ncount += d+10;
                    }
                    
                    break;
                case '\t':
                case ' ':
                    break;
                default:
                    continue;
            }
            
        }
        // }}}
        
        if(reverse)
            std::reverse(genome.begin(),genome.end());

        contig_intervals[i].push_back(ncount+ccount+tcount+acount+gcount);
        genomes.push_back(genome);
        
        if(1)
        {
            float percent = 0.0;
            percent =  float(gcount)/float(genome.size())*100;
            percent = float(ccount)/float(genome.size())*100;
            percent = float(tcount)/float(genome.size())*100;
            percent = float(acount)/float(genome.size())*100;
  
            gcCount.push_back(float(gcount)+float(ccount));
            atCount.push_back(float(acount)+float(tcount));
            cout << fname << ",Len:" << genome.size() << ",GC:" << ((float(gcount)+float(ccount))/float(genome.size()-ncount))*100 << endl;
        }

        genome.erase();
        file.close();
        i++;
    }
    
    // }}}
    
    
    
    string mfiled =outdir;
    mfiled+="/parsnpAligner.log";
    ofstream mfile ( mfiled.c_str());
    
    
    cerr << "\n*****************************************************" << endl;
    cerr << "\nparsnpAligner:: rapid whole genome SNP typing" << endl;
    cerr << "\n*****************************************************\n" << endl;
    time (&end);
    cerr << "ParSNP: Preparing to construct global multiple alignment framework"<< endl;
    cerr << "\nStep 1: Preparing to verify and process input sequences..." << endl;
    dif = difftime (end,start);
    printf("        Finished processing input sequences, elapsed time: %.0lf seconds\n\n", dif );
    
    Aligner align( genomes, files, c, d, q, p, anchors, mums, random, clustalparams, fasta,factor,harsh,gcCount,atCount,shustring,doAlign,gridRun,cores,extendmums, header_to_index,pos_to_header,headers,calc_mumi,diag_diff,prefix,outdir,recomb_filter);
    for ( ssize i = 0; i < align.n; i ++ )
    {
        align.mumlayout.push_back(mumrow);
        align.mumlayout[i].assign(align.genomes.at(i).size()+1,false);
        align.mumlayout[i][align.genomes.at(i).size()] = true;
    }
    time ( &start);
    if (! calc_mumi)
        cerr << "Step 2: Searching for initial MUM anchors..." << endl;
    else
        cerr << "Calculating mumi distances.." << endl;
    
    bool mumsfound = 0;
    
    if( anchorfile.size() )
    {
        cerr << "*Status: MUM anchor file found " << anchorfile << endl;
        mumsfound = align.setInitialClusters(anchorfile);
        if (mumsfound < 0)
        {
            //returning from calc_mumi
            return 0;
        }
    }
    else if ( ! mumfile.size() )
        mumsfound = align.setInitialClusters();
    
    if (calc_mumi)
        exit(0);
    
    time ( &end );
    
    dif = difftime (end,start);
    align.anchorTime = dif;    
    time ( &start);
    if ( ! anchorsOnly && ! mumfile.size() && ! shustring)
    {
        cerr << "Step 3: Performing recursive MUM search between MUM anchors..." << endl;
        mumsfound = align.doWork();
    }
    time ( &end);
    
    if ( !mumsfound && !mumfile.size())
    {
        
        mfile << "NO MUMS FOUND" << endl;
        mfile.close();
        return 0;
    }
    else
    {
        mfile << "MUMS FOUND" << endl;
        mfile.close();
    }
    
    dif = difftime (end,start);
    printf("        Finished recursive MUM search, elapsed time: %.0lf seconds\n\n", dif );
    align.coarsenTime = dif;
    
    if ( random && ! mumfile.size() )
    {
        cerr << "Step 4: Filtering spurious matches..." << endl;
        time ( &start);
        align.random = random;
        align.filterRandom1(random);
        
        //align.filterRandom();
        time ( &end);
        dif = difftime(end,start);
        printf("        Finished filtering spurious matches, elapsed time: %.0lf seconds\n\n",dif);
        align.randomTime = dif;
    }
    
    
    time ( &start);
    
    cerr << "Step 5: Creating and verifying final LCBs..." << endl;
    if( mumfile.size())
        align.setFinalClusters(mumfile);
    else
        align.setFinalClusters();
    
    align.filterRandomClustersSimple1();
    
    if( mumfile.size())
        align.setFinalClusters(mumfile);
    else
        align.setFinalClusters();
    
    time ( &end);
    dif = difftime(end,start);
    align.clustersTime = dif;
    printf("        Final LCBs verified, elapsed time: %.0lf seconds\n\n",dif);
    time ( &start);
    cerr << "Step 6: Calculating and creating inter-LCB regions..." << endl;
    
    align.setInterClusterRegions();
    time ( &end);
    dif = difftime(end,start);
    align.iclustersTime = dif;
    printf("        LCB regions created, elapsed time: %.0lf seconds\n\n",dif);
    /**/
    vector < float> coverager;

    if (!doAlign)
        cerr << "Step 7: Writing output files..." << endl;
    time ( &start);
    align.writeOutput("parsnpAligner",coverager);
    //align.setUnalignableRegions();
    time ( &end);
    
    dif = difftime (end,start);
    printf("        Output files updated, elapsed time: %.0lf seconds\n\n", dif );
    
    
    time ( &tend);
    dif = difftime (tend,tstart);
    
    cerr << "ParSNP: Finished core genome alignment" << endl;
    printf("        See log file for futher details. Total processing time: %.0lf seconds \n\n", dif );
    exit(0);
}

// }}}
