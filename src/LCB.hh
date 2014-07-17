/////////////////////////////////////////
// LCR.cpp
// IntraMUM region class
// If collinear across all genomes, will be aligned
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.

#ifndef LCB_H
#define LCB_H

#include "LCR.hh"
#include "TMum.hh"

class Cluster{
    
    
    static long counter;
    long id;
    
public:
    int type;
    vector<TRegion> regions;
    vector<TMum>    mums;
    vector<long>   start;
    vector<long>   end;
    vector<string> alignment;
    long length;
    Cluster();
    Cluster( TMum, int type=1);

    ~Cluster();
    void addMum ( TMum);
    
    void clear ( void );
    long getId ( void );
    void print();
    friend int operator<(const Cluster &c1, const Cluster &c2);

	
private:
	
    
};

#endif // LCB_H
