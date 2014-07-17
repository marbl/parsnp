// See the LICENSE file included with this software for license information.

#include "LCB.hh"
#include <iostream>
#include <fstream>


long Cluster::counter = 1;

Cluster::Cluster()
{
    this->id = counter++;
    type = 1;
}

long Cluster::getId(void)
{
	return this->id;
}

Cluster::Cluster( TMum mum, int type)
: start(mum.start), end(mum.end), length(mum.length)
{
    this->mums.push_back(mum);
    this->id = counter++;
    this->type = type;

}

void Cluster::addMum ( TMum mum )
{
    this->end = mum.end;
    this->length += mum.length;
    this->mums.push_back(mum);
    
}


void Cluster::clear ( void )
{
    this->mums.clear();
    this->end.clear();
    this->length = 0;
    this->start.clear();
    
    
}
void Cluster::print (void)
{
    
    
}
int operator<( const Cluster &c1, const Cluster &c2)
{
    return  c1.start.at(0) < c2.start.at(0);
    
}
Cluster::~Cluster()
{}
