/////////////////////////////////////////
// LCR.cpp
// IntraMUM region class
// If collinear across all genomes, will be aligned
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.

#include "LCR.hh"

long TRegion::counter = 1;
TRegion::TRegion()
{
  this->id = counter++;
}
TRegion::TRegion( vector<long>& start, vector<long>& end)
{
  long slength = 500000000;
  long llength = 0;
  
  this->id = counter++;
  for( vector<long>::iterator it = start.begin(); it != start.end(); it++)
    this->start.push_back(*it);
  for( vector<long>::iterator it = end.begin(); it != end.end(); it++)
    this->end.push_back(*it);

  for ( int i = 0; i < int(start.size()); i ++)
  {
    this->length.push_back( end.at(i) - start.at(i) );
    if ( this->length.at(i) < slength )
      slength = this->length.at(i);
    if ( this->length.at(i) > llength )
      llength = this->length.at(i);
  }
  this->slength = slength;
  this->llength = llength;
}

TRegion::~TRegion()
{}

int operator<( const TRegion &r1, const TRegion &r2)
{
  return  r1.start.at(0) < r2.start.at(0);

}

int operator==( const TRegion &r1, const TRegion &r2)
{
  for ( int i = 0; i < int(r1.start.size()); i++)
    {
      if ( r1.start.at(i) != r2.start.at(i) || r1.end.at(i) != r2.end.at(i) )
	return 0;

    }
  return 1;

}
int TRegion::getId( void )
{
  return this->id;
}// end TRegion::getId

void TRegion::getLength( void )
{
  //return this->length;
}// end TRegion::getLength

void TRegion::setLength( long l )
{
  //this->length = l;
  
}// end TRegion::setLength

void TRegion::getStart( vector<long>& start )
{
  start = this->start;
}// end TRegion::getStart

void TRegion::setStart( vector<long> start )
{
  for( vector<long>::iterator it = start.begin(); it != start.end(); it++)
    this->start.push_back(*it);

  //this->start = start;
  
}// end TRegion::setStart

void TRegion::getEnd( vector<long>& end )
{
  end = this->end;
}// end TRegion::getEnd

void TRegion::setEnd( vector<long> end )
{
  //this->end = end;
  for( vector<long>::iterator it = end.begin(); it != end.end(); it++)
    this->end.push_back(*it);
}// end TRegion::setEnd

