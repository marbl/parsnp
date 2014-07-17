// See the LICENSE file included with this software for license information.

#include "TMum.hh"
#include "csgmum/types.h"
#include <iostream>
long TMum::counter = 1;

TMum::TMum()
{
  this->id = counter++;
  this->hasregion = false;
}
TMum::TMum( vector<long>& start, long length, vector<int>& forward, vector<string> &r, bool &ok, bool adjust)
  :length(length)

{
  
  //need to convert (r-startpos)-length=startpos
  //                endpos = length+startpos
  //                ---->     <----
  //                          ---->
  //	                          ---->

  int j = 0;
  for( vector<int>::iterator bt = forward.begin(); bt != forward.end(); bt++)
  {  
     this->isforward.push_back(*bt);

     if ( *bt or !adjust)
     {
      this->start.push_back(start.at(j));
     }
     else
     {
      this->start.push_back(r.at(j).size() - (start.at(j)+length));
     }
     j++;
  }

  this->hasregion = false;

  slength = 0;
  int size = this->start.size();
  
  for (int i = 0; i < size; i++)
  {
    if ( this->start.at(i) + length > r.at(i).size() )
	{

       ok = 0;

	}
	else if ( this->start.at(i)  < 0 )
	{
       ok = 0;
	}

	else
	{

        if(this->isforward.at(i) || 1)
			this->end.push_back(this->start[i] + length);
		else
            this->end.push_back(this->start[i] - length);
        ok = 1;

	}
   
  }
  this->id = counter++;
  
}

TMum::TMum( vector<long>& start, long length)
  :length(length)

{

  for( vector<long>::iterator it = start.begin(); it != start.end(); it++)
    this->start.push_back(*it);
    
  for( int i =0; i < int(start.size()) ; i++)
    this->isforward.push_back(true);
  
  this->hasregion = false;
  slength = 0;
  int size = this->start.size();
  
  for (int i = 0; i < size; i++)
  {
    this->end.push_back(this->start[i] + length);
  
  }
  this->id = counter++;
}

TMum::~TMum()
{}

long TMum::getid( void )
{
	return this->id;
}
void TMum::trimleft( void )
{
  bool nforward = false;
  for (uint i = 0; i < this->start.size(); i++)
  {
	this->start.at(i) += 1;
	if ( 0 && !this->isforward.at(i) )
	{
        nforward = true;
	}
   
  }
  if ( 0 && nforward )
  {

  for (uint i = 0; i < this->start.size(); i++)
  {
     this->end.at(i) -=1;

  }
  this->length-=1;
  this->slength-=1;
  }
  this->length -=1;
  this->slength -=1;
}

void TMum::trimright( void)
{
  bool nforward = false;
  for (uint i = 0; i < this->start.size(); i++)
  {

    this->end.at(i) -= 1;
	if ( 0 && !this->isforward.at(i) )
	{
        nforward = true;

	}
   
  }

  this->length -=1;
  this->slength -=1;
}


int operator<( const TMum &m1, const TMum &m2)
{
  return  m1.start.at(0) < m2.start.at(0);

}

int operator==( const TMum &m1, const TMum &m2)
{
  for ( int i = 0; i < int(m1.start.size()); i++)
    {

      return (m1.id == m2.id);

    }
  return 1;

}
