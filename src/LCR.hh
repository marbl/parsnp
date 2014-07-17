// See the LICENSE file included with this software for license information.


#ifndef LCR_H
#define LCR_H

#include <vector>
#include <string>

using namespace std;


class TRegion{
  static long counter;
  long id;
public:
  long slength, llength;

  vector<long> start,end, length;
  TRegion();
  TRegion  ( vector<long>& , vector<long>&) ;
  ~TRegion();
  int  getId       ( void          );
  void getLength   ( void          );
  void getStart    ( vector<long>& );
  void getEnd      ( vector<long>& );
  void getRanking  ( vector<long>& );

  friend int operator<(const TRegion &r1, const TRegion &r2);
  friend int operator==(const TRegion &r1, const TRegion &r2);
  	  	
private:
  
  void setStart    ( vector<long>   );
  void setEnd      ( vector<long>   );
  void setLength   ( long          );
  void setRanking  ( vector<long>    );
};

#endif // LCR_H
