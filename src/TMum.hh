// See the LICENSE file included with this software for license information.

#ifndef TMUM_H
#define TMUM_H

#include <vector>
#include <string>
#include "LCR.hh"

using namespace std;
class TMum{
  

  static long counter;
  long id;

 
public:
  vector<long>   start;
  vector<long>   end;
  vector<int>   isforward;
  //vector<string> bps;
  long           length;
  bool hasregion;
  //smallest length of the parent region where the mum was 
  //originally found, used to determine later on if mum is random
  long          slength;
 
  TMum();
  TMum(vector<long>&, long );    
  TMum(vector<long>&, long, vector<int>& , vector<string> &, bool &, bool adjust=true);
  long getid(void);
  void trimleft();
  void trimright();
  ~TMum();

  friend int operator<(const TMum &m1, const TMum &m2);
  friend int operator==(const TMum &m1, const TMum &m2);
  // isRandom, checkRandom, determineRandom
};

#endif // MUM_H
