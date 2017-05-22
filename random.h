//------------------------------------------------------------------------------
//	
//	random.h
//	
//------------------------------------------------------------------------------
#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <iostream>

#define _GA_RND             ran2
#define _GA_RND_SEED        sran2

#define STA_CAST(type,x) ((type)(x))
#define DYN_CAST(type,x) ((type)(x))

void sran2(unsigned int seed=1);
float ran2();

inline int RandomInt(){ return _GA_RND() > 0.5 ? 1 : 0; }
inline int RandomInt(int low, int high){ 
  float val=STA_CAST(float,high)-STA_CAST(float,low)+(float)1; 
  val*=_GA_RND(); 
  return (STA_CAST(int,val)+low);
}

inline double RandomDouble(){ return _GA_RND(); }
inline double RandomDouble(double low, double high){
  double val=high-low; val*=_GA_RND(); return val+low;
}

inline float RandomFloat(){ return _GA_RND(); }
inline float RandomFloat(float low, float high){
  float val=high-low; val*=_GA_RND(); return val+low;
}

int RandomBit();

inline bool FlipCoin(float p){
    return((p == 1.0) ? true : (p == 0.0) ? false :
	   ((RandomFloat() <= p) ? true : false));
}

#endif	// __RANDOM_H_
