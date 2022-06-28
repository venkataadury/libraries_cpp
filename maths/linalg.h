#ifndef INCLUDED_LINALG
#define INCLUDED_LINALG 1
#include "maths/maths.h"
using maths::Vector;
using maths::Matrix;

//Basic Functions
inline bool hasAcuteAngle(const Vector& v1,const Vector& v2) {return v1.dot(v2)>=0;}
inline bool isParallel(const Vector& v1,const Vector& v2,double cut=1e-6) {return ::abs(::abs(v1.dot(v2)/(v2.getMagnitude()*v1.getMagnitude()))-1)<cut;}
bool isInSpace(const std::vector<Vector>& vecs,const Vector& v)
{
  Vector rem=v;
  for(const Vector& vc : vecs)
  {
    rem=rem-(vc.dot(v)/vc.getMagnitude())*vc;
    if(rem.isZero()) return true;
  }
  return false;
}

#endif INCLUDED_LINALG
