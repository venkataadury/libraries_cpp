#ifndef INCLUDED_VECMATHS
#define INCLUDED_VECMATHS 1
//#include "maths/maths.h"

namespace fnvsmath //field n vector space math
{
  //Fields are any classes which +,* defined on them, and have ZERO,ONE elements
  class DoubleField
  {
    double v;
    static DoubleField *ZERO,*ONE;
  public:
    DoubleField() : DoubleField(0) {}
    DoubleField(double d) {v=d;}
    operator double() const {return v;}
    void operator+=(const DoubleField& df) {v+=df.v;}
    void operator-=(const DoubleField& df) {v-=df.v;}
    void operator*=(const DoubleField& df) {v*=df.v;}
    void operator/=(const DoubleField& df) {v/=df.v;}

    static const DoubleField& getZero() {return *ZERO;}
    static const DoubleField& getOne() {return *ONE;}
  };
  DoubleField* DoubleField::ZERO=new DoubleField(0);
  DoubleField* DoubleField::ONE=new DoubleField(1);
  /*template<class F> class VectorSpace
  {
    VectorSpace<F> operator+(const VectorSpace<F>& el) const =0;
    VectorSpace<F> operator-() const =0;
    inline VectorSpace<F> operator-(const VectorSpace<F>& el) {return operator+(el.operator-());}
    VectorSpace<F> operator*(const F& el) const =0;
  };*/
}
using namespace fnvsmath;
namespace vcalc
{

}
#endif
