#ifndef INCLUDED_NUMBERS
#define INCLUDED_NUMBERS 1
#define USEBITS 28
#include <iostream>
#include <cmath>
namespace bignums
{
  static const unsigned long long FIRSTBIT=1,ZEROBITS=0,LASTBIT=(1<<(USEBITS-1)),HEAD=(0xF<<USEBITS-4),TAIL=0xF,FULLBITS=0xFFFFFFF,OVERFLOW=(1<<(USEBITS));
  //static const unsigned long long ZERO=0,ONE=1,TWO=2,THREE=3,FOUR=4,FIVE=5,SIX=6,SEVEN=7,EIGHT=8,NINE=9,TEN=10,ELEVEN=11,TWELVE=12,THIRTEEN=13,FOURTEEN=14,FIFTEEN=15,SIXTEEN=16;
  enum Digit {ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,TEN,ELEVEN,TWELVE,THIRTEEN,FOURTEEN,FIFTEEN};
  class BigInt
  {
    unsigned long long* number=nullptr;
  public:
    int sz=0;
    bool nega=false;
  public:
    BigInt() : BigInt(ZEROBITS) {}
    BigInt(unsigned long long n)
    {
      number = new unsigned long long[1];
      number[0]=n;
      sz=1;
    }
    BigInt(unsigned long long* n,int cs=1)
    {
      number = new unsigned long long[cs];
      sz=cs;
      for(int i=0;i<cs;i++) number[i]=n[i];
    }
    BigInt(const BigInt& bi)
    {
      number=new unsigned long long[bi.sz];
      for(int i=0;i<sz;i++) bi.number[i]=number[i];
      sz=bi.sz;
    }
    //~BigInt() {delete[] number;}

    unsigned long long operator[](int n) const {return number[n];}
    BigInt operator-() const
    {
      BigInt ret(*this);
      ret.nega=!nega;
      return ret;
    }
  };
}
namespace simple
{
  template<class T=int,class E=double> class Fraction
  {
    T p,q;
  public:
    Fraction() : Fraction(0,1) {}
    Fraction(E n,bool autosimp=true,int b=10,int p=-3);

    Fraction(const T& n,const T& d) {p=n; q=d;}

    inline const T& getNumerator() const {return p;}
    inline const T& getDenominator() const {return q;}
    inline operator E() const {return ((E)p)/((E)q);}

    inline void setNumerator(T n) {p=n;}
    inline void setDenominator(T n) {q=n;}
  };
  static void simplify(Fraction<int,double>& f,bool force=false,int p=-3)
  {
    if(p>0) p=-p;
    if(f.getNumerator()==0) {f.setDenominator(1); return;}
    if(force)
    {
      double number=f;
      double sc=::pow(10,-p);
      int num=(int)(number*sc);
      int den=(int)sc;
      f.setNumerator(num); f.setDenominator(den);
      simplify(f,false);
    }
    bool upd=true;
    int num=f.getNumerator(),den=f.getDenominator();
    if(num%den==0) {num/=den; den=1;}
    while(upd)
    {
      upd=false;
      int stop=den/2+1;
      for(int i=2;i<=stop;i++) {while(num%i==0 &&  den%i==0) {num/=i; den/=i; upd=true;} if(upd) break;}
    }
    //if(num%den==0) num/=den; den=1;
    f.setNumerator(num); f.setDenominator(den);
  }
}
template<class T,class E> simple::Fraction<T,E>::Fraction(E n,bool autosimp,int b,int p)
{
  if(p>0) p=-p;
  double sc=::pow(b,-p);
  E eff=n*sc;
  T trunceff=(T)eff;
  this->p=trunceff;
  q=T(sc);
  if(autosimp) simplify(*this);
}
static std::ostream& operator<<(std::ostream& os,const bignums::BigInt& n)
{
  std::string hex="";
  unsigned int hold=0;
  for(int i=n.sz-1;i>=0;i--)
  {
    unsigned long long bits=n[i];
    for(int i=0;i<USEBITS/4;i++)
    {
      hold=(bits>>USEBITS-4)&(bignums::TAIL);
      bits<<=4;
      switch(hold)
      {
        case 0:
          if(hex!="") hex+="0";
          break;
        case 10: hex+='A'; break;
        case 11: hex+='B'; break;
        case 12: hex+='C'; break;
        case 13: hex+='D'; break;
        case 14: hex+='E'; break;
        case 15: hex+='F'; break;
        default: hex+=std::to_string(hold); break;
      }
    }
  }
  if(hex=="") hex="0";
  return (os << hex);
}
static bignums::BigInt operator+(const bignums::BigInt& n1,const bignums::BigInt& n2)
{
  if(n2.sz>n1.sz) return n2+n1;
  int nsz=n1.sz;
  //if(n1.sz==n2.sz && (n1[n1.sz-1]&bignums::LASTBIT) && (n2[n2.sz-1]&bignums::LASTBIT)) nsz+=1;
  unsigned long long* retv=new unsigned long long[nsz+1]; retv[nsz]=0;
  bool carry=false;
  for(int i=0;i<n1.sz;i++)
  {
    if(i>=n2.sz)
    {
      if(carry && n1[i]==bignums::FULLBITS) {retv[i]=0; carry=true;}
      else
      {
        if(carry) {retv[i]=n1[i]+1;}
        else retv[i]=n1[i];
        carry=false;
      }
      continue;
    }
    long long int p1=n1[i],p2=n2[i];
    if(carry)
    {
      if(p1==bignums::FULLBITS && p2==bignums::FULLBITS) {carry=true; retv[i]=bignums::FULLBITS; continue;}
      else if(p1!=bignums::FULLBITS) p1+=1;
      else p2+=1;
    }
    //if((p1&bignums::LASTBIT) && (p2&bignums::LASTBIT)) {carry=true; p1&=!bignums::LASTBIT; p2&=!bignums::LASTBIT;}
    retv[i]=p1+p2;
    carry=retv[i]&bignums::OVERFLOW;
  }
  if(carry) retv[nsz]+=1;
  if(retv[nsz]) return bignums::BigInt(retv,n1.sz+1);
  else return bignums::BigInt(retv,n1.sz);
}
static bignums::BigInt operator-(const bignums::BigInt& n1,const bignums::BigInt& n2)
{
  return n1+(-n2);
}
typedef simple::Fraction<int,double> Frac;
template<class T,class E> inline static std::ostream& operator<<(std::ostream& os,const simple::Fraction<T,E>& f) {return (f.getDenominator()!=1)?(os << f.getNumerator()<<"/"<<f.getDenominator()):(os << f.getNumerator());}
template<class T,class E> inline static simple::Fraction<T,E> operator*(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return simple::Fraction<T,E>(f1.getNumerator()*f2.getNumerator(),f1.getDenominator()*f2.getDenominator());}
template<class T,class E> inline static simple::Fraction<T,E> operator/(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return simple::Fraction<T,E>(f1.getNumerator()*f2.getDenominator(),f1.getDenominator()*f2.getNumerator());}
template<class T,class E> inline static simple::Fraction<T,E> operator+(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return simple::Fraction<T,E>(f1.getNumerator()*f2.getDenominator()+f1.getDenominator()*f2.getNumerator(),f1.getDenominator()*f2.getDenominator());}
template<class T,class E> inline static simple::Fraction<T,E> operator-(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return simple::Fraction<T,E>(f1.getNumerator()*f2.getDenominator()-f1.getDenominator()*f2.getNumerator(),f1.getDenominator()*f2.getDenominator());}
template<class T,class E> inline static simple::Fraction<T,E> operator>(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()>f2.getNumerator()*f1.getDenominator();}
template<class T,class E> inline static simple::Fraction<T,E> operator>=(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()>=f2.getNumerator()*f1.getDenominator();}
template<class T,class E> inline static simple::Fraction<T,E> operator<(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()<f2.getNumerator()*f1.getDenominator();}
template<class T,class E> inline static simple::Fraction<T,E> operator<=(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()<=f2.getNumerator()*f1.getDenominator();}
template<class T,class E> inline static simple::Fraction<T,E> operator==(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()==f2.getNumerator()*f1.getDenominator();}
template<class T,class E> inline static simple::Fraction<T,E> operator!=(const simple::Fraction<T,E>& f1,const simple::Fraction<T,E>& f2) {return f1.getNumerator()*f2.getDenominator()!=f2.getNumerator()*f1.getDenominator();}
#endif
