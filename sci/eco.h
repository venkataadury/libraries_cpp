#include "maths/maths.h"

namespace eco
{
  //typedef double EnvironmentVariable;
  class EcoSystem
  {
    double scale=1;
    int season=0;

  public:
    EcoSystem() {}

  };
  typedef EcoSystem Ecosystem;
  class EnvironmentVariable
  {
    double ov=0;
    EcoSystem es;
    bool system=false;
  public:
    EnvironmentVariable() {}
    EnvironmentVariable(const double& d) {ov=d;}
    EnvironmentVariable(const double& d,const EcoSystem& e) : EnvironmentVariable(d) {es=e; system=true;}

    virtual double getModified(const EcoSystem& es) {return ov;}
    void setEcosystem(const EcoSystem& e) {es=e; system=true;}
    double getBaseValue() {return ov;}

    //Operator overloading
    double& operator=(const double& v) {ov=v; return ov;}
    operator double() {if(system) return getModified(es); else return ov;}
  };
  class LeslieMatrix
  {
  public:
    Matrix P,G;
    Vector F;
    int size=0;
    Matrix L;

  public:
    LeslieMatrix() {P=Matrix(),F=Vector();}
    LeslieMatrix(int s) {size=s; P=Matrix(s,s); F=Matrix(s,s); G=Matrix(s,s);}
    LeslieMatrix(const Vector& f)
    {
      size=f.getSize();
      P=Matrix(size,size);
      F=Matrix(size,size);
      for(int i=0;i<size;i++)
        F[0][i]=f[i][0];
    }
    LeslieMatrix(const Vector& f,const Vector& p) : LeslieMatrix(f)
    {
      for(int i=0;i<size-1;i++)
        P[i+1][i]=p[i][0];
      pack();
    }

    inline void pack() {L=P+F+G;}

    MatrixRow& operator[](int x) {return P[x];}
    const MatrixRow& operator[](int x) const {return P[x];}

  };
  class GrowthModel
  {
  protected:
    LeslieMatrix LM;
  public:
    virtual double growthRate(const Vector& pop) {return 1;} //Per-Capita growth rate
  };
  class LogisticModel : public GrowthModel
  {
    EnvironmentVariable r;
    EnvironmentVariable k;
  public:
    LogisticModel(double rf=1,unsigned long int kv=1000) {r=rf;k=kv;}
    double growthRate(const Vector& pop) {return r*(1-sum(pop)/k);}
  };
}
