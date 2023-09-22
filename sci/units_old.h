#ifndef UNITS_INCLUDED
#define UNITS_INCLUDED 1
#include "commons/commons.h"

class DerivedUnit;
class PrimitiveUnit;
class Unit
{
public:
  double v;
  double pow=1;
  std::vector<const Unit*> basis;

public:
  virtual double toSI(const double& v) const =0;
  virtual double fromSI(const double& v) const =0;

  inline double getVal() const {return v;}
  inline double getValSI() const {return toSI(v);}
  inline DerivedUnit& reciprocal() const;

  inline void setFromSI(const double& val) {v=fromSI(val);}
  inline void setDefault(const double& val) {v=val;}

  //Operator overloading
  DerivedUnit operator*(const Unit& du) const;
  DerivedUnit operator/(const Unit& du) const;
  operator double() const ;
};

class PrimitiveUnit : public Unit
{
public:
  PrimitiveUnit(double n=1) {pow=n; basis.push_back(this);}
  PrimitiveUnit(const Unit& pu) {pow=pu.pow; v=pu.v; basis.push_back(this);}

  double toSI(const double& d) const
  {
    if(pow<0) return fromSIInt(d);
    else return toSIInt(d);
  }
  double toSIInt(const double& d) const
  {
    double r=d;
    for(int i=0;i<::abs(pow);i++)
      r=conversion(r);
    return r;
  }

  double fromSI(const double& d) const
  {
    if(pow>=0) return fromSIInt(d);
    else return toSIInt(d);
  }
  double fromSIInt(const double& d) const
  {
    double r=d;
    for(int i=0;i<::abs(pow);i++)
      r=anticonversion(r);
    return r;
  }
  virtual inline double conversion(const double& d) const =0;
  virtual inline double anticonversion(const double& d) const =0;

  //void operator+=(double n) {pow+=n;}
  double getPower() {return pow;}

};

class DerivedUnit : public Unit
{
public:
  explicit DerivedUnit(const Unit& u,double val=0) {basis.push_back(&u);v=val;}
  DerivedUnit(const Unit& u1,const Unit& u2,double val=0) : DerivedUnit(u1) {basis.push_back(&u2);v=val;}
  DerivedUnit(const std::vector<const Unit*>& uns,double val=0) {basis=uns;v=val;}
  DerivedUnit(const DerivedUnit& d) {basis=d.basis; v=d.v; pow=d.pow;}

  double toSI(const double& d) const override
  {
    if(pow<0) return fromSIInt(d);
    else return toSIInt(d);
  }
  double toSIInt(const double& v) const
  {
    double r=v;
    for(int i=0;i<::abs(pow);i++)
    {
      for(const Unit* u : basis)
        r=u->toSI(r);
    }
    return r;
  }
  double fromSI(const double& d) const override
  {
    if(pow>=0) return fromSIInt(d);
    else return toSIInt(d);
  }
  double fromSIInt(const double& v) const
  {
    double r=v;
    for(int i=0;i<::abs(pow);i++)
    {
      for(const Unit* u : basis)
        r=u->fromSI(r);
    }
    return r;
  }
};
namespace atomicunits
{
  //This namespace takes real-world (S.I) values as input, and (converts and) stores them in Atomic units.
  const class Length : public PrimitiveUnit //Angstroms
  {
  public:
    Length(double val=0) : Length(val,1) {}
    Length(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Length(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*1e-10;}
    inline double anticonversion(const double& d) const {return d*1e+10;}
  } LENGTH;

  const class Mass : public PrimitiveUnit // amu
  {
  public:
    Mass(double val=0) : Mass(val,1) {}
    Mass(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Mass(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*1.67e-27;}
    inline double anticonversion(const double& d) const {return d/1.67e-27;}
  } MASS;

  const class Time : public PrimitiveUnit // ms
  {
  public:
    Time(double val=0) : Time(val,1) {}
    Time(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Time(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d;}
    inline double anticonversion(const double& d) const {return d;}
  } TIME;
  const class Temperature : public PrimitiveUnit // K
  {
  public:
    Temperature(double val=0) : Temperature(val,1) {}
    Temperature(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Temperature(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d;}
    inline double anticonversion(const double& d) const {return d;}
  } TEMPERATURE;
  const class Charge : public PrimitiveUnit // Charge on electron
  {
  public:
    Charge(double val=0) : Charge(val,1) {}
    Charge(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Charge(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*1.602e-19;}
    inline double anticonversion(const double& d) const {return d/1.602e-19;}
  } CHARGE;
  const class Moles : public PrimitiveUnit // atoms
  {
  public:
    Moles(double val=0) : Moles(val,1) {}
    Moles(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Moles(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d/6.023e+23;}
    inline double anticonversion(const double& d) const {return d*6.023e+23;}
  } MOLES;
  //Basic Derived Units
  const Length AREA(0,2),VOLUME(0,3);

  //All derived units
  const class Speed : public DerivedUnit
  {
  public:
    Speed() : DerivedUnit(LENGTH,TIME.reciprocal()) {}
    Speed(double val) : Speed() {v=val;}
    Speed(double val,double p) : Speed(val) {pow=p;}
    Speed(const DerivedUnit& d) : DerivedUnit(d) {}
  } SPEED;
  const class Acceleration : public DerivedUnit
  {
  public:
    Acceleration() : DerivedUnit(SPEED,TIME.reciprocal()) {}
    Acceleration(double val) : Acceleration() {v=val;}
    Acceleration(const DerivedUnit& d) : DerivedUnit(d) {}
  } ACCELERATION;
  const class Force : public DerivedUnit
  {
  public:
    Force() : DerivedUnit(ACCELERATION,MASS) {}
    Force(double val) : Force() {v=val;}
    Force(const DerivedUnit& d) : DerivedUnit(d) {}
  } FORCE;
  const class Energy : public DerivedUnit
  {
  public:
    Energy() : DerivedUnit(FORCE,LENGTH) {}
    Energy(double val) : Energy() {v=val;}
    Energy(const DerivedUnit& d) : DerivedUnit(d) {}
  } ENERGY;
  const class Density : public DerivedUnit
  {
  public:
    Density() : DerivedUnit(MASS,VOLUME.reciprocal()) {}
    Density(double val) : Density() {setFromSI(val);}
    Density(const DerivedUnit& d) : DerivedUnit(d) {}
  } DENSITY;
  const class Pressure : public DerivedUnit
  {
  public:
    Pressure() : DerivedUnit(FORCE,AREA.reciprocal()) {}
    Pressure(double val) : Pressure() {setFromSI(val);}
    Pressure(const DerivedUnit& d) : DerivedUnit(d) {}
  } PRESSURE;

  //Free derived unit constants
  const class FreeDerived1 : public DerivedUnit
  {
  public:
    FreeDerived1() : DerivedUnit(ENERGY,TEMPERATURE.reciprocal()) {}
    FreeDerived1(double val) : FreeDerived1() {setFromSI(val);}
    FreeDerived1(const DerivedUnit& d) : DerivedUnit(d) {}
  } Kb(1.38e-23);
  const class FreeDerived2 : public DerivedUnit
  {
  public:
    FreeDerived2() : DerivedUnit(MOLES.reciprocal()) {}
    FreeDerived2(double val) : FreeDerived2() {setFromSI(val);}
    FreeDerived2(const DerivedUnit& d) : DerivedUnit(d) {}
  } NA(6.023e+23);
  const class FreeDerived3 : public DerivedUnit
  {
  public:
    FreeDerived3() : DerivedUnit(ENERGY,TIME) {}
    FreeDerived3(double val) : FreeDerived3() {setFromSI(val);}
    FreeDerived3(const DerivedUnit& d) : DerivedUnit(d) {}
  } hcut(0.527e-34);
}
namespace thermodynamicunits
{
  const class Length : public PrimitiveUnit //decimeters (because volume is dm^3)
  {
  public:
    Length(double val=0) : Length(val,1) {}
    Length(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Length(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*0.1;}
    inline double anticonversion(const double& d) const {return d*10;}
  } LENGTH;
  const class Area : public DerivedUnit
  {
  public:
    Area() : DerivedUnit(LENGTH,LENGTH) {}
    Area(double val) : Area() {v=val;}
    Area(const DerivedUnit& d) : DerivedUnit(d) {}
  } AREA(1);
  const class Volume : public DerivedUnit
  {
  public:
    Volume() : DerivedUnit(AREA,LENGTH) {}
    Volume(double val) : Volume() {setFromSI(val);}
    Volume(const DerivedUnit& d) : DerivedUnit(d) {}
  } VOLUME(1);
  //Basic Derived Units

  const class Temperature : public PrimitiveUnit // K
  {
  public:
    Temperature(double val=0) : Temperature(val,1) {}
    Temperature(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Temperature(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d;}
    inline double anticonversion(const double& d) const {return d;}
  } TEMPERATURE;

  const class Moles : public PrimitiveUnit // moles
  {
  public:
    Moles(double val=0) : Moles(val,1) {}
    Moles(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Moles(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d;}
    inline double anticonversion(const double& d) const {return d;}
  } MOLES;

  const class Mass : public PrimitiveUnit // amu
  {
  public:
    Mass(double val=0) : Mass(val,1) {}
    Mass(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Mass(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*1.67e-27;}
    inline double anticonversion(const double& d) const {return d/1.67e-27;}
  } MASS;

  const class Pressure : public PrimitiveUnit // atm
  {
  public:
    Pressure(double val=0) : Pressure(val,1) {}
    Pressure(double val,double p) : PrimitiveUnit(p) {setFromSI(val);}
    Pressure(const Unit& pu) : PrimitiveUnit(pu) {}

    inline double conversion(const double& d) const {return d*1.01325e+5;}
    inline double anticonversion(const double& d) const {return d*9.86923e-06;}
  } PRESSURE;
  const class MolarPressure : public DerivedUnit
  {
  public:
    MolarPressure() : DerivedUnit(PRESSURE,MOLES.reciprocal()) {}
    MolarPressure(double val) : MolarPressure() {setFromSI(val);}
    MolarPressure(const DerivedUnit& d) : DerivedUnit(d) {}
  } MOLARPRESSURE;
  const class MolarVolume : public DerivedUnit
  {
  public:
    MolarVolume() : DerivedUnit(VOLUME,MOLES.reciprocal()) {}
    MolarVolume(double val) : MolarVolume() {setFromSI(val);}
    MolarVolume(const DerivedUnit& d) : DerivedUnit(d) {}
  } MOLARVOLUME;
  const class Energy : public DerivedUnit
  {
  public:
    Energy() : DerivedUnit(PRESSURE,VOLUME) {}
    Energy(double val) : Energy() {setFromSI(val);}
    Energy(const DerivedUnit& d) : DerivedUnit(d) {}
  } ENERGY;
  typedef Energy Work;
  const Work WORK(0);
  const class Entropy : public DerivedUnit
  {
  public:
    Entropy() : DerivedUnit(ENERGY,TEMPERATURE.reciprocal()) {}
    Entropy(double val) : Entropy() {setFromSI(val);}
    Entropy(double val,double p) : Entropy(val) {pow=p;}
    Entropy(const DerivedUnit& d) : DerivedUnit(d) {}
  } ENTROPY;

  const class FreeDerived2 : public DerivedUnit
  {
  public:
    FreeDerived2() : DerivedUnit(MOLES.reciprocal()) {}
    FreeDerived2(double val) : FreeDerived2() {setFromSI(val);}
    FreeDerived2(const DerivedUnit& d) : DerivedUnit(d) {}
  } NA(1);
  const class FreeDerived3 : public DerivedUnit
  {
  public:
    FreeDerived3() : DerivedUnit(ENERGY,TEMPERATURE.reciprocal()) {}
    FreeDerived3(double val) : FreeDerived3() {setFromSI(val);}
    FreeDerived3(const DerivedUnit& d) : DerivedUnit(d) {}
  } R_gas(8.314);
  class FreeDerived4 : public DerivedUnit
  {
  public:
    FreeDerived4() : DerivedUnit(PRESSURE,VOLUME.reciprocal()) {}
    FreeDerived4(double val) : FreeDerived4() {setFromSI(val);}
    FreeDerived4(const DerivedUnit& d) : DerivedUnit(d) {}
  };

  Temperature room_temp(298);
  Temperature zero(273.15);
  Temperature abszero(0);
  Pressure atmos_press(1.01325e+5);
}

//Missing Methods
inline DerivedUnit& Unit::reciprocal() const
{
  DerivedUnit* du=new DerivedUnit(*this);
  du->pow*=-1;
  return *du;
}
DerivedUnit Unit::operator*(const Unit& du) const
{
  /*std::vector<const Unit*> der=basis;
  for(const Unit* u : du.basis)
    der.push_back(u);*/
  DerivedUnit dur(*this,du);
  dur.v=v*du.v;
  return dur;
}
DerivedUnit Unit::operator/(const Unit& du) const
{
  std::vector<const Unit*> der=basis;
  for(const Unit* u : du.basis)
    der.push_back(&(u->reciprocal()));
  DerivedUnit dur(der);
  dur.v=v/du.v;
  return dur;
}

//operator overloading
Unit::operator double() const {return getVal();}

static DerivedUnit operator/(double n,const Unit& u)
{
  DerivedUnit der(u.reciprocal());
  der.v=n/u.v;
  return der;
}
static DerivedUnit operator/(const Unit& u,double n)
{
  DerivedUnit der(u);
  der.v=u.v/n;
  return der;
}
static DerivedUnit operator*(const Unit& u,double n)
{
  DerivedUnit der(u);
  der.v=u.v*n;
  return der;
}
inline static DerivedUnit operator-(const Unit& u) {return u*-1;}

static DerivedUnit operator-(const Unit& u1,const Unit& u2)
{
  return DerivedUnit(u1,u1.v-u2.v);
}
static DerivedUnit operator+(const Unit& u1,const Unit& u2)
{
  return DerivedUnit(u1,u1.v+u2.v);
}
#endif
