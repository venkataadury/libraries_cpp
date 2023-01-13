#include "commons/commons.h"
#include "maths/maths.h"
#include "maths/stats.h"
#include "sci/sim.h"
#include "sci/units.h"

using namespace simtools;
using namespace randgen;
using maths::Vector;
using namespace atomicunits;
/** ALL ARE S.I. UNITS **/
//Conditions and constants
const Temperature T=Temperature(298);
const DerivedUnit R=Kb*NA;
const DerivedUnit BETA=1/(Kb*T);
Pressure P(101325);
//System
System Pt;
const Length a(1e-7)/*0.1mm*/,l(1.5e-10),cutoff(l-0.6e-10);
DerivedUnit V=a*a*l*6-(l*l*l*8);
DerivedUnit LAMBDA3(VOLUME,3.6E-31);
DerivedUnit MU=Kb*T*(log(LAMBDA3*2*P/(R*T)));
Energy E=577392/(NA*1000*2.71); /*J/atom*/ //Need to fix?


inline Vector getRandomVector() {return Vector(a*maths::random(),a*maths::random(),l*maths::random());}
Energy overlap(const Vector& v1,const Vector& v2)
{
  double dis=Vector(v1-v2).getMagnitude();
  if(dis<3.82e-9) /*Initially e-10*///2*rH + H-H sigma length
    return E*exp(4e-9/dis);
  else
    return 0;
}
Energy energyOf(const System& sys,Particle ext)
{

  if(ext.getZ()>cutoff) {return 1;}
  Energy e=sys.getEnergy()-E*exp(1-ext.getZ()/l);
  Energy o;
  for(const Particle& p : sys.getParticles())
  {
    o=overlap(p.getPosition(),ext.getPosition());
    if(o!=0)
    {
      e=e+o;
      break;
    }
  }
  return e;
}

//Simulation
void addParticle(System& sys)
{
  Particle p(getRandomVector());
  Energy de=energyOf(sys,p)-sys.getEnergy();
  //double ap=(V/(LAMBDA3*(sys.getParticleCount()+1)))*exp(BETA*(MU-de));
  double ap=((V*BETA*P)/((sys.getParticleCount()+1)))*exp(BETA*(MU-de));
  //cout << MU<<" "<<de << " "<<(MU-de)<<"\t"<<exp(BETA*(MU-de)) << "\t"<<(V*BETA*P)<<"\n";
  //cout << ap <<"\n";
  if(maths::random()<=ap) //Accepted
  {
    sys.addParticle(p);
    sys.setEnergy(sys.getEnergy()+de);
  }
}
void removeParticle(System& sys)
{
  if(sys.getParticleCount()<1) return;
  Particle& p=sys.getRandomParticle();
  Energy de=E;
  double o;
  for(Particle& sp : sys.getParticles())
  {
    if(&sp==&p) continue;
    o=overlap(p.getPosition(),sp.getPosition());
    if(o!=0)
    {
      de=de-o;
      break;
    }
  }
  //double ap=((LAMBDA3*sys.getParticleCount())/V)*exp(BETA*(MU-de));
  double ap=((sys.getParticleCount())/(V*BETA*P))*exp(BETA*(MU-de));
  //cout << "\t" << ap << "\n";
  if(maths::random()<=ap) //Accepted
  {
    sys.removeParticle(p);
    sys.setEnergy(Energy(sys.getEnergy())+de);
  }
}

const int STEPS=1e5;
int main()
{
  for(int i=1;i<=3000;i+=100)
  {
    P.setFromSI(i*101325);
    MU=Kb*T*(log(LAMBDA3*P/(Kb*T)));
    Pt=System();
  std::vector<double> nParts;
  for(int i=0;i<STEPS;i++)
  {
    if(maths::random()<0.5)
      addParticle(Pt);
    else
      removeParticle(Pt);
    nParts.push_back(Pt.getParticleCount());
  }
  cout <<i<< " "<<stats::mean(nParts)<<"\n";
  }
}
