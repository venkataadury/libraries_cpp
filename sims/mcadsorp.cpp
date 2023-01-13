#include "commons/commons.h"
#include "maths/maths.h"
#include "maths/stats.h"
#include "sci/sim.h"

using namespace simtools;
using namespace randgen;
using maths::Vector;
/** ALL ARE S.I. UNITS **/
//Conditions and constants
const double T=298,Kb=1.38e-23,NA=6.023e+23,R=Kb*NA,M=2,h=6.626e-34;
const double BETA=1/(Kb*T);
double P=101325;
//System
System Pt;
const double a=1e-7/*0.1mm*/,l=1.5e-10,V=6*a*a*l-(8*l*l*l),cutoff=l-0.6e-10;
double LAMBDA3=3.6E-31,MU=Kb*T*(log(LAMBDA3*2*P/(R*T))),E=577392/(NA*1000*2.71); /*J/atom*/ //Need to fix?


inline Vector getRandomVector() {return Vector(a*maths::random(),a*maths::random(),l*maths::random());}
double overlap(const Vector& v1,const Vector& v2)
{
  double dis=Vector(v1-v2).getMagnitude();
  if(dis<3.82e-9) /*Initially e-10*///2*rH + H-H sigma length
    return E*exp(4e-9/dis);
  else
    return 0;
}
double energyOf(const System& sys,Particle ext)
{

  if(ext.getZ()>cutoff) {return 1;}
  double e=sys.getEnergy()-E*exp(1-ext.getZ()/l);
  double o;
  for(const Particle& p : sys.getParticles())
  {
    o=overlap(p.getPosition(),ext.getPosition());
    if(o!=0)
    {
      e+=o;
      break;
    }
  }
  return e;
}

//Simulation
void addParticle(System& sys)
{
  Particle p(getRandomVector());
  double de=energyOf(sys,p)-sys.getEnergy();
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
  double de=E;
  double o;
  for(Particle& sp : sys.getParticles())
  {
    if(&sp==&p) continue;
    o=overlap(p.getPosition(),sp.getPosition());
    if(o!=0)
    {
      de-=o;
      break;
    }
  }
  //double ap=((LAMBDA3*sys.getParticleCount())/V)*exp(BETA*(MU-de));
  double ap=((sys.getParticleCount())/(V*BETA*P))*exp(BETA*(MU-de));
  //cout << "\t" << ap << "\n";
  if(maths::random()<=ap) //Accepted
  {
    sys.removeParticle(p);
    sys.setEnergy(sys.getEnergy()+de);
  }
}

const int STEPS=1e5;
int main()
{
  for(int i=1;i<=3000;i+=100)
  {
    P=i*101325;
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
