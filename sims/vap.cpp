#include "maths/maths.h"
#include "sci/sim.h"

using namespace simtools;

double a=10,V=a*a*a;
double T=298,k=1.38064e-3,NA=6.023e+23,PI=3.14159265359,BETA=1/(k*T),DX=0.4,Dr=1e-4,STEPD=0.1; //Lengths in A
const Vector DXV(Dr,0,0),DYV(0,Dr,0),DZV(0,0,Dr);
//Water
double SIGMA=2.725,EPSILON=0.49115,MASS=30.103e-27,DENSITY=0.032*2.652e-2;
//double SIGMA=3.166,EPSILON=78.2*k,MASS=30.103e-27,DENSITY=0.03*2.652e-2;
double LAMBDA=6.626070e-14/sqrt(2*PI*MASS*k*T),LAMBDA3=LAMBDA*LAMBDA*LAMBDA;

//Large System: rho is atoms per volume
double P=1.01325e-5*0.025;
double rho=P/(k*T),MU=k*T*log(LAMBDA3*rho); //Volume of container = 1nm^3

inline double LJPot(double d) {return 4*EPSILON*(pow(SIGMA/d,12)-pow(SIGMA/d,6));}
double energy(const System& sys)
{
  double e=0,d;
  for(int i=0;i<sys.getParticleCount();i++)
  {
    for(int j=i+1;j<sys.getParticleCount();j++)
    {
      try
      {
        d=sys.distance(sys.getParticleAt(i),sys.getParticleAt(j));
        e+=LJPot(d);
      }
      catch(exception e)
      {
        cout << i << ' ' << j << "\t"<<sys.getParticleCount()<<"\n";
        cout <<"---"<< sys.getParticleAt(i).getPosition().getSize() <<"---\n"<<sys.getParticleAt(j).getPosition()<<"\n\nEND\n";
        assert(1!=1);
      }
    }
  }
  return e;
}

void mcmoveDesc(PBCBox& sys)
{
  Particle& ptemp=sys.getRandomParticle();
  Vector op=ptemp; //Automatically extracts position
  ptemp.setPosition(simtools::randVect(DX)+op);
  double en=energy(sys);
  //cout << en << " "<<sys.getEnergy() << "\n";
  if(en<=sys.getEnergy()) {sys.setEnergy(en); return;}
  else {ptemp.setPosition(op); return;}
}
void equilibriate(PBCBox& sys)
{
  //cout << sys.getParticles().size() <<"\n";
  for(int i=0;i<1000;i++)
    mcmoveDesc(sys);
}

void mcad(PBCBox& sys)
{
  if(maths::random()<0.5) //insert
  {
    Particle* np=&sys.addRandomParticle();
    double en=energy(sys); //Fails
    double p=(V/(LAMBDA3*(sys.getParticles().size()+1)))*exp(-BETA*(en-sys.getEnergy()-MU));
    if(maths::random()<=p) {sys.setEnergy(en); return;}
    else {sys.removeParticle(*np); return;}
  }
  else //remove
  {
    Particle& p=sys.removeRandomParticle();
    double en=energy(sys);
    double pr=((LAMBDA3*sys.getParticles().size())/V)*exp(-BETA*(MU+en-sys.getEnergy()));
    //cout << MU<<" "<<en<<" "<<sys.getEnergy()<<"\t"<<pr<<"\n";
    if(maths::random()<=pr) {sys.setEnergy(en); delete &p;}
    else sys.addParticle(p);

  }
}
void mcmove(PBCBox& sys)
{
  Particle& ptemp=sys.getRandomParticle();
  Vector op=ptemp; //Automatically extracts position
  //ptemp.setPosition(sys.randVect());
  ptemp.setPosition(simtools::randVect(DX)+op);
  double en=energy(sys);
  if(maths::random()<exp(-BETA*(-sys.getEnergy()+en))) {sys.setEnergy(en); return;}
  else {ptemp.setPosition(op); return;}
}
void mdmove(PBCBox& box)
{
  double vo,vdx,vdy,vdz;
  Vector RV;
  for(int i=0;i<box.getParticleCount();i++)
  {
    vo=0; vdx=0; vdy=0; vdz=0;
    for(int j=0;j<box.getParticleCount();j++)
    {
      if(j==i) continue;
      RV=box.getParticleAt(i);
      vo+=LJPot(box.distance(box.getParticleAt(i),box.getParticleAt(j)));
      vdx+=LJPot(box.distance(RV+DXV,box.getParticleAt(j)));
      vdy+=LJPot(box.distance(RV+DYV,box.getParticleAt(j)));
      vdz+=LJPot(box.distance(RV+DZV,box.getParticleAt(j)));
    }
    vdx=(vdx-vo)/(Dr*MASS); vdy=(vdy-vo)/(Dr*MASS); vdz=(vdz-vo)/(Dr*MASS);
    box.getParticleAt(i).setVelocty(-(Vector(vdx,vdy,vdz).unitVector()));
    box.getParticleAt(i).move(STEPD);
    //cout << "\t" << box.getParticleAt(i).getVelocity().getMagnitude() << "\n";
  }
}
void mainX(double x,int cd=-2)
{
  STEPD=0.18;
  P=x; rho=P*BETA;
  MU=k*T*log(LAMBDA3*rho);
  int Nw=(int)((DENSITY*V)/MASS);
  double moveP=0.05; //probability to move existing particle (not adding or removing)
  double en=0;
  int nmoves=7.5e3;

  //cout << LAMBDA << "\n"; <!-- Verified -->
  //cout << Nw <<"\n"; <!-- Verified -->
  //cout << "\t"<<MU<< "\n";


  PBCBox box(a);
  //for(int i=0;i<Nw;i++) box.addRandomParticle();
  for(double i=10/4.0-0.4;i<=10;i+=10/4.0)
  {
    for(double j=10/3.0-0.3;j<=10;j+=10/3.0)
    {
      for(double k=10/3.0-0.3;k<=10;k+=10/3.0)
        box.addParticle(Vector(i,j,k));
    }
  }
  box.setEnergy(energy(box));
  //cout << box.getEnergy() << "\n";
  equilibriate(box);

  cout <<"P="<<"\t"<< P<<"\n";
  int printfreq=nmoves/10,K=0;
  std::vector<double> nC;
  for(int i=0;i<nmoves;i++)
  {
    //if(i%printfreq==0) cout << i <<" "<<box.getParticleCount()<< "\t"<<box.getEnergy()<<"\n";
    box.setEnergy(energy(box));
    if(maths::random()<moveP)
    {
      mdmove(box);
      i--;
    }
    else
    {
      mcad(box);
      //if(i>printfreq)
      nC.push_back(box.getParticleCount());
      cout << K++ <<" "<< box.getParticleCount() << "\n";
      if(box.getParticleCount()<=1) break;
    }
  }
  double mN=stats::mean(nC);
  if(box.getParticleCount()<=1)
  {
    if(cd==2)
      return;
    else
      mainX(x,cd+1);
  }
  else
  {
    cout << P << "\t";
    cout << mN<<"\n";
    cin >> mN;
  }
  /*cout << rho<<","<<MU << "\t";
  cout << mN << "\t";
  cout << stats::stddev(nC) << "\n";*/
}
int main()
{
  /*cout << "\t"<<k*T<<"\t"<<LAMBDA3<<" "<<rho<<"->"<<(LAMBDA3*rho)<<" log:"<<log(LAMBDA3*rho)<<"\n";
  cout << MU << "\n";*/
  maths::resetRandSeed();
  //Expecting 7.41 mbar = 0.0075 +/- 0.0005 bar
  for(double p=1.01325e-6;p<=10*1.01325e-6;p+=1.01325e-6) //~0.5 +/- 0.1 bar for eq.
  {
    mainX(p);
  }
}
