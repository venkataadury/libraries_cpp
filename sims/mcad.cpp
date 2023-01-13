#include "sci/sim.h"
#include "maths/stats.h"
using maths::Vector;
using namespace simtools;

const double ETAS=148.2*1.38e-3,k=1.38e-3,T=300;
const double SIGMA=3.817,ES=148.2*k,HSC=SIGMA/3.5,SEIVER=7.01;
const double a=SEIVER-HSC,PI=3.141592,DX=1.75;
const double LAMS=1.93*k,VS=33.49,V=SEIVER*SEIVER*SEIVER/sqrt(2);

double EX=0;
double uwp=0,uljp=0;
int accM=0,n=5;

inline static double l(double x) {return (1+12*x+25.2*x*x+12*x*x*x+pow(x,4))*pow(1-x,-10.0);}
inline static double m(double x) {return (1+x)*pow(1-x,-4.0);}
inline static double wallPot(double x)
{
  //cout << pow(VS/V,4)*l((x*x)/(a*a)) << "," << 2*(pow(VS/V,2))*m((x*x)/(a*a)) << "\t";
  if(x>a) return maths::inf();
  return LAMS*(pow(VS/V,4)*l((x*x)/(a*a))-2*(pow(VS/V,2))*m((x*x)/(a*a)));
}
bool HSP(const std::vector<Particle>& sys)
{
  Vector v1;
  for(int i=0;i<sys.size();i++)
  {
    for(int j=i+1;j<sys.size();j++)
    {
      if(Vector(sys[i].getPosition()-sys[j].getPosition()).getMagnitude()<HSC)
        return false;
    }
  }
  return true;
}
double LJPot(const std::vector<Particle>& sys)
{
  double ljp=0,d,wp=0;
  for(int i=0;i<sys.size();i++)
  {
    for(int j=i+1;j<sys.size();j++)
    {
      d=(Vector(sys[i].getPosition()-sys[j].getPosition())).getMagnitude();
      if(d<HSC) return maths::inf();
      d=4*ETAS*(pow(SIGMA/d,12)-pow(SIGMA/d,6));
      //cout <<"\t"<<d<<" "<< ::pow(SIGMA/d,12) << "\n";
      ljp+=d;
    }
    wp+=wallPot(sys[i].getPosition().getMagnitude());
  }
  /*if(wp>0)
    cout << wp << "\n";*/
  uljp+=ljp; uwp+=wp;
  ljp+=wp;
  return EX*ljp;
}
inline double energy(const std::vector<Particle>& sys) {return LJPot(sys);}

std::vector<Particle> randomConfig(int n)
{
  std::vector<Particle> ret;
  for(int i=0;i<n;i++)
    ret.push_back(Particle(a*maths::random()*randVect()));
  return ret;
}
std::vector<Particle> randomConfig() {return randomConfig(n);}
std::vector<Particle> randomConfigInclusive(int n)
{
  std::vector<Particle> ret;
  for(int i=0;i<n;i++)
    ret.push_back(Particle(SEIVER*maths::random()*randVect()));
  return ret;
}
std::vector<Particle> randomConfigInclusive() {return randomConfigInclusive(n);}
double mcmove(std::vector<Particle>& conf,double& en)
{
  Particle& pt=randgen::choice(conf);
  Vector mv=DX*maths::random()*randVect();
  pt.pos=pt.pos+mv;
  double e=LJPot(conf);
  /*if(e>1e5)
    cout << e << "\n";*/
  //cout << exp((en-e)/(k*T)) << "\n";
  /*if(e>1e5)
    cout << e <<" "<<en<< "\t" << (en-e)<<"\t"<< exp((en-e)/(k*T)) << "\n";*/
  if(maths::random()<=exp((en-e)/(k*T)))
  {
    accM++;
    return e;
  }
  else
  {
    pt.pos=pt.pos-mv;
    return en;
  }
}
int nmoves=1e5;
int main()
{
  maths::resetRandSeed();
  n=5;
  //auto cG=[n]() {return randomConfig(n);}
  double ef;
  for(EX=0.3;EX<=3;EX+=0.3)
  {
    ef=MonteCarlo::importance(randomConfig,energy,1/(k*T),DX);
    cout <<(int)(EX/0.29999)<<" "<<ef<<"\t";
    ef=MonteCarlo::bruteForce(randomConfig,energy,1/(k*T));
    cout <<ef<<"\n";
  }
}
int mainolder()
{
  maths::resetRandSeed();
  int nmoves=1e5;
  n=5;
  double Z0;
  std::vector<Particle> conf;
  EX=0;
  /*int acc=0;
  double f;
  for(int i=0;i<nmoves;i++)
  {
    conf=randomConfig(n);
    if(HSP(conf))
      acc++;
  }
  f=(double)acc/nmoves;
  Z0=pow((4/3.0)*3.14159*a*a*a*1e-30,n)*f;
  cout <<"Z0 = "<< Z0 << "kgA^2/s^2\n";*/


  double e=maths::inf(),eo=maths::inf();
  double pc=4,it=12;
  //std::vector<double> Z0;
  std::vector<Particle> CONFD;
  for(EX=0.3;EX<=3;EX+=0.3)
  {
    e=maths::inf();
    while(e>0)
    {
      conf=randomConfig(n);
      e=energy(conf);
    }
    accM=0;
    std::vector<double> energies;
    for(int i=0;i<nmoves;i++)
    {
      energies.push_back(e);
      e=mcmove(conf,e);
    }
    //f=((double)energies.size())/nmoves;
    cout <<(int)(EX/0.29999)<<"\t"<< stats::mean(energies) <<"\t"<<accM<<"\n"; //<<f<<"\n";
  }

}
int mainold2()
{
  std::vector<Particle> conf;
  double e,s=0,nf=0;
  EX=1;
  for(n=1;n<7;n++)
  {
    std::vector<double> energies;
    uwp=0;
    uljp=0;
    for(int i=0;i<nmoves;i++)
    {
      conf=randomConfigInclusive(n);
      energies.push_back(exp(-LJPot(conf)/(k*T)));
      //cout << energies[energies.size()-1]<<"\n";
    }
    e=stats::mean(energies);
    s=pow((4.0/3.0)*3.141592*a*a*a,n)*e;
    cout << n << "\t"<<s << "\n";
    /*s=0; nf=0;
    for(double& ev : energies)
    {
      s+=ev*exp(-ev/(k*T));
      nf+=exp(-ev/(k*T));
    }
    cout <<n<<"\t"<< (double)(s/nf)<<"\t"<<uljp<<" "<<uwp<<"\n";*/
  }
}
int mainold()
{

  int acc=0;
  //int n=4;
  //std::vector<double> energies;
  std::vector<Particle> conf;
  for(int n=1;n<2;n++)
  {
    acc=0;
    for(int i=0;i<nmoves;i++)
    {
      conf=randomConfig(n);
      if(HSP(conf)) acc++;
    }
    cout <<n<<"\t"<< acc <<"/"<<nmoves<<"\n";
  }
  //n=5;
  EX=1;
  double e=1e9;
  /*for(n=5;n<9;n++)
  {
    std::vector<double> energies;
    //conf=randomConfig(n);
    //e=LJPot(conf);
    for(int i=0;i<nmoves;i++)
    {
      conf=randomConfig(n);
      e=LJPot(conf);
      energies.push_back(e);
    }
    double s=0,nf=0;
    for(double& en : energies)
    {
      //cout << en << "\n";
      s+=en*exp(-en/(k*T));
      nf+=exp(-en/(k*T));
    }
    cout << n << "\n";
    cout <<"\t"<< s<<"\t"<<nf <<"\n";
    cout <<"\t"<< ((s/nf)*1e-3)/k << "\n";
  }*/
  for(n=5;n<9;n++)
  {
    std::vector<double> energies;
    e=1e9;
    while(e>0)
    {
      conf=randomConfig(n);
      e=LJPot(conf);
    }
    for(int i=0;i<nmoves;i++)
    {
      e=mcmove(conf,e);
      energies.push_back(e);
    }
    cout <<n<<"\t"<< stats::mean(energies) << "\n";
  }
}
