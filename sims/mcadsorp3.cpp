#include "sci/sim.h"
#include "maths/stats.h"
using maths::Vector;
using namespace simtools;

//rCH4,rZeo,LAMS,VS,V,Kb,SIGMA
//const double rCH4=1.994,rZeo=3.0;
double EX=1;
const double k=1.38e-3,T=300,DX=0.01;
const double LAMS=1.38e-3*1.93,VS=33.49,V=243.6,a=pow(sqrt(2)*V,1/3.0);
const double ETAS=1.38e-3*148.2,SIGMA=3.817,CUTOFF=SIGMA/5;//SIGMA/20;
inline static double l(double x) {return (1+12*x+25.2*x*x+12*x*x*x+pow(x,4))*pow(1-x,-10.0);}
inline static double m(double x) {return (1+x)*pow(1-x,-4.0);}
inline static double wallPot(double x)
{
  //cout << pow(VS/V,4)*l((x*x)/(a*a)) << "," << 2*(pow(VS/V,2))*m((x*x)/(a*a)) << "\t";
  return ((x<=CUTOFF)?1e9:0)+LAMS*(pow(VS/V,4)*l((x*x)/(a*a))-2*(pow(VS/V,2))*m((x*x)/(a*a)));
}
inline static double LJPot(double d) //Between CH4s
{/*cout << EX*(4*ETAS*(pow(SIGMA/d,12)-pow(SIGMA/d,6))) << ",";*/ return ((d<=CUTOFF)?maths::inf():0)+EX*(4*ETAS*(pow(SIGMA/d,12)-pow(SIGMA/d,6)));}

static double energy(const std::vector<Particle>& parts)
{
  double E=0;
  for(int i=0;i<parts.size();i++)
  {
    for(int j=i+1;j<parts.size();j++)
      E+=LJPot(Vector(parts[i].getPosition()-parts[j].getPosition()).getMagnitude());
    //cout << parts[i].getPosition().getMagnitude() << "\n";
    E+=wallPot(parts[i].getPosition().getMagnitude());
  }
  return E;
}
std::vector<Particle> randomConfig(int n)
{
  std::vector<Particle> ret;
  for(int i=0;i<n;i++)
    ret.push_back(Particle(a*maths::random()*randVect()));
  return ret;
}
double mcmove(std::vector<Particle>& conf,double& en)
{
  Vector mv=DX*maths::random()*randVect();
  Particle& p=randgen::choice(conf);
  p.pos=p.pos+mv;
  double e=energy(conf);
  //cout << (en-e)/(k*T) << "\n";
  if(maths::random()<=exp((en-e)/(k*T)))
    return e;
  else
  {
    p.pos=p.pos-mv;
    return en;
  }
}
int main()
{
  int nmoves=1e5;
  int n=5;
  double Z0;
  std::vector<Particle> conf;
  EX=0;
  int acc=0;
  double f;
  for(int i=0;i<nmoves;i++)
  {
    conf=randomConfig(n);
    if(energy(conf)>0)
      acc++;
  }
  f=(double)acc/nmoves;
  Z0=pow((4/3.0)*3.14159*a*a*a,n)*f;
  cout <<"Z0 = "<< Z0 << "\n";


  double e=maths::inf();
  double pc=4,it=12;
  //std::vector<double> Z0;

  for(int I=0;I<=it;I++)
  {
    EX=I/pc;
    //cout << EX << "\n";
    e=maths::inf();
    while(e>1)
    {
      conf=randomConfig(n);
      e=energy(conf);
    }
    std::vector<double> energies;
    for(int i=0;i<nmoves;i++)
    {
      energies.push_back(e);
      e=mcmove(conf,e);
      //cout << e << "\n";
      //cout <<I<<"; "<< e << "\n";
    }
    //f=((double)energies.size())/nmoves;
    cout <<I<<"\t"<< stats::mean(energies)*Z0 <<"\n"; //<<f<<"\n";
  }

  /*std::vector<Particle> parts;
  parts.push_back(Particle(Vector(-15,0,0)));
  parts.push_back(Particle(Vector(15,0,0)));*/
  /*double D=30;
  double dX=2e-2;
  double et;
  bool f1=false;
  cout << a<<"\n";
  for(int i=0;i<nmoves;i++)
  {
    et=e;
    e=LJPot(D);
    if(!f1 && e>et)
    {
      //cout << D <<" "<< e <<"\n";
      f1=true;
    }
    if(e>0)
      break;
    D-=dX;
  }
  cout << D <<" "<<e<< "\n";*/
}
