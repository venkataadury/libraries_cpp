#ifndef INCLUDED_SIM
#define INCLUDED_SIM 1
#include "commons/commons.h"
#include "maths/maths.h"
#include "maths/stats.h"

namespace simtools
{
  //HashMap<std::vector<int>,std::vector<double>> LJInteractions;
  class NoSuchParticleException : public exception
  {
  public:
    NoSuchParticleException() {}
  };

  /*static void addLJInteraction(const std::vector<int>& subt,const std::vector<double>& se) {LJInteractions.append(subt,se);}
  static void addLJInteraction(int p1,int p2,double e,double s)
  {
    std::vector<int> pt; pt.push_back(p1); pt.push_back(p2);
    std::vector<double> inv; pt.push_back(e); pt.push_back(s);
    addLJInteraction(pt,inv);
  }*/
  class Particle
  {
  public:
    double mass=1;
    Vector pos;
    Vector vel;
    void* data;

    Particle(double m=1) {mass=m;}
    Particle(double m,const Vector& r) : Particle(m) {pos=r;}
    Particle(double m,const Vector& r,const Vector& v) : Particle(m,r) {vel=v;}
    Particle(const Vector& r) : Particle(1,r) {}

    inline const double& getX() const {return pos[0][0];}
    inline double& getX() {return pos[0][0];}
    inline const double& getY() const {return pos[1][0];}
    inline double& getY() {return pos[1][0];}
    inline const double& getZ() const {return pos[2][0];}
    inline double& getZ() {return pos[2][0];}
    inline const double& getCoordinate(int i) const {return pos[i][0];}
    inline double& getCoordinate(int i) {return pos[i][0];}
    inline const Vector& getPosition() const {return pos;}
    inline Vector& getPosition() {return pos;}

    double getMass() {return mass;}
    const Vector& getVelocity() const {return vel;}
    Vector& getVelocity() {return vel;}
    void setVelocty(const Vector& v) {vel=v;}

    operator Vector() {return pos;}

    double move(double t=1) {pos=pos+(t*vel);}
    void setPosition(const Vector& v) {pos=v;}
  };
  /*class LJParticle : public Particle
  {
    int ptype=-1;
  public:
    LJParticle(int pt) {ptype=pt;}

    double LJPotential(const LJParticle& p)
    {
      std::vector<int> pair;
    }
  };*/

  class System
  {
  protected:
    std::vector<Particle> particles;
    double energy=0;
  public:
    System() {}
    System(std::vector<Particle> v) {particles=v;}

    const std::vector<Particle>& getParticles() const {return particles;}
    std::vector<Particle>& getParticles() {return particles;}
    const Particle& getRandomParticle() const {return randgen::choice(particles);}
    Particle& getRandomParticle() {return randgen::choice(particles);}
    Particle& getParticleAt(int i) {return particles[i];}
    const Particle& getParticleAt(int i) const {return particles[i];}

    //Particle& addParticle(Particle& p) {particles.push_back(p); return particles[particles.size()-1];}
    Particle& addParticle(Particle p) {particles.push_back(p); return particles[particles.size()-1];}
    Particle& removeParticle(int ind)  {return pop(particles,ind);}
    Particle& removeParticle(const Particle& p)
    {
      for(int i=0;i<particles.size();i++)
      {
        if(&particles[i]==&p)
          return pop(particles,i);
      }
      throw NoSuchParticleException();
    }
    Particle& addRandomParticle() {return addParticle(Particle(randVect()));}
    Particle& removeRandomParticle() {return removeParticle(randint(particles.size()));}

    inline int getParticleCount() const {return particles.size();}
    double getEnergy() const {return energy;}
    void setEnergy(double E) {energy=E;}


    inline virtual Vector getVectorFrom(const Vector& p1,const Vector& p2) const {return Vector(p2-p1);}
    inline virtual double distance(const Vector& p1,const Vector& p2) {return getVectorFrom(p1,p2).getMagnitude();}
    inline double distance(const Particle& p1,const Particle& p2) const {return distance(p1.getPosition(),p2.getPosition());}
    inline virtual Vector randVect() const;
  };
  class PBCCuboidBox : public System
  {
    double w=0,b=0,h=0;
  public:
    PBCCuboidBox() {}
    PBCCuboidBox(const Vector& v,const std::vector<Particle>& p=std::vector<Particle>()) : System(p) {w=v[0][0]; b=v[1][0]; h=v[2][0];}
    PBCCuboidBox(double r,const std::vector<Particle>& p=std::vector<Particle>()) : System(p) {w=b=h=r;}

    inline Vector getVectorFrom(const Vector& p1,const Vector& p2) const override
    {
      Vector dr=p2-p1;
      dr[0][0]=min(maths::abs(dr.getX()),w-maths::abs(dr.getX()));
      dr[1][0]=min(maths::abs(dr.getY()),b-maths::abs(dr.getY()));
      dr[2][0]=min(maths::abs(dr.getZ()),h-maths::abs(dr.getZ()));
      //cout << dr <<"\n";
      return dr;
    }

    inline Vector randVect() const override {return Vector(maths::random()*w,maths::random()*b,maths::random()*h);}
  };
  typedef PBCCuboidBox PBCBox;

  Vector randVect(double d=1)
  {
    Vector ret((maths::toss()?1:-1)*rand(),(maths::toss()?1:-1)*rand(),(maths::toss()?1:-1)*rand());
    ret=ret/ret.getMagnitude();
    return d*ret;
  }
  inline Vector System::randVect() const {return simtools::randVect();}

  class MonteCarlo
  {
  public:
    static double BETA;
    static double bruteForce(std::vector<Particle> (*genF)(),double (*ener)(const std::vector<Particle>& sys),double B=MonteCarlo::BETA,int nmoves=1e5)
    {
      double ef=0,nf=0;
      std::vector<Particle> conf;
      double en;
      for(int i=0;i<nmoves;i++)
      {
        conf=genF();
        en=ener(conf);
        if(isinf(en))
        {
          i--;
          continue;
        }
        //cout << "\t"<<en<<" "<<exp(-B*en)<<" "<<en*(exp(-B*en)) << "\n";
        ef+=(en*exp(-B*en));
        //cout <<en<<" "<<B<<"\t"<<(-B*en)<<" "<< (exp(-B*en)) << "\n";
        nf+=exp(-B*en);
      }
      return (ef/nf);
    }
    static double importance(std::vector<Particle> (*genF)(),double (*energy)(const std::vector<Particle>& sys),double B=MonteCarlo::BETA,double DX=0.1,int nmoves=1e5)
    {
      std::vector<double> energies;
      std::vector<Particle> conf;
      double e=maths::inf(),en;
      Vector mv;
      while(e>0)
      {
        conf=genF();
        e=energy(conf);
      }
      for(int i=0;i<nmoves;i++)
      {
        energies.push_back(e);
        Particle& p =randgen::choice(conf);
        mv=DX*maths::random()*randVect();
        p.pos=p.pos+mv;
        en=energy(conf);
        if(maths::random()<=exp((e-en)*B))
          e=en;
        else
          p.pos=p.pos-mv;
        /*if(i%100==0)
          cout <<i<<"\t"<< e << "\n";*/
      }
      return stats::mean(energies);
    }
  } MC;
}
double simtools::MonteCarlo::BETA=1/(1.38e-3*300);
namespace special
{
  class Layer : public simtools::System
  {
    int limit=100;
  public:
    Layer(int lim=100) : simtools::System() {limit=lim;}
    Layer(std::vector<simtools::Particle> v) : simtools::System(v) {}

    inline bool isFull() {return getParticleCount()>=limit;}

  };
}
#endif
