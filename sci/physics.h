#ifndef INCLUDED_PHYSICS_GENERAL
#define INCLUDED_PHYSICS_GENERAL 1
#include "maths/Eigen.h"
#include "commons/commons.h"
namespace physics
{
  class Particle
  {
    Eigen::Vector3d pos,vel,acc;
    double mass=1,rad=0;
    Eigen::Vector3d passiveacc;
  public:
    void* data=nullptr;
    Particle(const Eigen::Vector3d& r,const Eigen::Vector3d& v=Eigen::Vector3d(0,0,0))
    {
      pos=r;
      vel=v;
    }

    inline void setRadius(double r) {rad=r;}
    inline void setMass(double m) {mass=m;}
    inline double getMass() const {return mass;}
    inline double getRadius() const {return rad;}

    inline virtual void move(double time) //Simple Newtonian integration
    {
      Eigen::Vector3d facc=acc+passiveacc;
      pos+=vel*time+facc*((time*time)/2);
      vel+=facc*time;
    }
    inline void addForce(const Eigen::Vector3d& F) {acc+=F/mass;}
    inline void addPassiveForce(const Eigen::Vector3d& F) {passiveacc+=F/mass;}
    inline void addAcceleration(const Eigen::Vector3d& a) {acc+=a;}
    inline void addPassiveAcceleration(const Eigen::Vector3d& a) {passiveacc+=a;}
    inline void setVelocity(const Eigen::Vector3d& v) {vel=v;}
    inline void addVelocity(const Eigen::Vector3d& v) {vel+=v;}
    inline void setPosition(const Eigen::Vector3d& r) {pos=r;}
    inline const Eigen::Vector3d& getPosition() const {return pos;}
    inline const Eigen::Vector3d& getVelocity() const {return vel;}
    inline Eigen::Vector3d getNetForce() const {return acc*mass;}
    inline void clearAcceleration() {acc=Eigen::Vector3d(0,0,0);}
  };

  class Simulator
  {
    std::vector<Particle> particles;
    std::vector<Eigen::Vector3d> uniformacc;
    bool pbc=true;
    Eigen::Vector3d boxdims;
  public:
    Simulator() {}
    Simulator(const std::vector<Particle>& parts,const std::vector<Eigen::Vector3d>& acc,const Eigen::Vector3d& box)
    {
      uniformacc=acc;
      particles=parts;
      boxdims=box;
      for(Particle& part : particles) {for(const Eigen::Vector3d& v : acc) part.addPassiveAcceleration(v);}
    }

    inline void setPBC(bool b=true) {pbc=b;}
    inline bool isPBC() {return pbc;}
  };
}
#endif
