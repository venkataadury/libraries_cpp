#ifndef INCLUDED_EVO
#define INCLUDED_EVO 1
#include "maths/stats.h"
#include <map>
#include "Eigen337/Eigen/Dense"
#include "Eigen337/Eigen/Core"
#ifndef DISABLE_RENDERER
#include "draw/graphics.h"
#endif

class NotImplementedException : public std::exception {};
template<class T,class U> std::ostream& operator<<(std::ostream& os,const std::map<T,U>& mp)
{
  os << "{";
  for(const auto& k : mp) os << get<0>(k) <<":"<< get<1>(k) <<", ";
  os << "}";
  return os;
}
template<class T> static void sortedInsert(std::vector<T>& v, T val)
{
  int i;
  for(i=0;i<v.size();i++)
  {
    if(v[i]>=val) break;
  }
  if(i==v.size()) v.push_back(val);
  else v.insert(v.begin()+i,val);
}

namespace evolution
{
  template<class T> std::ostream& operator<<(std::ostream& os,const std::vector<T>& v)
  {
    os << "{";
    for(const T& el : v) os << el << ", ";
    os <<"}";
  }
  //Excpetions
  class Organism;
  class Environment;
  class InvalidTraitValueException : public exception {};
  class EndofFramesException : public exception {};
  class NoSuchTraitException : public exception
  {
  public:
    NoSuchTraitException() {}
    NoSuchTraitException(const std::string& tn) {std::cerr << "'"<<tn<<"' - Trait not found\n";}
  };

  static Eigen::Vector3d getRandomUnitVector(bool use2d=false)
  {
    Eigen::Vector3d randdir(randgen::throwarandompoint(-1,1),randgen::throwarandompoint(-1,1),(use2d)?0:randgen::throwarandompoint(-1,1));
    randdir/=randdir.norm();
    return randdir;
  }

  class TraitSpecifics
  {
  protected:
    double min,max;
    std::string name;
    double mutationStd=0; // By default -Stable inheritence
  public:
    TraitSpecifics(std::string traitname="NullTrait",double minv=-1.0/0.0,double maxv=1.0/0.0)
    {
      min=minv;
      max=maxv;
      name=traitname;
    }

  public:
    inline double getMaxValue() const {return max;}
    inline double getMinValue() const {return min;}
    inline void setMutationRate(double v) {mutationStd=v;}
    inline double getMutationRate() const {return mutationStd;}
    inline const std::string& getName() const {return name;}

    double generateVariation(double forcemean=0.0) const {return randgen::throwarandompoint_normal(forcemean,mutationStd);}
    //virtual double generateVariation(double forcemean=0.0) const {return 0.0;}
  };

  static std::map<std::string,TraitSpecifics> traits;

  class Trait
  {
  public:
    const TraitSpecifics& specs;
    double myval;
  public:
    Trait(const TraitSpecifics& tsp,double val,bool cast=false) : specs(tsp)
    {
      if(!cast && (val<tsp.getMinValue() || val>tsp.getMaxValue())) throw InvalidTraitValueException();
      myval=val;
      if(myval>specs.getMaxValue()) myval=specs.getMaxValue();
      else if(myval<specs.getMinValue()) myval=specs.getMinValue();
    }
    Trait(const Trait& t) : specs(t.specs) {myval=t.myval;}

    inline Trait replicate() const {return Trait(specs,myval+specs.generateVariation(),true);}
    Trait& operator=(const Trait& t) {myval=t.myval; return *this;}
    inline const std::string& getName() const {return specs.getName();}
  };
  class StateVariable
  {
  public:
    std::string name;
    double val;
    bool heritableTrait=false;
    StateVariable(const std::string& n="None",double v=0,bool herit=false)
    {
      name=n;
      val=v;
      heritableTrait=herit;
    }

    const std::string& getName() const {return name;}
    double getData() const {return val;}
    void setData(double d) {val=d;}
  };

  static bool addTraitSpecifics(std::string name,double minv=-1.0/0.0,double maxv=1.0/0.0,double mutrate=0.0)
  {
    TraitSpecifics newtrait(name,minv,maxv);
    newtrait.setMutationRate(mutrate);
    traits.emplace(std::make_pair(name,newtrait));
    return true;
  }
  inline static bool addTraitSpecifics(const TraitSpecifics& tsp) {traits.emplace(tsp.getName(),tsp);}
  inline static bool addConstantTraitSpecifics(std::string name,double val,double err=1e-3)
  {
    TraitSpecifics newtrait(name,val-err,val+err);
    newtrait.setMutationRate(0.0d);
    traits.emplace(std::make_pair(name,newtrait));
    return true;
  }
  inline static TraitSpecifics& getTraitSpecifics(const std::string& name) {return traits[name];}

  /*class GenericTraitGenerator
  {
  protected:

  public:
    GenericTraitGenerator(const std::string& name) : GenericTraitGenerator(getTraitSpecifics(name)) {}
    GenericTraitGenerator(const TraitSpecifics& tsp) : specs(tsp) {}
    GenericTraitGenerator(const GenericTraitGenerator& tgen) : specs(tgen.specs) {}

    virtual Trait generate() const {throw NotImplementedException();}

  };*/
  class TraitGenerator
  {
    const TraitSpecifics& specs;
    double mean,std;
  public:
    TraitGenerator(const std::string& name,double mv,double sv) : TraitGenerator(getTraitSpecifics(name),mv,sv) {}
    TraitGenerator(const TraitSpecifics& tsp,double mv,double sv) : specs(tsp) {mean=mv; std=sv;}
    TraitGenerator(const TraitGenerator& tgen) : specs(tgen.specs) {mean=tgen.mean; std=tgen.std;}

    Trait generate() const {return Trait(specs,randgen::throwarandompoint_normal(mean,std),true);}
    inline double getMean() const {return mean;}
    inline double getStd() const {return std;}
    inline const std::string& getName() const {return specs.getName();}
    virtual std::string toString() const {return getName()+"("+to_string(getMean())+","+to_string(getStd())+")";}
  };
  typedef std::vector<TraitGenerator> TraitGroup;

  static std::vector<TraitGenerator> createTraitGenerators(const std::vector<std::string>& names, const std::vector<double>& means,const std::vector<double>& stds)
  {
    std::vector<TraitGenerator> tgens;
    int K=0;
    for(const std::string& n : names)
    {
      tgens.push_back(TraitGenerator(getTraitSpecifics(n),means[K],stds[K]));
      K++;
    }
    return tgens;
  }

  class OrganismTemplate
  {
    std::vector<TraitGenerator> tgens;
  public:
    std::string ID="None";
    std::string speedtrait="None";
  public:
    OrganismTemplate() {}
    OrganismTemplate(const std::string& IDv,const std::vector<std::string>& traitnames,const std::vector<double>& means,const std::vector<double>& stds) : OrganismTemplate(IDv,createTraitGenerators(traitnames,means,stds)) {}
    OrganismTemplate(const std::string& IDv,const std::vector<TraitGenerator>& tgensin)
    {
      ID=IDv;
      for(const auto& gen : tgensin) tgens.push_back(gen);
    }

    const std::vector<TraitGenerator>& getTraitGenerators() const {return tgens;}
    void setSpeedTrait(const std::string& st) {speedtrait=st;}
    const std::string& getSpeedTrait() const {return speedtrait;}
    std::vector<Trait> drawTraits() const
    {
      std::vector<Trait> ret;
      for(const TraitGenerator& tgen : tgens) ret.push_back(tgen.generate());
      return ret;
    }
  };

  std::map<std::string,OrganismTemplate> organisms;
  inline static bool addOrganismTemplate(std::string name,const std::vector<std::string>& traitnames,const std::vector<double>& means,const std::vector<double>& stds)
  {
    organisms.emplace(make_pair(name,OrganismTemplate(name,traitnames,means,stds)));
    return true;
  }
  inline static bool addOrganismTemplate(std::string name,const std::vector<TraitGenerator>& tgensin)
  {
    organisms.emplace(make_pair(name,OrganismTemplate(name,tgensin)));
    return true;
  }
  inline static bool addOrganismTemplate(const OrganismTemplate& otemp) {organisms.emplace(make_pair(otemp.ID,otemp)); return true;}
  inline static const OrganismTemplate& getOrganismTemplate(const std::string& nm) {return organisms[nm];}
  inline static OrganismTemplate& getEditableOrganismTemplate(const std::string& nm) {return organisms[nm];}

  class Organism
  {
  protected:
    const OrganismTemplate& temp;
    std::vector<Trait> traits;
    std::map<std::string,StateVariable> vars;
  public:
    Organism(const OrganismTemplate& t) : temp(t) {traits=temp.drawTraits();}
    Organism(const OrganismTemplate& t,const std::vector<Trait>& mytraits) : temp(t) {traits=mytraits;}
    /*~Organism()
    {
      for(auto& mv : vars)
      {
        delete get<1>(mv).getDataAs<>();
      }
    }*/

    inline bool isType(const std::string& tp) const {return temp.ID==tp;}
    inline void addVariable(const std::string& svn,double dat,bool herit=false) {vars.emplace(make_pair(svn,StateVariable(svn,dat,herit)));}
    inline bool hasVariable(const std::string& n) const {return vars.find(n)!=vars.end();}
    inline void setVariable(const std::string& svn,double dat) {vars[svn].setData(dat);}
    inline double getVariable(const std::string& n) const {return vars.at(n).getData();}

    inline const std::string& getName() const {return temp.ID;}
    inline const std::vector<Trait>& getTraits() const {return traits;}
    inline bool hasTrait(const std::string& tn) const
    {
      for(const Trait& t : traits)
      {
        if(t.getName()==tn) return true;
      }
      return false;
    }
    inline double getTrait(const std::string& tn) const {return getTraitObj(tn).myval;}
    inline const Trait& getTraitObj(const std::string& tn) const
    {
      for(const Trait& t : traits)
      {
        if(t.getName()==tn) return t;
      }
      throw NoSuchTraitException();
    }
    inline Trait& getTraitObj(const std::string& tn)
    {
      for(Trait& t : traits)
      {
        if(t.getName()==tn) return t;
      }
      throw NoSuchTraitException();
    }

    inline std::vector<Trait> getProgenyTraits() const
    {
      std::vector<Trait> fintraits;
      for(const Trait& t : traits) fintraits.push_back(t.replicate());
      return fintraits;
    }

    virtual Organism* replicate() const {return new Organism(temp,getProgenyTraits());}
  };

  class Interaction
  {
  public:
    Interaction() {}

    virtual void resolve(Environment& env,double t=1.) const =0;
  };

  class PostitionedOrganism : public Organism
  {
    Eigen::Vector3d position,velocity;
  public:
    PostitionedOrganism(const OrganismTemplate& t,double x=0,double y=0,double z=0) : Organism(t) {position=Eigen::Vector3d(x,y,z); velocity=getRandomUnitVector(::abs(z)<1e-3);}
    PostitionedOrganism(const OrganismTemplate& t,const std::vector<Trait>& mytraits,double x=0,double y=0,double z=0) : Organism(t,mytraits) {position=Eigen::Vector3d(x,y,z); velocity=getRandomUnitVector(::abs(z)<1e-3);}
    PostitionedOrganism(const OrganismTemplate& t,const Eigen::Vector3d& p) : Organism(t) {position=p; velocity=getRandomUnitVector(::abs(p(2))<1e-3);}
    PostitionedOrganism(const OrganismTemplate& t,const std::vector<Trait>& mytraits,const Eigen::Vector3d& p) : Organism(t,mytraits) {position=p; velocity=getRandomUnitVector(::abs(p(2))<1e-3);}

    inline double getX() const {return position(0);}
    inline double getY() const {return position(1);}
    inline double getZ() const {return position(2);}
    inline void setX(double v) {position(0)=v;}
    inline void setY(double v) {position(1)=v;}
    inline void setZ(double v) {position(2)=v;}
    inline double getXVel() const {return velocity(0);}
    inline double getYVel() const {return velocity(1);}
    inline double getZVel() const {return velocity(2);}
    const std::string& getSpeedTrait() const {return temp.getSpeedTrait();}

    inline const Eigen::Vector3d& getPosition() const {return position;}
    inline Eigen::Vector3d& getPosition() {return position;}
    inline const Eigen::Vector3d& getVelocity() const {return velocity;}
    inline Eigen::Vector3d& getVelocity() {return velocity;}
    inline void setVelocity(const Eigen::Vector3d& v)
    {
      velocity=v;
      if(this->hasTrait(temp.getSpeedTrait()))
      {
        double st=this->getTrait(temp.getSpeedTrait());
        if(velocity.norm() > st) velocity*=st/velocity.norm();
      }
    }
    inline Eigen::Vector3d& updateVelocity(const Eigen::Vector3d& delta) {velocity+=delta; return velocity;}

    double distaceTo(PostitionedOrganism* p2) const {return (this->position-p2->position).norm();}

    inline void updatePosition(double time=1) {position+=velocity*time;}
    virtual inline Eigen::Vector2d get2DProjection() const {return Eigen::Vector2d(position(0),position(1));}
    virtual PostitionedOrganism* replicate() const override
    {
      Eigen::Vector3d randdir=getRandomUnitVector(::abs(position(2))<1e-3);
      PostitionedOrganism* ret=new PostitionedOrganism(PostitionedOrganism(temp,getProgenyTraits(),position));
      if(ret->hasTrait(temp.getSpeedTrait())) ret->velocity=randdir*ret->getTrait(temp.getSpeedTrait());
      else ret->velocity=randdir*velocity.norm();
      return ret;
    }
  };

  class Environment
  {
  protected:
    std::vector<PostitionedOrganism*> organisms;
    std::vector<Interaction*> interactions;
    double timestep=1;
  public:
    Environment(double ts=1,const std::vector<PostitionedOrganism*> orgs=std::vector<PostitionedOrganism*>()) {timestep=ts;organisms=orgs;}

    inline int getOrganismCount() const {return organisms.size();}
    inline void killOrganism(int ind)
    {
      delete organisms[ind];
      organisms.erase(organisms.begin()+ind);
    }
    inline void addOrganism(PostitionedOrganism* org) {organisms.push_back(org);}
    inline PostitionedOrganism* getOrganism(int no) {return organisms[no];}
    inline void addInteraction(Interaction* intr) {interactions.push_back(intr);}
    inline const std::vector<Interaction*>& getInteractions() const {return interactions;}
    inline double getTimestep() const {return timestep;}
    inline const std::vector<PostitionedOrganism*>& getAllOrganisms() const {return organisms;}
    virtual void fixPositions() {}
    void updateAllPositions(double t=-1)
    {
      if(t<0) t=timestep;
      for(PostitionedOrganism* porg : organisms) porg->updatePosition(t);
      fixPositions();
    }
    inline void solveAllInteractions(double t=-1)
    {
      for(Interaction* intr : interactions) intr->resolve(*this,t);
    }
    virtual void simulate(double t=-1) {updateAllPositions(t); solveAllInteractions(t);}
    virtual void runSimulation(const std::string& outf,double t0=0,double tf=100,double dt=1.,int writeInterval=1)
    {
      ofstream of; of.open(outf);
      runSimulation(of,t0,tf,dt,writeInterval);
      of.close();
    }
    virtual void runSimulation(ostream& outp,double t0=0,double tf=100,double dt=1.,int writeInterval=1)
    {
      double st=t0;
      int K=0;
      while(st<tf)
      {
        simulate(dt);
        st+=dt;
        if(K%writeInterval==0) dumpToXYZ(outp,true);
        K++;
        cout << K << "\n";
      }
    }

    virtual void addPopulation(const OrganismTemplate& otemp,int N,double xmin=0,double ymin=0,double zmin=0,double xmax=1,double ymax=1,double zmax=1)
    {
      for(int i=0;i<N;i++)
      {
        addOrganism(new PostitionedOrganism(otemp,randgen::throwarandompoint(xmin,xmax),randgen::throwarandompoint(ymin,ymax),randgen::throwarandompoint(zmin,zmax)));
        if(organisms[organisms.size()-1]->hasTrait(otemp.getSpeedTrait())) organisms[organisms.size()-1]->setVelocity(organisms[organisms.size()-1]->getVelocity()*organisms[organisms.size()-1]->getTrait(otemp.getSpeedTrait()));
      }
      fixPositions();
    }

    virtual void dumpToXYZ(std::ostream& of,bool writevel=true) const
    {
      of << "Environment\n";
      of<<" " << this->getOrganismCount()<<"\n";
      for(PostitionedOrganism* po : this->getAllOrganisms())
      {
        of << po->getName()<<" "<<po->getX()<<" "<<po->getY()<<" "<<po->getZ();
        if(writevel) of << " "<<po->getXVel()<<" "<<po->getYVel()<<" "<<po->getZVel();
        of << "\n";
      }
    }
    virtual void dumpToXYZ(const std::string& fn,bool writevel=true) const
    {
      ofstream of; of.open(fn);
      this->dumpToXYZ(of,writevel);
      of.close();
    }

    std::vector<double> getTraitDistribution(const std::string& tname,const std::string& type="All") const
    {
      std::vector<double> ret;
      for(PostitionedOrganism* porg : organisms)
      {
        if(type!="All" && type!=porg->getName()) continue;
        if(porg->hasTrait(tname)) ret.push_back(porg->getTrait(tname));
      }
      return ret;
    }
    int getOrganismCount(const std::string& oname) const
    {
      int N=0;
      for(PostitionedOrganism* porg : organisms)
      {
        if(porg->getName()==oname) N++;
      }
      return N;
    }
  };

  class Environment2DBox : public Environment
  {
    double wid,hei;
  public:
    Environment2DBox(double ts,double width,double height,const std::vector<PostitionedOrganism*> orgs=std::vector<PostitionedOrganism*>()) : Environment(ts,orgs) {wid=width,hei=height; this->fixPositions();}

    inline double getWidth() const {return wid;}
    inline double getHeight() const {return hei;}

    virtual void fixPositions() override
    {
      for(PostitionedOrganism* porg : organisms)
      {
        porg->setZ(0);
        /*porg->getPosition()(0)%=wid;
        porg->getPosition()(1)%=hei;*/
        while(porg->getX()>wid) porg->setX(porg->getX()-wid);
        while(porg->getX()<0) porg->setX(wid+porg->getX());
        while(porg->getY()>hei) porg->setY(porg->getY()-hei);
        while(porg->getY()<0) porg->setY(hei+porg->getY());
      }
    }

    virtual void dumpToXYZ(std::ostream& of,bool writevel=true) const override
    {
      of << "Environment 2D\n";
      of<<" " << this->getOrganismCount()<<" "<<this->getWidth()<<" "<<this->getHeight()<<"\n";
      for(PostitionedOrganism* po : this->getAllOrganisms())
      {
        of << po->getName()<<" "<<po->getX()<<" "<<po->getY()<<" "<<po->getZ();
        if(writevel) of << " "<<po->getXVel()<<" "<<po->getYVel()<<" "<<po->getZVel();
        of << "\n";
      }
    }
    virtual void dumpToXYZ(const std::string& fn,bool writevel=true) const override
    {
      ofstream of; of.open(fn);
      this->dumpToXYZ(of,writevel);
      of.close();
    }
    //virtual void dumpToXYZ(const std::string& fn,bool writevel=true) const override {}
  };


  //Output methods
  std::ostream& operator<<(std::ostream& os,const TraitSpecifics& tsp)
  {
    os << tsp.getName()<<" - ["<<tsp.getMinValue()<<","<<tsp.getMaxValue()<<"]";
    return os;
  }

  std::ostream& operator<<(std::ostream& os,const Trait& t)
  {
    os << t.specs.getName() << "="<<t.myval;
    return os;
  }

  inline std::ostream& operator<<(std::ostream& os,const TraitGenerator& tgen) {return (os << tgen.toString());}
  std::ostream& operator<<(std::ostream& os,const OrganismTemplate& otemp)
  {
    os <<"Organism(ID="<< otemp.ID<<") with traits [";
    for(const auto& tgen : otemp.getTraitGenerators()) os << tgen <<", ";
    os << "]";
    return os;
  }
  std::ostream& operator<<(std::ostream& os,const Organism& o)
  {
    os << o.getName()<<" - [";
    for(const Trait& t : o.getTraits()) os << t<<",";
    os << "]";
    return os;
  }
  std::ostream& operator<<(std::ostream& os,const PostitionedOrganism& o)
  {
    os << o.getName()<<"("<<o.getX()<<","<<o.getY()<<","<<o.getZ()<<") - [";
    for(const Trait& t : o.getTraits()) os << t<<",";
    os << "]";
    return os;
  }
  std::ostream& operator<<(std::ostream& os,const Environment& env)
  {
    os <<" "<<env.getOrganismCount()<<"\n";
    os << "Environment\n";
    for(PostitionedOrganism* po : env.getAllOrganisms()) os << po->getName()<<" "<<po->getX()<<" "<<po->getY()<<" "<<po->getZ()<<"\n";
    return os;
  }
  /*std::ostream& operator<<(std::ostream& os,const Environment2DBox& env)
  {
    os <<" "<<env.getOrganismCount()<<" "<<env.getWidth()<<" "<<env.getHeight()<<"\n";
    os << "Environment 2D\n";
    for(PostitionedOrganism* po : env.getAllOrganisms()) os << po->getName()<<" "<<po->getX()<<" "<<po->getY()<<" "<<po->getZ()<<"\n";
    return os;
  }*/

  static void saveState(std::ostream& os) //Only supported for TraitGenerator traits (not generic generators)
  {
    int ntraits=traits.size();
    os << ntraits << "\n";
    for(const auto& t : traits) os << get<0>(t)<<" "<<get<1>(t).getMinValue()<<" "<<get<1>(t).getMaxValue()<<" "<<get<1>(t).getMutationRate() << "\n";
    int norg=organisms.size();
    os << norg << "\n";
    for(const auto& o : organisms)
    {
      int tgs=get<1>(o).getTraitGenerators().size();
      os << get<0>(o)<<" "<<tgs << "\n";
      for(const TraitGenerator& txg : get<1>(o).getTraitGenerators())
      {
        const auto& tg = (const TraitGenerator&)txg;
        os << tg.getName()<<" "<<tg.getMean()<<" "<<tg.getStd()<<"\n";
      }
    }
  }
  static void saveState(const std::string& fn)
  {
    std::ofstream of; of.open(fn);
    saveState(of);
    of.close();
  }

  static void loadState(std::istream& in)
  {
    int ntrait;
    in >> ntrait;
    std::string tname,oname;
    double mv,mV,mR;
    std::string curline;
    for(int i=0;i<ntrait;i++)
    {
      //getline(in,curline);
      in >> tname >> mv >> mV >> mR;
      cout << tname << " "<<mv<<" "<<mV<<" "<<mR<<"\n";
      addTraitSpecifics(tname,mv,mV,mR);
    }
    int norg;
    in >> norg;
    int nf;
    for(int i=0;i<norg;i++)
    {
      in >> oname >> nf;
      cout << oname <<"\t";
      std::vector<std::string> tnames;
      std::vector<double> means,stds;
      double mn,st;
      for(int j=0;j<nf;j++)
      {
        in >> tname >> mn >> st;
        cout << tname << " "<<mn<<" "<<st<<";";
        means.push_back(mn);
        stds.push_back(st);
        tnames.push_back(tname);
      }
      cout << "\n";
      addOrganismTemplate(oname,tnames,means,stds);
    }
  }
  static void loadState(const std::string& fn)
  {
    std::ifstream is; is.open(fn);
    loadState(is);
    is.close();
  }

  namespace interactions
  {
    static std::string PREY_SPEED="MaxSpeed",PRED_SPEED="MaxSpeed",PREY_VISION="Vision",PRED_VISION="Vision",PREY_WEIGHTAGE="R-Weight",PREY_ACC="MaxAcceleration"; //Can be altered
    class PreyAvoidance : public Interaction
    {
      std::string pred,prey;
    public:
      double defaultSpeed=0.,defaultWeight=0.;
    public:
      PreyAvoidance(const std::string& preyv="Prey",const std::string& predator="Predator") : Interaction() {prey=preyv; pred=predator;}

      void setDefaultSpeed(double ds) {defaultSpeed=ds;}
      void setDefaultWeight(double wt) {defaultWeight=wt;}
      virtual void resolve(Environment& env,double time) const override
      {
        if(time<0) time=1.;
        std::vector<PostitionedOrganism*> preds;
        std::vector<PostitionedOrganism*> preys;
        std::vector<std::vector<bool>> preyfeat;
        for(PostitionedOrganism* org : env.getAllOrganisms())
        {
          if(org->isType(pred)) preds.push_back(org);
          if(org->isType(prey)) // Predator and prey can be same???
          {
            preyfeat.push_back(std::vector<bool>());
            preyfeat[preyfeat.size()-1].push_back(org->hasTrait(PREY_SPEED));
            preyfeat[preyfeat.size()-1].push_back(org->hasTrait(PREY_VISION));
            preyfeat[preyfeat.size()-1].push_back(org->hasTrait(PREY_WEIGHTAGE));
            preyfeat[preyfeat.size()-1].push_back(org->hasTrait(PREY_ACC));
            preys.push_back(org);
          }
        }
        for(int i=0;i<preys.size();i++)
        {
          Eigen::Vector3d acc(0,0,0);
          double spd=(preyfeat[i][0])?preys[i]->getTrait(PREY_SPEED):defaultSpeed;
          double maxspd=(preyfeat[i][0])?preys[i]->getTraitObj(PREY_SPEED).specs.getMaxValue():spd;
          double accn=(preyfeat[i][3])?preys[i]->getTrait(PREY_ACC):maxspd/2;
          for(PostitionedOrganism* predn : preds)
          {
            Eigen::Vector3d rundir=preys[i]->getPosition()-predn->getPosition();
            double vsn=(preyfeat[i][1])?preys[i]->getTrait(PREY_VISION):spd;
            double wt=(preyfeat[i][2])?preys[i]->getTrait(PREY_WEIGHTAGE):defaultWeight;
            double r=rundir.norm();
            if(i>vsn) continue;
            rundir/=::pow(r,1+wt);
            acc+=rundir;
          }
          if(acc.norm()<1e-6) continue;
          acc*=accn/acc.norm();
          preys[i]->setVelocity(preys[i]->getVelocity()+acc);
        }
      }
    };


    class DefaultBirthAndDeath : public Interaction
    {
      std::string birth,death,ondeath="None";
      double K=0;
      bool ucc=false;
    public:
      DefaultBirthAndDeath(double carrying_capacity=1.0/0.0,bool universal_cc=false,const std::string& b="Birth",const std::string& d="Death",const std::string& dieto="None")
      {
        birth=b;
        death=d;
	ondeath=dieto;
        K=carrying_capacity;
        ucc=universal_cc;
      }

      virtual void resolve(Environment& env, double time) const override
      {
        if(time<0) time=1.;
        std::vector<int> dels;
        int maxs=env.getOrganismCount();
        std::map<std::string,int> popsizes;
        for(PostitionedOrganism* porg : env.getAllOrganisms())
        {
          if(popsizes.find(porg->getName())==popsizes.end()) popsizes.emplace(make_pair(porg->getName(),1));
          else popsizes[porg->getName()]++;
        }
        //operator<<(cout,(std::map<std::string,int>)popsizes);
        for(int i=0;i<maxs;i++)
        {
          PostitionedOrganism* corg=env.getOrganism(i);
          bool hb=corg->hasTrait(birth);
          if(hb)
          {
            if(randgen::throwarandompoint(0,1)<corg->getTrait(birth)*time*(1.-(double)popsizes[corg->getName()]/(2.0*K)))
            {
              PostitionedOrganism* norg=corg->replicate();
              env.addOrganism(norg);
            }
          }
          if((corg->hasTrait(death) && randgen::throwarandompoint(0,1)<corg->getTrait(death)*time) || (hb && randgen::throwarandompoint(0,1)<corg->getTrait(birth)*time*((double)popsizes[corg->getName()]/(2.0*K))))
          {
            dels.push_back(i);
            if(ondeath!="None")
            {
              PostitionedOrganism* norg=new PostitionedOrganism(getOrganismTemplate(ondeath),corg->getPosition());
              //cout << "Adding population of "<<ondeath<<"using speed-trait "<<corg->getSpeedTrait()<<"\n";
              //cout << norg->getVelocity().norm()<<" is the speed of new organism: "<<norg->getName()<<"\n";
              norg->setVelocity(Eigen::Vector3d(0,0,0));
              env.addOrganism(norg);
            }
          }
        }
        int SHIFT=0;
        for(int n : dels) {env.killOrganism(n-SHIFT++);}
      }
    };

    static std::string PRED_EAT="EatCount",PREY_EAT="Eaten",PRED_LIVES="Lives";
    static double LIVES=5;
    class PredatorEating : public Interaction //Includes birth?
    {
      std::string pred,prey;
      double rad=4,prob=1,boost=0;
      bool autodeath=false,autobirth=false;
      double birthmod=0,deathmod=0;
    public:
      PredatorEating(const std::string& prd,const std::string& pry,double r,double p=1,double bst=0) {pred=prd; prey=pry; rad=r;prob=p; boost=bst;}

      PredatorEating& allowBirth(double bm) {birthmod=bm;autobirth=true; return *this;}
      PredatorEating& allowDeath(double dm) {deathmod=dm; autodeath=true; return *this;}

      virtual void resolve(Environment& env,double time) const override
      {
        std::vector<int> eaters;
        std::vector<int> dels;
        for(int i=0;i<env.getOrganismCount();i++)
        {
          PostitionedOrganism* porg=env.getOrganism(i);
          if(porg->isType(pred))
          {
            if(!porg->hasVariable(PRED_EAT)) porg->addVariable(PRED_EAT,0.5+boost);
            if((autodeath || autobirth) && !porg->hasVariable(PRED_LIVES)) porg->addVariable(PRED_LIVES,LIVES);
            eaters.push_back(i);
          }
          if(porg->isType(prey)) {if(!porg->hasVariable(PREY_EAT)) porg->addVariable(PREY_EAT,0);}
        }
        for(int ix : eaters)
        {
          PostitionedOrganism* porg=env.getOrganism(ix);
          for(int i=0;i<env.getOrganismCount();i++)
          {
            PostitionedOrganism* porg2=env.getOrganism(i);
            if(!porg2->isType(prey)) continue;
            if(porg2->getVariable(PREY_EAT)<0.5 && porg->distaceTo(porg2)<=rad && (prob>=1 || randgen::throwarandompoint(0,1)<prob))
            {
              sortedInsert(dels,i);
              porg2->setVariable(PREY_EAT,1);
              porg->setVariable(PRED_EAT,porg->getVariable(PRED_EAT)+1.0);
              cout << "Eat Predator! "<<porg->getVariable(PRED_EAT)<<"\n";;
            }
          }
          if(autodeath) {/*cout << "\t"<<deathmod << " on "<< porg->getVariable(PRED_EAT) << "\n";*/ porg->setVariable(PRED_EAT,max(0.0,porg->getVariable(PRED_EAT)-deathmod));}
          cout << "Predator foods: "<< porg->getVariable(PRED_EAT) << "\n";
        }

        if(autodeath)
        {
          for(int ix : eaters)
          {
            PostitionedOrganism* porg=env.getOrganism(ix);
            cout << porg->getVariable(PRED_LIVES) << " lives left\n";
            if(porg->getVariable(PRED_EAT)<=0.5) porg->setVariable(PRED_LIVES,porg->getVariable(PRED_LIVES)-1);
            else porg->setVariable(PRED_LIVES,min(porg->getVariable(PRED_LIVES)+3,LIVES));
            if(porg->getVariable(PRED_LIVES)<=0) {cout << "Die Predator! " << porg->getVariable(PRED_LIVES)<<" lives left\n"; sortedInsert(dels,ix);}
          }
        }
        if(autobirth)
        {
          for(int ix : eaters)
          {
            PostitionedOrganism* porg=env.getOrganism(ix);
            if(porg->getVariable(PRED_EAT)>=1.5 && porg->getVariable(PRED_LIVES)>=2)
            {
              if(randgen::throwarandompoint(0,1)<birthmod*porg->getVariable(PRED_EAT))
              {
                PostitionedOrganism* norg=porg->replicate();
                env.addOrganism(norg);
                porg->setVariable(PRED_EAT,porg->getVariable(PRED_EAT)-1);
              }
            } //porg->setVariable(PRED_LIVES,porg->getVariable(PRED_LIVES)-1);
          }
        }
        //cout << dels << "\n";
        int SHIFT=0; for(int n : dels) {env.killOrganism(n-SHIFT++);}
      }
    };

    namespace generic
    {
      class ConversionInteraction : public Interaction
      {
      protected:
        std::string source,targ;
        double prob=1;
      protected:
        ConversionInteraction() {source="None"; targ="None";}
      public:
        ConversionInteraction(const std::string& src,const std::string& trg,double p=1) {source=src; targ=trg; prob=p;}

        //ConversionInteraction& allowBirth(double bm) {birthmod=bm;autobirth=true; return *this;}
        //PredatorEating& allowDeath(double dm) {deathmod=dm; autodeath=true; return *this;}

        virtual bool conditionSatisfied(const Environment& env,const PostitionedOrganism* porg) const {return false;}

        virtual void resolve(Environment& env,double time) const override
        {
          //if(!conditionSatisfied(env)) return;
          std::vector<int> converters;
          std::vector<int> dels;
          for(int i=0;i<env.getOrganismCount();i++)
          {
            PostitionedOrganism* porg=env.getOrganism(i);
            if(porg->isType(source) && conditionSatisfied(env,porg) && (prob==1 || randgen::throwarandompoint(0,1)<prob)) {converters.push_back(i); dels.push_back(i);}
          }
          for(int ix : converters)
          {
            PostitionedOrganism* porg=env.getOrganism(ix);
            PostitionedOrganism* norg=new PostitionedOrganism(getOrganismTemplate(targ),porg->getPosition());
            norg->setVelocity(porg->getVelocity());
            env.addOrganism(norg);
          }
          //cout << dels << "\n";
          int SHIFT=0; for(int n : dels) {env.killOrganism(n-SHIFT++);}
        }
      };
    }
  }

  namespace statistics
  {
    class GeneralPositionalStatistic
    {
    public:
      double mean,std,other;
      std::string name="Statistic";
    public:
      GeneralPositionalStatistic(const std::string& n="Statistic") {name=n;}

      virtual void compute(const std::vector<PostitionedOrganism*>& frame) =0;
      double getMean() const {return mean;}
      double getStd() const {return std;}
      double getOther() const {return other;}

    };
    typedef GeneralPositionalStatistic GPS;

    class TraitStatistic : public GeneralPositionalStatistic
    {
      std::string type,statname;
    public:
      TraitStatistic(const std::string tn,const std::string& t="All") : GeneralPositionalStatistic(tn+"_stat") {type=t; statname=tn;}

      void compute(const std::vector<PostitionedOrganism*>& porgs) override
      {
        std::vector<double> statdata;
        for(PostitionedOrganism* p : porgs)
        {
          if(!p->isType(type)) continue;
          if(p->hasTrait(statname)) statdata.push_back(p->getTrait(statname));
        }
        mean=stats::mean(statdata);
        std=stats::stddev(statdata);
        other=statdata.size();
        if(statdata.size()==0) cout << "WARN: Statistic couldn't be computed. Trait not found: "<<name<<"\n";
      }
    };
  }

  #ifndef DISABLE_RENDERER
  namespace renderer
  {
    static std::vector<sf::Color> COLORS={sf::Color::Green,sf::Color::Blue, sf::Color::Red, sf::Color::Yellow};
    using namespace graphics;
    class DrawingBoard
    {
      std::vector<std::vector<PostitionedOrganism*>> frames;
      std::vector<std::pair<double,double>> shapes;
    public:
      int cF=0;
    public:
      DrawingBoard() {}
      DrawingBoard(const Environment& env) {frames=std::vector<std::vector<PostitionedOrganism*>>(1,env.getAllOrganisms()); deduceShapes();}
      DrawingBoard(const std::string& fn)
      {
        std::ifstream is; is.open(fn);
        frames=std::vector<std::vector<PostitionedOrganism*>>();
        appendFramesFrom(is);
        cout << frames.size() << " frames loaded: DrawingBoard\n";
        is.close();
      }

      void appendFramesFrom(std::istream& is)
      {
        std::string dump;
        int N;
        double Nx;
        double w,h;
        while(true)
        {
          getline(is,dump); //Frame name
          getline(is,dump);
          dump=stringfx::trim(dump);
          stringstream nss(dump);
          nss >> N >> w >> h;
          cout << N << " in "<<frames.size()<<"\n";
          shapes.push_back(make_pair(w,h));
          std::vector<PostitionedOrganism*> orgs;
          std::string oname;
          for(int i=0;i<N;i++)
          {
            Eigen::Vector3d pos(0,0,0),vel(0,0,0);
            getline(is,dump);
            stringstream ss(dump);
            ss >> oname >> pos(0) >> pos(1)>>pos(2);
            if(ss.tellg()!=-1) ss >> vel(0) >> vel(1) >> vel(2);
            PostitionedOrganism* norg = new PostitionedOrganism(organisms[oname],pos);
            norg->setVelocity(vel);
            orgs.push_back(norg);
          }
          this->frames.push_back(orgs);
          if(is.eof()) break;
        }
      }

      inline std::pair<std::vector<PostitionedOrganism*>,std::pair<double,double>> getNextFrame()
      {
        if(cF>=frames.size()) throw EndofFramesException();
        cF++;
        return make_pair(frames[cF-1],shapes[cF-1]);
      }
      inline std::pair<std::vector<PostitionedOrganism*>,std::pair<double,double>> getCurrentFrame()
      {
        cout << cF << " requested for frame count "<<this->frames.size()<<"\n";
        return make_pair(this->frames[cF],shapes[cF]);
      }
      inline void resetFrame() {cF=0;}

      void deduceShapes()
      {
        for(int i=0;i<frames.size();i++)
        {
          double mvw=1.0/0.0,mVw=-1.0/0.0;
          double mvh=1.0/0.0,mVh=-1.0/0.0;
          for(PostitionedOrganism* porg : frames[i])
          {
            if(mvw>porg->getPosition()(0)) mvw=porg->getPosition()(0);
            if(mVw<porg->getPosition()(0)) mVw=porg->getPosition()(0);
            if(mvh>porg->getPosition()(1)) mvh=porg->getPosition()(1);
            if(mVh<porg->getPosition()(1)) mVh=porg->getPosition()(1);
          }
          if(mvw>0) mvw=0; if(mvh>0) mvh=0;
          shapes.push_back(make_pair(mVw-mvw,mVh-mvh));
        }
      }
    };
    class EnvironmentRender : public graphics::Frame
    {
      std::map<std::string,int> colIDs;
      DrawingBoard brd;
      int lastId=-1;
      std::string radius_feature="Radius";
      double defrad=3.;
      bool pause=false;
      std::string logf="renderlog.log";
      std::ofstream of;
      std::vector<evolution::statistics::GPS*> stats;
    public:
      EnvironmentRender(int w,int h,const std::string& title="Environment",const std::string& lf="renderlog.log") : Frame(w,h,title) {this->autoclear=true; logf=lf; of.open(lf);}
      ~EnvironmentRender() {of.close();}

      EnvironmentRender& addStatistic(evolution::statistics::GeneralPositionalStatistic& gps) {stats.push_back(&gps); return *this;}
      void eraseMemory() {colIDs.empty();}
      void addDrawingBoard(DrawingBoard b) {brd=b;}
      inline void setDefaultRadius(double r) {defrad=r;}
      inline void setRadiusTrait(const std::string& s) {radius_feature=s;}

      void onMouseRelease(const sf::Event& e) override
      {
        cout << ((sf::Event::MouseButtonEvent&)e).x << ","<<((sf::Event::MouseButtonEvent&)e).y << "\n";
        pause=!pause;
      }
      void onDraw() override
      {
        if(pause)
        {
          auto recv=brd.getCurrentFrame();
          this->render(get<0>(recv),get<1>(recv),radius_feature,defrad);
        }
        else
        {
          try
          {
            auto recv=brd.getNextFrame();
            this->render(get<0>(recv),get<1>(recv),radius_feature,defrad);
            writeLog(of);
          }
          catch(EndofFramesException e) {return;}
        }
      }
      virtual void writeLog(std::ostream& os)
      {
        os << brd.cF << " ";
        std::vector<PostitionedOrganism*> orgs=get<0>(brd.getCurrentFrame());
        Environment n(1,orgs);
        for(const auto& ot : evolution::organisms)  os << n.getOrganismCount(get<0>(ot))<<" ";
        if(stats.size()) os << " | ";
        for(auto& stat : stats) {stat->compute(orgs); os << stat->mean << " "<<stat->std<<" "<<stat->other<<" ";}
        os << "\n";
        os.flush();
      }
      void render(std::vector<PostitionedOrganism*> state, std::pair<double,double> shape, std::string rad="Radius",double raddef=3.)
      {
        //this->clear();
        double scaleX=wid/get<0>(shape),scaleY=hei/get<1>(shape);
        for(PostitionedOrganism* porg : state)
        {
          double r=raddef;
          if(porg->hasTrait(rad)) r=porg->getTrait(rad);
          sf::Color c;
          if(colIDs.find(porg->getName())!=colIDs.end()) c=COLORS[colIDs[porg->getName()]];
          else
          {
            lastId++;
            colIDs.emplace(make_pair(porg->getName(),lastId));
            c=COLORS[lastId];
          }
          Eigen::Vector2d proj=porg->get2DProjection();
          //auto sp=objects::shapes::getCircle((float)r,c,c,proj(0)*scaleX,proj(1)*scaleY);
          sf::CircleShape sp((float)r);
          sp.setFillColor(c);
          sp.setOutlineColor(c);
          sp.setPosition(proj(0)*scaleX,proj(1)*scaleY);
          this->draw(sp);
        }
      }
    };
  }
  #endif
}

//evolution::traits=std::map<std::string,evolution::TraitSpecifics>();
#endif
