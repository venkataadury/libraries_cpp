#ifndef INCLUDED_THERMO
#define INCLUDED_THERMO 1
#include "sci/units.h"
#include "maths/maths.h"
using namespace thermodynamicunits;
namespace thermo
{
  class System;
  class NoSuchSystemException : public std::exception {};
  class InvalidPathException : public std::exception {};
  class CorruptSystemException : public std::exception {};
  namespace states
  {
    class InsufficientStateException : public std::exception {};
    class IncompleteStateException : public std::exception {};
    class InvalidStateException : public std::exception {};
    class NoMoleculesInSystemException : public InvalidStateException {};
    struct State
    {
      Pressure p;
      Volume V;
      Moles n;
      Temperature T;
      bool lp=false,lV=false,ln=false,lT=false;
      State() {}
      State(const Pressure& pe) {p=pe; lp=true;}
      State(const Volume& Ve) {V=Ve; lV=true;}
      State(const Moles& ne) {n=ne; ln=true;}
      State(const Temperature& Te) {T=Te; lT=true;}
      State(const Pressure& pe,const Volume& Ve,const Moles& ne) : State(pe) {V=Ve; n=ne; ln=true; lV=true;}
      State(const Pressure& pe,const Volume& Ve,const Temperature& Te) : State(pe) {V=Ve; T=Te; lT=true; lV=true;}
      State(const Pressure& pe,const Temperature& Te,const Moles& ne) : State(pe) {T=Te; n=ne; lT=true; ln=true;}
      State(const Volume& Ve,const Temperature& Te,const Moles& ne) : State(Ve) {T=Te; n=ne; lT=true; ln=true;}
      State(const Pressure& pe,const Volume& Ve,const Moles& ne,const Temperature& Te) : State(pe,Ve,ne) {T=Te;}

      inline bool hasPressure() const {return lp;}
      inline bool hasVolume() const {return lV;}
      inline bool hasMoles() const {return ln;}
      inline bool hasTemperature() const {return lT;}

      inline const Pressure& getPressure() const
      {
        if(!lp) throw IncompleteStateException();
        return p;
      }
      inline const Volume& getVolume() const
      {
        if(!lV) throw IncompleteStateException();
        return V;
      }
      inline const Temperature& getTemperature() const
      {
        if(!lT) throw IncompleteStateException();
        return T;
      }
      inline const Moles& getMoles() const
      {
        if(!ln) throw IncompleteStateException();
        return n;
      }

      inline void setPressure(const Pressure& pe) {p=pe; lp=true;}
      inline void setVolume(const Volume& Ve) {V=Ve; lV=true;}
      inline void setTemperature(const Temperature& Te) {T=Te; lT=true;}
      inline void setMoles(const Moles& ne) {n=ne; ln=true;}

      inline void unsetPressure() {lp=false;}
      inline void unsetVolume() {lV=false;}
      inline void unsetTemperature() {lT=false;}
      inline void unsetMoles() {ln=false;}
    };
    struct StateEquations
    {
      Pressure (*pf)(const Volume& V,const Moles& n,const Temperature& T)=nullptr;
      Volume (*Vf)(const Pressure& p,const Moles& n,const Temperature& T)=nullptr;
      Temperature (*Tf)(const Pressure& p,const Volume& V,const Moles& n)=nullptr;
      Moles (*nf)(const Pressure& p,const Volume& V,const Temperature& T)=nullptr;
      std::string label;
    public:
      //Order: p,V,T,n
      StateEquations(Pressure (*pif)(const Volume& V,const Moles& n,const Temperature& T),Volume (*Vif)(const Pressure& p,const Moles& n,const Temperature& T),Temperature (*Tif)(const Pressure& p,const Volume& V,const Moles& n),Moles (*nif)(const Pressure& p,const Volume& V,const Temperature& T),const std::string& eqlab="")
      {
        pf=pif;
        Vf=Vif;
        Tf=Tif;
        nf=nif;
        label=eqlab;
      }

      inline Pressure getPressure(const Volume& V,const Moles& n,const Temperature& T) const {return pf(V,n,T);}
      inline Volume getVolume(const Pressure& p,const Moles& n,const Temperature& T) const {return Vf(p,n,T);}
      inline Temperature getTemperature(const Pressure& p,const Volume& V,const Moles& n) const {return Tf(p,V,n);}
      inline Moles getMoles(const Pressure& p,const Volume& V,const Temperature& T) const {return nf(p,V,T);}
      inline const std::string& getLabel() const {return label;}
      inline void setLabel(const std::string& l) {label=l;}
      void completeState(State& s) const
      {
        if(s.hasPressure() && s.hasTemperature() && s.hasVolume() && s.hasMoles()) return;
        if(s.hasPressure() && s.hasVolume() && s.hasMoles()) s.setTemperature(getTemperature(s.getPressure(),s.getVolume(),s.getMoles()));
        else if(s.hasPressure() && s.hasVolume() && s.hasTemperature()) s.setMoles(getMoles(s.getPressure(),s.getVolume(),s.getTemperature()));
        else if(s.hasPressure() && s.hasMoles() && s.hasTemperature()) s.setVolume(getVolume(s.getPressure(),s.getMoles(),s.getTemperature()));
        else if(s.hasVolume() && s.hasMoles() && s.hasTemperature()) s.setPressure(getPressure(s.getVolume(),s.getMoles(),s.getTemperature()));
        else throw InsufficientStateException();
      }
      template<class SV1,class SV2,class SV3>
      inline State generateState(const SV1& a1,const SV2& a2,const SV3& a3) const {State s(a1,a2,a3); completeState(s); return s;}
    };
    namespace idealgas //pV=nRT, U=3/2nRT
    {
      static inline Pressure p(const Volume& V,const Moles& n,const Temperature& T) {return (n*R_gas*T)/V;}
      static inline Volume V(const Pressure& p,const Moles& n,const Temperature& T) {return (n*R_gas*T)/p;}
      static inline Temperature T(const Pressure& p,const Volume& V,const Moles& n) {return (p*V)/(n*R_gas);}
      static inline Moles n(const Pressure& p,const Volume& V,const Temperature& T) {return (p*V)/(R_gas*T);}
      static inline Energy U_mono(const State& s) {return 1.5*R_gas*s.getTemperature();}
      static inline Energy U_di(const State& s) {return 2.5*s.getMoles()*R_gas*s.getTemperature();}
      static inline Temperature T_mono(const Energy& U,const Pressure& p,const Volume& V,const Moles& n) {return (2/3.0)*U/(n*R_gas);}
      static inline Temperature T_di(const Energy& U,const Pressure& p,const Volume& V,const Moles& n) {return ((2/5.0)*U)/(n*R_gas);}
    }
    struct ThermodynamicEquations
    {
      Energy (*Uf)(const State& s)=nullptr;
      Temperature (*Tf)(const Energy& s,const Pressure& p,const Volume& V,const Moles& n)=nullptr;
    public:
      ThermodynamicEquations(Energy (*Uif)(const State& s),Temperature (*Tif)(const Energy& s,const Pressure& p,const Volume& V,const Moles& n)) {Uf=Uif; Tf=Tif;}

      Energy getInternalEnergy(const State& s) const {return Uf(s);}
      Temperature getTemperature(const Energy& s,const Pressure& p,const Volume& V,const Moles& n) const {return Tf(s,p,V,n);}
    };

    static const StateEquations STATE_IDEALGAS(idealgas::p,idealgas::V,idealgas::T,idealgas::n,"Ideal Gas Equations");
    static const ThermodynamicEquations THERMO_MONOATOMIC_IDEALGAS(idealgas::U_mono,idealgas::T_mono);
    static const ThermodynamicEquations THERMO_DIATOMIC_IDEALGAS(idealgas::U_di,idealgas::T_di);
    static const states::State STP(atmos_press,Volume(0.0224127),zero);
  }
  static Pressure (*pf)(const Volume& V,const Moles& n,const Temperature& T)=nullptr;
  class Component //A component in a thermodynamic system
  {
    states::StateEquations equations;
    std::string name="Unknown";
    states::ThermodynamicEquations thermoeq;
    Moles mol;
    bool permeable=true;

  public:
    Component(const Moles& amt,const states::StateEquations& eq=states::STATE_IDEALGAS,const states::ThermodynamicEquations& teq=states::THERMO_MONOATOMIC_IDEALGAS) : equations(eq),thermoeq(teq) {mol=amt;}
    Component(const Moles& amt,const std::string& n,const states::StateEquations& eq=states::STATE_IDEALGAS,const states::ThermodynamicEquations& teq=states::THERMO_MONOATOMIC_IDEALGAS) : Component(amt,eq,teq) {name=n;}
    //Component(const Moles& amt,const std::string& n,const System& sys,const states::StateEquations& eq=states::STATE_IDEALGAS,const states::ThermodynamicEquations& teq=states::THERMO_MONOATOMIC_IDEALGAS);

    inline const Moles& getMoles() const {return mol;}
    inline void setMoles(double n) {mol=Moles(n);}
    inline void setPermeable(bool b=true) {permeable=b;}
    inline bool isPermeable() const {return permeable;}
    inline const std::string& getName() const {return name;}
    inline bool operator==(const Component& c2) const {return (name==c2.name);}
    inline const states::StateEquations& getStateEquations() const {return equations;}
    inline const states::ThermodynamicEquations& getThermodynamicEquations() const {return thermoeq;}

    Volume getOccupiedVolume(const states::State& v) const //Automatically omit reading vol,n from state
    {
      states::State stt=v; stt.unsetVolume(); stt.setMoles(mol);
      equations.completeState(stt);
      return stt.getVolume();
    }
    Pressure getExertedPressure(const states::State& v) const //Automatically omit reading p,n from state
    {
      states::State stt=v; stt.unsetPressure(); stt.setMoles(mol);
      equations.completeState(stt);
      return stt.getPressure();
    }
    Temperature getExpectedTemperature(const states::State& v) const //Automatically omit reading T,n from state
    {
      states::State stt=v; stt.unsetTemperature(); stt.setMoles(mol);
      equations.completeState(stt);
      return stt.getTemperature();
    }

    static inline Component createMonoatomicIdealGas(const std::string& n,const Moles& m=Moles(1)) {return Component(m,n,states::STATE_IDEALGAS,states::THERMO_MONOATOMIC_IDEALGAS);}
    static inline Component createDiatomicIdealGas(const std::string& n,const Moles& m=Moles(1)) {return Component(m,n,states::STATE_IDEALGAS,states::THERMO_DIATOMIC_IDEALGAS);}
    /*inline void setState(const states::State& s) {mystate=s;}
    inline const states::State& getState() const {return mystate;}*/

  };
  class PolytropicPath;
  //inline static Volume workIntegral(const Path& path,const Volume& vi,const Volume& vf);
  class PolytropicPath //Entails the definition of a polytropic process
  {
    //PV^r const
    double r;
  public:
    PolytropicPath() {r=0;} //Isobaric
    PolytropicPath(double ri) {r=ri;}

    DerivedUnit intf(Volume v) const {Volume volr; volr.setDefault(pow(v,-r)); return volr;}
    Volume workIntegral(const Volume& ll,const Volume& ul,const Volume& acc=Volume(1e-5)) const
    {
      if(::abs((double)(intf(ul)-intf(ll)))<acc) return intf((ul+ll)/2)*(ul-ll);
      return workIntegral(ll,(ll+ul)/2,acc)+workIntegral((ll+ul)/2,ul,acc);
    }

    Work getWork(const Pressure& p1,const Pressure& p2,const Volume& v1,const Volume& v2=Volume(0)) const
    {
      Volume vf;
      if(v2!=0)
      {
        //cout << ::abs(((double)(p1*pow(v1,r)))-((double)(p2*pow(v2,r)))) << "\n";
        if(::abs(((double)(p1*pow(v1,r)))-((double)(p2*pow(v2,r))))>1e-4) throw InvalidPathException();
        vf=v2;
      }
      else vf=pow(p1*maths::pow(v1,r)/p2,1/r);
      return (-p1*pow(v1,r)*workIntegral(v1,vf));
    }
  };
  //inline static Volume workIntegral(const Path& path,const Volume& vi,const Volume& vf) {return integrate<Volume,DerivedUnit>(&thermo::Path::intf,vi,vf);}
  typedef PolytropicPath Process;
  inline static Process generateIsobaricProcess() {return Process();}
  inline static Process generateIsochoricProcess() {return Process(1.0/0.0);}
  inline static Process generateIsothermalProcess(const states::StateEquations& eq,const states::State& s)
  {
    if(!s.hasTemperature()) throw InvalidPathException();
    states::State ts=s;
    Pressure p1; Volume v1; Temperature T=ts.getTemperature();
    if(ts.hasPressure()) p1=ts.getPressure();
    else p1.setDefault(1);
    if(ts.hasVolume()) v1=ts.getVolume();
    else v1.setDefault(22.4);
    ts.setVolume(v1); ts.setPressure(p1);
    if(!ts.hasMoles()) eq.completeState(ts);
    states::State s2(Volume(v1*2),ts.getTemperature(),ts.getMoles()); s2.unsetPressure();
    eq.completeState(s2);
    cout << p1/s2.getPressure() << "\t" << s2.getVolume()/v1 << "\n";
    double r=log(p1/s2.getPressure())/log(s2.getVolume()/v1);
    cout << "Got r = "<<r<<"\n";
    return Process(r);
  }
  inline static Process generateGuessedProcess(states::State s1,states::State s2,const states::StateEquations* steq=nullptr)
  {
    Pressure p1,p2; Volume v1,v2;
    if(s1.hasPressure()) p1=s1.getPressure();
    else {if(!steq) throw InvalidPathException(); else steq->completeState(s1); p1=s1.getPressure();}
    if(s2.hasPressure()) p2=s2.getPressure();
    else {if(!steq) throw InvalidPathException(); else steq->completeState(s2); p2=s2.getPressure();}
    if(s1.hasVolume()) v1=s1.getVolume();
    else {if(!steq) throw InvalidPathException(); else steq->completeState(s1); v1=s1.getVolume();}
    if(s2.hasVolume()) v2=s2.getVolume();
    else {if(!steq) throw InvalidPathException(); else steq->completeState(s2); v2=s2.getVolume();}
    double ratv=v2/v1,ratp=p1/p2;
    double r= log(ratp)/log(ratv);
    cout << "Guessed r="<<r<<"\n";
    return Process(r);
  }
  namespace pathfunctions
  {
    class Path
    {
      Energy Q; Work w;
      bool qset=false,wset=false,undet=false;
      const states::ThermodynamicEquations* syseq=nullptr;
      const states::StateEquations* steq=nullptr;
      states::State si,sf;
    public:
      Path() {undet=true;}
      Path(const states::State& s1,const states::State& s2,const states::ThermodynamicEquations* eq=nullptr,bool guess=true)
      {
        si=s1;
        sf=s2;
        syseq=eq;
        if(guess) guessPath();
      }
      Path(const states::State& s1,const states::State& s2,const states::ThermodynamicEquations& eq,const Energy& Qproc) : Path(s1,s2,&eq,false)
      {
        Q=Qproc;
        qset=true;
        solvePath();
      }
      Path(const states::State& s1,const states::State& s2,const Process& prc,const states::StateEquations* seq=nullptr)
      {
        si=s1;
        sf=s2;
        if(seq) addStateEquations(*seq);
        w=prc.getWork(si.getPressure(),sf.getPressure(),si.getVolume(),sf.getVolume());
        wset=true;
        solvePath();
      }

      void addStateEquations(const states::StateEquations& seq)
      {
        steq=&seq;
        steq->completeState(si);
        steq->completeState(sf);
      }

      void solvePath()
      {
        if(qset && wset) return;
        Energy dU(syseq->getInternalEnergy(sf)-syseq->getInternalEnergy(si));
        if(qset) {w=Work(dU-Q); wset=true;}
        else {Q=Energy(dU-w); qset=true;}
      }
      void guessPath()
      {
        Process gP=generateGuessedProcess(si,sf,steq);
        w=gP.getWork(si.getPressure(),sf.getPressure(),si.getVolume(),sf.getVolume());
        wset=true;
        solvePath();
      }

      inline bool isUndetermined() const {return undet;}
    } UNDETERMINEDPATH;
    inline static Path generateIsobaricProcess(const states::State& s1,const states::State& s2) {return Path(s1,s2);}
    inline static Path generateIsochoricProcess(const states::State& s1,const states::State& s2) {return Path(s1,s2);}
    inline static Path generateIsothermalProcess(const states::State& s1,const states::State& s2,const states::StateEquations& nse) {return Path(s1,s2,thermo::generateIsothermalProcess(nse,s1),&nse);}
    inline static Path generateAdiabaticProcess(const states::State& s1,const states::State& s2,const states::ThermodynamicEquations& teq) {return Path(s1,s2,teq,Energy(0));}
  }

  class System
  {
    std::vector<Component> components;
    states::State state;
    bool vpar=0; //0 - vol, 1-press
    std::string label="System";
  public:
    System(const states::State& s=states::STP) {state=s;}
    System(const std::string& n,const states::State& s=states::STP) : System(s) {label=n;}

    inline const std::string& getLabel() const {return label;}
    inline void setLabel(const std::string& s) {label=s;}
    inline void changePressureOnInsertion(bool c=true) {vpar=c;}
    inline void changeVolumeOnInsertion(bool c=true) {vpar=!c;}
    void addComponent(const Component& c)
    {
      bool found=false;
      for(int i=0;i<components.size();i++){if(components[i]==c) {components[i].setMoles(components[i].getMoles()+c.getMoles()); found=true; break;}}
      if(!found) components.push_back(c);
      adjustState(vpar);
    }
    void adjustState(bool varpar=0)
    {
      Volume nv(0); //=state.getVolume();
      Pressure np(0); //=state.getPressure();
      if(!varpar)
      {
          for(const Component& c : components) {nv=nv+c.getOccupiedVolume(state);}
          state.setVolume(nv);
      }
      else
      {
        for(const Component& c : components) {np=np+c.getExertedPressure(state);}
        state.setPressure(np);
      }
    }
    inline const std::vector<Component>& getComponents() const {return components;}
    inline std::vector<Component>& getComponents() {return components;}

    bool hasComponent(const Component& c) const
    {
      for(const Component& ci : components) {if(ci==c) return true;}
      return false;
    }

    inline Temperature getTemperature() const {return state.getTemperature();}
    inline Volume getVolume() const {return state.getVolume();}
    inline Pressure getPressure() const {return state.getPressure();}
    inline Moles getMolesOf(const Component& c) const {for(const Component& cmp : components) {if(c==cmp) return cmp.getMoles();} return Moles(0);}
    inline Moles getAllMoles() const {Moles ret(0); for(const Component& cmp : components) ret=ret+cmp.getMoles(); return ret;}
    inline states::State getState() const {return state;}

    inline void setTemperature(const Temperature& T) {state.setTemperature(T);}
    inline void setPressure(const Pressure& p) {state.setPressure(p);}
    inline void setVolume(const Volume& V) {state.setVolume(V);}
    inline Energy getInternalEnergy() const {return getInternalEnergyAt(this->getTemperature());}
    Energy getInternalEnergyAt(Temperature T) const
    {
      Energy ret(0);
      for(const Component& cmp : components)
      {
        states::State stt=state;
        stt.setMoles(cmp.getMoles());
        stt.setTemperature(T);
        ret=ret+cmp.getThermodynamicEquations().getInternalEnergy(stt);
      }
      return ret;
    }
    inline Pressure getPressureAt(Volume V) const
    {
      Pressure ret(0);
      for(const Component& cmp : components)
      {
        states::State stt=state; stt.setVolume(V);
        stt.setMoles(cmp.getMoles());
        ret=ret+cmp.getExertedPressure(stt);
      }
      return ret;
    }
  };

  static states::State STATE_SURR(atmos_press,Volume(1.0/0.0),zero);
  static System ENV("Environment",STATE_SURR);

  class Wall
  {
    bool rigid=true,diathermic=false,open=false;
  public:
    Wall(bool r=true,bool d=false,bool o=false)
    {
      rigid=r; diathermic=d; open=o;
      if(open) {rigid=true; diathermic=true;}
    }

    inline bool isAdiabatic() const {return !diathermic;}
    inline bool isDiathermic() const {return diathermic;}
    inline bool isRigid() const {return rigid;}
    inline bool isOpen() const {return open;}
    inline bool isClosed() const {return !open;}
  };
  namespace stdwalls
  {
    class AdiabaticWall : public Wall
    {
    public:
      AdiabaticWall(bool r=false) : Wall(r,false,false) {}
    };
    class DiathermalWall : public Wall
    {
    public:
      DiathermalWall(bool r=false) : Wall(r,true,false) {}
    };
    class OpenWall : public Wall {public: OpenWall() : Wall(true,true,true) {}};
    typedef OpenWall NoWall;
    const AdiabaticWall ISOLATINGWALL(true);

  }
  class Environment
  {
    std::vector<System*> systems;
    std::vector<std::vector<const Wall*>> interactions; //0 - Nothing (isolated), 1 - Exchange heat, 2-> exchange both
  public:
    Environment() {systems.push_back(&ENV); interactions.push_back(std::vector<const Wall*>(1,&stdwalls::ISOLATINGWALL));}

    inline std::vector<System*> getSystems() const {return systems;}
    inline std::vector<std::vector<const Wall*>> getWalls() const {return interactions;}
    void addSystem(System& s,std::vector<const Wall*> intv=std::vector<const Wall*>())
    {
      std::vector<const Wall*> ng;
      if(intv.size()) {for(int i=0;i<systems.size();i++) {if(!intv[i]) intv[i]=&stdwalls::ISOLATINGWALL; ng.push_back(intv[i]); interactions[i].push_back(intv[i]);}}
      else {for(int i=0;i<systems.size();i++) {ng.push_back(&stdwalls::ISOLATINGWALL); interactions[i].push_back(&stdwalls::ISOLATINGWALL);}}
      for(int i=0;i<systems.size();i++)
      {
        for(Component& cmp : systems[i]->getComponents()) {if(!s.hasComponent(cmp)) s.addComponent(Component(Moles(0),cmp.getName(),cmp.getStateEquations(),cmp.getThermodynamicEquations()));}
        for(Component& cmp : s.getComponents()) {if(!systems[i]->hasComponent(cmp)) systems[i]->addComponent(Component(Moles(0),cmp.getName(),cmp.getStateEquations(),cmp.getThermodynamicEquations()));}
      }
      ng.push_back(&stdwalls::ISOLATINGWALL);
      systems.push_back(&s);
      interactions.push_back(ng);
    }
    void removeSystem(const System& s)
    {
      int i=0;
      for(;i<systems.size();i++) {if(systems[i]==&s) break;}
      if(i>=systems.size()) throw NoSuchSystemException();
      systems.erase(systems.begin()+i);
      interactions.erase(interactions.begin()+i);
      for(int j=0;j<interactions.size();j++) interactions[j].erase(interactions[j].begin()+i);
    }
    inline int indexOf(const System* s) const
    {
      for(int i=0;i<systems.size();i++) {if(systems[i]==s) return i;}
      return -1;
    }
    void setWallBetween(const System* s1,const System* s2,const Wall* wall)
    {
      int i1=indexOf(s1),i2=indexOf(s2);
      if(i1==-1 || i2==-1) throw NoSuchSystemException();
      interactions[i1][i2]=wall;
      interactions[i2][i1]=wall;
    }

    void equilibiriate(double ptol=2.5e-5,double vtol=1e-3,double ttol=0.1,double ntol=0.001)
    {
      bool changed=true;
      while(changed)
      {
        changed=false;
        for(int i=0;i<systems.size();i++)
        {
          for(int j=i+1;j<systems.size();j++)
          {
            if(interactions[i][j]->isDiathermic())
            {
              //std::cout << "Inital temperatures:\t"<<systems[i]->getTemperature()<<"K at "<<i<<", "<<systems[j]->getTemperature()<<"K at "<<j<<"\n";
              Temperature T1=systems[i]->getTemperature(),T2=systems[j]->getTemperature();
              Temperature deltaT=T1-T2;
              if(::abs(deltaT)<ttol) {}
              else
              {
                //dT/dU = 1/(dU/dT).
                double der1=1/((systems[i]->getInternalEnergyAt(systems[i]->getTemperature()+ttol*0.1)-systems[i]->getInternalEnergy())/(ttol*0.1));
                double der2=1/((systems[j]->getInternalEnergyAt(systems[j]->getTemperature()+ttol*0.1)-systems[j]->getInternalEnergy())/(ttol*0.1));
                if(der1!=der1) der1=0;
                if(der2!=der2) der2=0;
                double q=deltaT/(der1+der2);
                systems[i]->setTemperature(systems[i]->getTemperature()-der1*q);
                systems[j]->setTemperature(systems[j]->getTemperature()+der2*q);
                Temperature T1d=systems[i]->getTemperature(),T2d=systems[j]->getTemperature();
                if(T1!=T1d) systems[i]->adjustState(1);
                if(T2!=T2d) systems[j]->adjustState(1);
                //std::cout <<"Adjusted:"<<systems[i]->getPressure()<<","<<systems[j]->getPressure()<<"\n";
                changed=true;
              }
            }
            if(!interactions[i][j]->isRigid())
            {
              Pressure p1=systems[i]->getPressure(),p2=systems[j]->getPressure();
              Volume V1=systems[i]->getVolume(),V2=systems[j]->getVolume();
              Pressure deltap=Pressure(p1-p2);
              //std::cout << i<<","<<j<<":\t"<<p1<<","<<p2<<"\t"<<deltap << "\n";
              double cc=((deltap>0)?(double)deltap:-(double)deltap);
              if(cc<ptol) continue;
              //dp/dV=?
              FreeDerived4 der1=(systems[i]->getPressureAt(systems[i]->getVolume()+Volume(vtol))-p1)/(Volume(vtol));
              FreeDerived4 der2=(systems[j]->getPressureAt(systems[j]->getVolume()+Volume(vtol))-p2)/(Volume(vtol));
              if(der1!=der1) der1.setDefault(0);
              if(der2!=der2) der2.setDefault(0);
              Volume dV=deltap/(der1+der2);
              /*systems[i]->setPressure(systems[i]->getPressure()-der1*dV);
              systems[j]->setPressure(systems[j]->getPressure()+der2*dV);*/
              systems[i]->setVolume(systems[i]->getVolume()-dV);
              systems[j]->setVolume(systems[j]->getVolume()+dV);
              //Pressure p1d=systems[i]->getPressure(),p2d=systems[j]->getPressure();
              Volume V1d=systems[i]->getVolume(),V2d=systems[j]->getVolume();
              if(V1!=V1d) systems[i]->adjustState(1);
              if(V2!=V2d) systems[j]->adjustState(1);
              changed=true;
            }
          }
        }
      }
    }
  };


  inline static void setAirEnvironment()
  {
    ENV.addComponent(Component(Moles(1.0/0.0),"Air",states::STATE_IDEALGAS,states::THERMO_DIATOMIC_IDEALGAS));
  }
  inline static void setEnvironmentTemperature(const Temperature& tp)
  {
    ENV.setTemperature(tp);
  }

  namespace commonprocs
  {
    enum StandardPath {UNDETERMINED=0, ADIABATIC,ISOTHERMAL,ISOCHORIC,ISOBARIC};
    namespace closed
    {
      enum ChangingVariable {VOLUME,PRESSURE,TEMPERATURE};
      enum ChangeType {CONSTANT,ADDITIVE,MULTIPLICATIVE};
      static pathfunctions::Path changeSystem(System& sys,const states::State& fstate,StandardPath path=UNDETERMINED)
      {
        if(!path) return thermo::pathfunctions::UNDETERMINEDPATH;

      }
    }
  }
}
using namespace thermo::commonprocs::closed;
static std::ostream& operator<<(std::ostream& os,const thermo::states::State& s)
{
  os << "State description (Begin)\n";
  os << "\tPressure: ";
  if(s.hasPressure()) os << (double)s.getPressure() << " atm\n";
  else os << "Not defined\n";
  os << "\tVolume: ";
  if(s.hasVolume()) os << (double)s.getVolume() << " L\n";
  else os << "Not defined\n";
  os << "\tMoles: ";
  if(s.hasMoles()) os << (double)s.getMoles() << " moles\n";
  else os << "Not defined\n";
  os << "\tTemperature: ";
  if(s.hasTemperature()) os << (double)s.getTemperature() << " K\n";
  else os << "Not defined\n";
  os << "State description (End)\n";
  return os;
}
static std::ostream& operator<<(std::ostream& os,const thermo::System& s)
{
  os << "System description (Begin)\n";
  os << "\tPressure: ";
  os << (double)s.getPressure() << " atm\n";
  os << "\tVolume: ";
  os << (double)s.getVolume() << " L\n";
  os << "\tTemperature: ";
  os << (double)s.getTemperature() << " K\n";
  os << "\tComponents:\n";
  for(const thermo::Component& cmp : s.getComponents())
  {
    if(!cmp.getMoles()) continue;
    os << "\t\t"<<cmp.getName();
    if((cmp.getMoles()-cmp.getMoles())!=0) {os <<" (infinite)\n"; return os;}
    if(cmp.getStateEquations().getLabel().length()) os << " (following "<<cmp.getStateEquations().getLabel()<<")";
    os<<": "<<cmp.getMoles().getVal()<<" mol\n";
  }
  os << "System description (End)\n";
  return os;
}
static std::ostream& operator<<(std::ostream& os,const thermo::Wall& w)
{
  os << "Wall(";
  if(w.isOpen()) return (os << "Open/permeable)");
  else os << "Non-permeable,";
  os << (w.isRigid()?"Rigid,":"Free,");
  return (os << (w.isAdiabatic()?"Adiabatic)":"Diathermic)"));
}
static std::ostream& operator<<(std::ostream& os,const thermo::Environment& s)
{
  os << "----------------------------------------\nEnvironment description:\n";
  for(thermo::System* sp : s.getSystems()) os << (*sp)<<"\n";
  os << "\n\nWalls:\n";
  auto intr=s.getWalls();
  int sysrefs[500]; for(int i=0;i<500;i++) sysrefs[i]=0;
  unsigned int K=1;
  for(int i=0;i<intr.size();i++)
  {
    if(s.getSystems()[i]->getLabel()!="System") os << s.getSystems()[i]->getLabel()<<" with:\n";
    else
    {
      if(sysrefs[i]==0) sysrefs[i]=K++;
      os << "System "<<(sysrefs[i])<< " with:\n";
      //os <<"System "<<(i+1)<<" with:\n";
    }
    for(int j=0;j<intr[i].size();j++)
    {
      os << "\t";
      if(s.getSystems()[j]->getLabel()!="System") os << s.getSystems()[j]->getLabel()<<":\t";
      else
      {
        if(sysrefs[j]==0) sysrefs[j]=K++;
        os << "System "<<(sysrefs[j])<< ":\t";
      }
      if(intr[i][j]) os<<(*intr[i][j])<<"\n";
      else os << "N/A\n";
    }
    os << "\t.\n";
  }
  return (os << "----------------------------------------\n");
}
#endif
