#ifndef INCLUDED_OPTIM
#define INCLUDED_OPTIM 1
#include "general/general.h"
#include "maths/numbers.h"
namespace optim
{
  class UnsolvedLinearOptimizationException : public std::exception {};
  template<class T=double> class LinearOptimizer
  {
  protected:
    GeneralVector<T> obj;
    std::vector<GeneralVector<T>> conds;
    GeneralVector<T> bounds;

  public:
    LinearOptimizer(std::vector<T> objective,std::vector<std::vector<T>> constr,std::vector<T> vals)
    {
      obj=GeneralVector<T>(objective);
      for(const auto& con : constr) conds.push_back(GeneralVector<T>(con));
      bounds=GeneralVector<T>(vals);
    }

    inline void setObjective(std::vector<T> objective) {obj=GeneralVector<T>(objective);}
    inline void addConstraint(const std::vector<T>& v,T l)
    {
      std::vector<T> bn=bounds;
      bn.push_back(l);
      bounds=GeneralVector<T>(bn);
      conds.push_back(GeneralVector<T>(v));
    }
    inline void addConstraints(std::vector<std::vector<T>> constr,std::vector<T> vals)
    {
      int K=0;
      for(const std::vector<T>& v : constr) addConstraint(v,vals[K++]);
    }
    inline void clearConstraints() {conds=std::vector<GeneralVector<T>>(); bounds=GeneralVector<T>(std::vector<T>());}
    inline const GeneralVector<T>& getObjective() const {return obj;}
    inline const std::vector<GeneralVector<T>>& getConstraints() const {return conds;}
    inline const GeneralVector<T>& getBoundaries() const {return bounds;}

    inline virtual GeneralVector<T> solve(bool quick=true,bool write=false,std::ostream& os=std::cout) const =0;

  };

  template<class T> class Maximize : public LinearOptimizer<T>
  {
  public:
    Maximize(std::vector<T> objective,std::vector<std::vector<T>> constr,std::vector<T> vals) : LinearOptimizer<T>(objective,constr,vals) {}

    GeneralVector<T> solve(bool quick=true,bool write=false,std::ostream& os=std::cout) const override;
  };
}

template<class T> inline static std::ostream& operator<<(std::ostream& os,const optim::LinearOptimizer<T>& f)
{
  os << "Max/Min (linear): ";
  for(int i=0;i<f.getObjective().size();i++) os <<"+ "<< f.getObjective()(i)<<"x"<<(i+1)<<" "; os<<"\n";
  os<< "s.t.\n";
  int L=0;
  for(const auto& cex : f.getConstraints())
  {
    int K=1;
    std::vector<T> c=cex;
    for(const T& coef : c)
    {
      if(coef!=0) os <<"+ "<< coef<<"x"<<K<<" ";
      K++;
    }
    os << "<?/?> "<<f.getBoundaries()(L++)<<"\n";
  }
  return os;
}

//Completing the simplex algorithms
template<class T> GeneralVector<T> optim::Maximize<T>::solve(bool quick,bool write,std::ostream& os) const
{
  os << *this << "\n\n"; T ZERO(0),ONE(1);
  std::vector<int> basic;
  std::vector<T> solvarstemp;
  std::vector<T> allvarstemp=this->obj; //for(T t: allvarstemp) solvarstemp.push_back(ZERO);
  int sI=allvarstemp.size(); for(int i=0;i<this->conds.size();i++) {allvarstemp.push_back(ZERO); basic.push_back(sI++); solvarstemp.push_back(this->bounds(i));}
  GeneralVector<T> objcoeff=allvarstemp;
  GeneralVector<T> solution=solvarstemp;
  std::vector<GeneralVector<T>> expandedconds;
  int appends=this->conds.size();
  int K=0;
  for(GeneralVector<T> cv : this->conds)
  {
    std::vector<T> temp=cv;
    for(int i=0;i<appends;i++)
    {
      if(i==K) temp.push_back(ONE);
      else temp.push_back(ZERO);
    }
    K++;
    expandedconds.push_back(temp);
  }
  bool found=true,found2=true;
  int iters=0;
  while(found)
  {
    iters++;
    found=false;
    if(write)
    {
      os << "\n\nIteration "<<iters<<" start\n";
      os << "-------------------\n";
      os << objcoeff.transpose() << "\n";
      int Kd=0;
      for(const auto& v : expandedconds) os << v.transpose() <<" = "<<solution(Kd++) <<"\n";
    }
    T mC; int ind=-1;
    for(int i=0;i<objcoeff.getSize();i++)
    {
      if(objcoeff(i)>ZERO)
      {
        if(!found)
        {
          found=true;
          mC=objcoeff(i);
          ind=i;
          if(!quick) break;
        }
        else
        {
          if(objcoeff(i)>mC)
          {
            mC=objcoeff(i);
            ind=i;
          }
        }
      }
    }
    if(!found) break;
    int lInd=-1; T lim;
    found2=false;
    cout << "Solution: "<<solution.transpose() << "\n";
    for(int i=0;i<expandedconds.size();i++)
    {
      if(expandedconds[i](ind)<=ZERO) continue;
      T est=solution(i)/expandedconds[i](ind);

      if(!found2)
      {
        found2=true;
        lim=est;
        lInd=i;
      }
      else
      {
        if(lim>est)
        {
          lim=est;
          lInd=i;
        }
      }
    }
    if(!found2) {cerr << "Error encountered variable. Entering variable present but no suitable leaving variable \n"; throw optim::UnsolvedLinearOptimizationException();}
    if(write) os << "Entering: "<<ind<<" and leaving: "<<basic[lInd]<<"\n";
    basic[lInd]=ind;
    T cf=expandedconds[lInd](ind); //os << expandedconds[lInd].transpose() << " as coefficients in row\n"; //os << cf << " as coefficient.\n";
    for(int i=0;i<expandedconds[lInd].getSize();i++) expandedconds[lInd](i)=expandedconds[lInd](i)/cf;
    solution(lInd)=solution(lInd)/cf;
    for(int i=0;i<expandedconds.size();i++)
    {
      if(i==lInd) continue;
      T sc=expandedconds[i](ind);
      expandedconds[i]=expandedconds[i]-expandedconds[lInd]*sc;
      solution(i)=solution(i)-solution(lInd)*sc;
    }
    T ocf=objcoeff(ind);
    objcoeff=objcoeff-expandedconds[lInd]*ocf;
    if(write)
    {
      os << "Iteration "<<iters<<" complete\n";
      os << "-------------------\n";
      os << objcoeff.transpose() << "\n";
      int Kd=0;
      for(const auto& v : expandedconds) os << v.transpose() <<" = "<<solution(Kd++) <<"\n";
    }
  }
  if(write) os << "All coefficients negative. Algorithm STOP.\n";
  os << "With vars: "<<basic<<"\n";
  os << "And solution: "<<solution.transpose()<<"\n";
  std::vector<T> finalsol;
  for(int i=0;i<expandedconds[0].getSize();i++)
  {
    int find=-1;
    for(int j=0;j<basic.size();j++)
    {
      if(basic[j]==i) {find=j; break;}
    }
    if(find!=-1) finalsol.push_back(solution(find));
    else finalsol.push_back(ZERO);
  }
  cout << finalsol << "\n";
  return GeneralVector<T>(finalsol);
}
#endif
