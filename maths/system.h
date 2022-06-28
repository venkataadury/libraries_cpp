#ifndef INCLUDED_SYSTEM
#define INCLUDED_SYSTEM 1
//#include "commons/parsers.h"
#include "commons/commons.h"
#include "maths/maths.h"

namespace mathsystem
{
  class SuperFrame;
  static std::ostream& operator<<(std::ostream& os,const SuperFrame& sf);
  class SuperFrame
  {
    void* data;
    int typeID=-1;
    int subID=-1;
    //Type IDs: 0=double,1=vector,2=matrix,3.x=geometric figures
  public:
    SuperFrame() {data=nullptr;}
    SuperFrame(double x) {typeID=0; data=new double(x);}
    SuperFrame(maths::Vector x) {typeID=1; data=new Vector(x);}
    SuperFrame(maths::Matrix x) {typeID=2; data=new Matrix(x);}
    SuperFrame(geom3D::Point3D x) {typeID=3;subID=0; data=new geom3D::Point3D(x);}

    bool sameTypeAs(const SuperFrame& f2) const {return (typeID==f2.typeID && subID==f2.subID);}
    friend std::ostream& operator<<(std::ostream& os,const SuperFrame& sf);
    const int& getTypeID() const {return typeID;}
    /*bool operator<(const SuperFrame& sf)
    {

    }*/
  };
  static std::ostream& operator<<(std::ostream& os,const SuperFrame& sf)
  {
    switch(sf.typeID)
    {
      case 0:
        return os<<*((double*)sf.data);
      case 1:
        return os<<*((maths::Vector*)sf.data);
      case 2:
        return os<<*((maths::Matrix*)sf.data);
      case 3:
        return os<<*((geom3D::Point3D*)sf.data);
      default:
        return os;
    }
  }

  class Variable
  {
  public:
    String name;
    SuperFrame* data;

  public:
    Variable() {name="";data=nullptr;}
    Variable(const String& str) {name=str;data=nullptr;}
    Variable(const String& str,SuperFrame* d) : Variable(str) {data=d; cout << "Con1\n";}
    Variable(const String& str,SuperFrame& d) : Variable(str) {data=&d; cout << "Con2\n";}
    Variable(const String& str,SuperFrame d) : Variable(str) {data=&d; data=new SuperFrame(d);}

    SuperFrame& getData() {return *data;}
    const SuperFrame& getData() const {return *data;}
    String getName() const {return name;}

    //Operator overloading
    bool operator>(const Variable& v) const
    {
      if(v.name.getLength()==name.getLength()) return (v.name<name);
      else return (v.name.getLength()<name.getLength());
    }
    bool operator<(const Variable& v) const
    {
      if(v.name.getLength()==name.getLength()) return (v.name>name);
      else return (v.name.getLength()>name.getLength());
    }

    Variable& operator+(const Variable& v2) const
    {
      if(getData().sameTypeAs(v2.getData())) {throw exception();}
        //throw TypeMismatchException();
      switch(getData().getTypeID())
      {
        case 0:
          return *(new Variable("res",*((double*)data)+*((double*)v2.data)));
        case 1:
          return *(new Variable("res",*((maths::Vector*)data)+*((maths::Vector*)v2.data)));
        case 2:
          return *(new Variable("res",*((maths::Matrix*)data)+*((maths::Matrix*)v2.data)));
        case 3:
          return *(new Variable("res",*((geom3D::Point3D*)data)+*((geom3D::Point3D*)v2.data)));
        default:
          return *(new Variable("res"));
      }
    }

  };

  //static operator<(const Variable* v1,const Variable* v2) {return (*v1<*v2);}
  class MathSystem
  {
  public:
    std::vector<Variable> vars;
    bool sorted=false;

    void addVar(const Variable& var)
    {
      for(int i=0;i<vars.size();i++)
      {
        if(vars[i].getName()==var.getName())
        {
          vars[i].data=var.data;
          return;
        }
      }
      vars.push_back(var);
      //cout << "Post construction: "<<vars[vars.size()-1].getData() << "\n";
      sorted=false;
    }
    void sortVars()
    {
      if(sorted) return;
      vars=sort(vars);
      sorted=true;
    }
    void replaceVars(String& str)
    {
      //sortVars();
      std::string strv;
      for(int i=vars.size()-1;i>=0;i--)
      {
        std::ostringstream oss;
        oss << vars[i].getData(); //Not to stdout
        strv=oss.str();
        str=str.replace(vars[i].getName(),strv);
      }
    }
  };
}

static std::ostream& operator<<(std::ostream& os,const mathsystem::Variable& v)
{
  os << v.getData();
  return os;
}
#endif
