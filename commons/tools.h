#ifndef INCLUDED_TOOLS
#define INCLUDED_TOOLS 1
#include "commons.h"
#include <cassert>
namespace tools
{
  class UnboundString
  {
    std::string str;
  public:
    UnboundString() {str="";}
    UnboundString(int s) {str=""; for(int i=0;i<s;i++) str+=std::string(1,(char)0);}
    UnboundString(int s,char c) {str=""; for(int i=0;i<s;i++) str+=std::string(1,c);}
    UnboundString(char c) {str=std::string(1,c);}
    UnboundString(const char* c) {str=std::string(c);}
    UnboundString(const char* c,int s) {str=std::string(c,s);}
    UnboundString(const std::string& st) {str=st;}

    inline int length() const {return str.length();}
    inline int getLength() const {return str.length();}
    inline int getSize() const {return str.length();}

    inline char operator[](int i) const
    {
      if(i>=0 && i<str.length()) return str[i];
      else return (char)0;
    }
    inline char& charAt(int i)
    {
      if(i>=0 && i<str.length()) return str[i];
      else throw commons::IndexOutOfBoundsException(i,str.length());
    }
    inline int indexOf(char c) const
    {
      for(int i=0;i<str.length();i++) {if(str[i]==c) return i;}
      return -1;
    }
    inline UnboundString substring(int s,int e=-1) const
    {
      if(e==-1) e=str.length();
      assert(s>=0);
      return str.substr(s,e-s);
    }
    inline UnboundString substringByLength(int s,int l) const
    {
      assert(s>=0);
      return str.substr(s,l);
    }
    inline std::string toString() const {return str;}
    inline operator std::string() const {return toString();}
  };
}
static std::istream& operator>>(std::istream& is,tools::UnboundString& ubs)
{
  std::string temp;
  getline(is,temp);
  ubs=tools::UnboundString(temp);
  return is;
}
static std::ostream& operator<<(std::ostream& os,const tools::UnboundString& ubs)
{
  os << ubs.toString();
  return os;
}
#endif
