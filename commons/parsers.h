#ifndef INCLUDED_PARSERS
#define INCLUDED_PARSERS 1
#include "commons/commons.h"
#include "commons/strfx.h"
#include "maths/maths.h"
#include "maths/system.h"
using commons::String;
using mathsystem::Variable;
namespace parsers
{
  static long int curcount=0l;
  class ImproperStringArithmeticException : public exception
  {
    String str=String("String unknown");
  public:
    ImproperStringArithmeticException() {}
    ImproperStringArithmeticException(const String& st) {str=st;}

    virtual const char* what() throw() {return ("Bad string: "<<(str.toString())).c_str();}

  };
  const int OPERATOR_COUNT=6;
  const char BRACKET_OPEN='{',BRACKET_CLOSE='}';
  const char OPERATORS[OPERATOR_COUNT]={'^','/','%','*','+','-'};

  static bool checkParsable(const String& str)
  {
    for(char ch : str.toString())
    {
      if(isdigit(ch) || iswhitespace(ch) || ch==BRACKET_OPEN || ch==BRACKET_CLOSE || ch=='.' || indexOf(ch,OPERATORS,OPERATOR_COUNT)!=-1)
        continue;
      return false;
    }
    return true;
  }
  static double stringArithmetic(const String& str);
  static void processOperation(String& s,char op,int ind)
  {
    String o1=s.substring(0,ind).trim(),o2=s.substring(ind+1).trim();
    if(op=='-' && (o1.isEmpty() || !isdigit(o1[o1.size()-1])))
      return;
    double n1=stringArithmetic(o1),n2=stringArithmetic(o2);
    if(isnan(n1) || isnan(n2))
      throw ImproperStringArithmeticException();
    switch(op)
    {
      case '^':
        s=String("")<<(pow(n1,n2));
        return;
      case '%':
        s=String("")<<(((int)(n1*1e6)%(int)(n2*1e6))/1e6);
        return;
      case '/':
        s=String("")<<(n1/n2);
        return;
      case '*':
        s=String("")<<(n1*n2);
        return;
      case '+':
        s=String("")<<(n1+n2);
        return;
      case '-':
        s=String("")<<(n1-n2);
        return;
    }
  }
  static double stringArithmetic(const String& stri)
  {
      String str=stri;
      int ind,ind2;
      ind=str.indexOf(BRACKET_OPEN);
      if(ind!=-1)
      {
        String ss=getBracketContent(str,'{','}',ind);
        return stringArithmetic(str.replace(str.substring(ind,ss.size()+ind+2),String((std::string("")<<stringArithmetic(ss)))));
      }
      int opi=OPERATOR_COUNT-1;
      char cop;
      for(int i=opi;i>=0;i--)
      {
        cop=OPERATORS[i];
        ind=str.indexOf(cop);
        if(ind==-1) continue;
        if(cop=='-')
        {
          ind2=str.indexOf('+');
          if(ind2>ind) {ind=ind2;cop='+';}
        }
        else if(cop=='*')
        {
          ind2=str.indexOf('/');
          if(ind2>ind) {ind=ind2;cop='/';}
        }
        processOperation(str,cop,ind);
      }
      return std::stod(str.trim());
  }
}
#endif
