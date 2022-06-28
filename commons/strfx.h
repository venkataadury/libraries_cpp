#ifndef STRINGS
#define STRINGS 1
#include "commons/commons.h"

//Provides some easy methods to use with strings.
// ### NO NAMESPACE ### //

using commons::String;
int getMatchingBracket(const String& str,char ob,char cb,int obi)
{
	if(str[obi]!=ob)
		return -1; //Invalid
	int tc=0;
	for(int i=obi+1;i<str.getLength();i++)
	{
		if(str[i]==cb)
		{
			tc--;
			if(tc==-1)
				return i;
		}
		else if(str[i]==ob)
			tc++;
	}
	return -1;
}

String getBracketContent(const String& str,char ob,char cb,int obi) //Returns without the brackets
{
	if(str[obi]!=ob)
		return nullptr; //Invalid
	int eI=getMatchingBracket(str,ob,cb,obi);
	return str.substring(obi+1,eI);
}
String* extractBracketContent(const String& str,char ob,char cb,int obi)
{
	if(str[obi]!=ob)
		return nullptr; //Invalid
	int eI=getMatchingBracket(str,ob,cb,obi);
	std::cout << eI << "\n";
	String* ret=new String[2];
	ret[0]=*(new String(str.substring(obi,eI+1)));
	ret[1]=*(new String(str.substring(0,obi)));
	ret[1]=ret[1]+*(new String(str.substring(eI+1)));
	return ret;
}
inline bool iswhitespace(char ch) {return (ch=='\n' || ch=='\r' || ch==' ' || ch=='\t');}
inline bool isCaps(char ch) {return (ch>=65 && ch<=90);}
inline bool isUpper(char ch) {return isCaps(ch);}
inline bool isLower(char ch) {return (ch>=97 && ch<=122);}
inline bool isAlphabet(char ch) {return (isCaps(ch) || isLower(ch));}
inline bool isNumeric(char ch) {return (ch>=48 && ch<=57);}
inline bool isBracket(char ch) {return(ch=='{' || ch=='}' || ch=='(' || ch==')' || ch=='[' || ch==']');} //not counting <>
inline char toUpper(char ch)
{
	if(isLower(ch))
		return ch-32;
	return ch;
}
inline char toLower(char ch)
{
	if(isUpper(ch))
		return ch+32;
	return ch;
}
#endif
