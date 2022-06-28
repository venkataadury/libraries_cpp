#ifndef INCLUDED_X
#define INCLUDED_X 1

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <commons/flex.h>
/*
 * Commons is a package for most common functions including basic logic, simple strings, file io,etc
 */
//Implementing complex java analogies
class Object; //Refer end
//Simple overrides !!!Note; this overrides << operator for the std::string class when used with integers and floats!!!
static inline std::string operator<<(std::string str,int num) {return std::string(str+std::to_string(num));}
static inline std::string operator<<(int num,std::string str) {return std::string(std::to_string(num)+str);}
static inline std::string operator<<(std::string str,double num) {return std::string(str+std::to_string(num));}
static inline std::string operator<<(double num,std::string str) {return std::string(std::to_string(num)+str);}
static inline std::string operator<<(std::string s1,std::string s2) {return s1+s2;}

//static inline std::ostream& operator<<(std::ostream& os,bool bl) {os << ((bl)?"true":"false");}
//static inline std::string operator<<(const char* ptr,double num) {return operator<<(std::string(ptr),num);}

//Common access functions (And template functions)
//template<class T> static inline String toString(T obj);
static inline bool isPrintable(char ch) {return (ch>32 && ch<128);}
template<class S> static std::vector<S> sort(std::vector<S> v)
{
	if(v.size()==0)
		return v;
	S min;
	int ind;
	S temp;
	for(int i=0;i<v.size();i++)
	{
		min=v[i];
		ind=i;
		for(int j=i;j<v.size();j++)
		{
			if(v[j]<min)
			{
				ind=j;
				min=v[j];
			}
		}
		if(i!=ind)
		{
			temp=v[i];
			v[i]=v[ind];
			v[ind]=temp;
		}
	}
	return v;
}
template<class T> static T* sort(T* taex,int l)
{
	if(l==0)
		return taex;
	T* ta=new T[l];
	for(int i=0;i<l;i++)
		ta[i]=taex[i];
	T min;
	int ind;
	T temp;
	for(int i=0;i<l;i++)
	{
		min=ta[i];
		ind=i;
		for(int j=i;j<l;j++)
		{
			if(ta[j]<min)
			{
				min=ta[j];
				ind=j;
			}
		}
		if(i!=ind)
		{
			temp=ta[i];
			ta[i]=ta[ind];
			ta[ind]=temp;
		}
	}
	return ta;
}
template<class T> static bool contains(std::vector<T> v,T el)
{
	for(T e : v)
	{
		if(e==el)
			return true;
	}
	return false;
}
template<class T> static int countUnique(std::vector<T> v)
{
	int c=0;
	bool flag;
	for(int i=0;i<v.size();i++)
	{
		flag=true;
		for(int j=i+1;j<v.size();j++)
		{
			if(v[i]==v[j])
			{
				//std::cout << i << j <<"\n";
				flag=false;
				break;
			}
		}
		if(flag)
			c++;
	}
	return c;
}
template<class T> static std::vector<T> set(std::vector<T> v)
{
	std::vector<T> ret;
	for(T el : v)
	{
		if(!contains(ret,el))
			ret.push_back(el);
	}
	return ret;
}
template<class T> static int find(std::vector<T> v,const T& o,int stI=0)
{
	for(int i=stI;i<v.size();i++)
	{
		if(v[i]==o)
			return i;
	}
	return -1;
}
//template<class T>
template<class T> static T* append(T* ta,T t,int l) //Generates append function to append a value to an array of ANY DATA TYPE
{
	T* ret=new T[l+1];
	for(int i=0;i<l;i++)
		ret[i]=ta[i];
	ret[l]=t;
	return ret;
}
template<class T> static T* merge(T* ta,T* tb,int a1,int a2=1)
{
	T* ret=new T[a1+a2];
	for(int i=0;i<a1;i++)
		ret[i]=ta[i];
	for(int i=0;i<a2;i++)
		ret[i+a1]=tb[i];
	return ret;
}
template<class T> static T* remove(T* ta,int el,int l) //Generates remove function to remove a value from any array BY INDEX
{
	T* ret=new T[l-1];
	for(int i=0;i<el;i++)
		ret[i]=ta[i];
	for(int i=el+1;i<l;i++)
		ret[i-1]=ta[i];
	return ret;
}
template<class T> static T& pop(std::vector<T>& v,int ind)
{
	T* temp=new T(v[ind]);
	v.erase(v.begin()+ind);
	return *temp;
}
template<class N> static N sum(std::vector<N> v)
{
	N ret=0;
	for(N no : v)
		ret=ret+no;
	return ret;
}
template<class N> static N sum(N* a,int n)
{
	N sum=0;
	for(int i=0;i<n;i++) sum+=a[i];
	return sum;
}
template<class N> static bool contains(const N* a,const N& o,int l)
{
	for(int i=0;i<l;i++)
	{
		if(a[i]==o)
			return true;
	}
	return false;
}
template<class N> static int indexOf(const N& o,const N* a,int l,int s=0)
{
	for(int i=s;i<l;i++)
	{
		if(a[i]==o)
			return i;
	}
	return -1;
}
template<class T,class U=T> static int indexOf(const T& o,const std::vector<U>& v,int l=-1,int s=0)
{
	if(l==-1) l=v.size();
	for(int i=s;i<l;i++)
	{
		if(v[i]==o)
			return i;
	}
	return -1;
}

namespace stringfx
{
	/**@brief Empty all remaining data from a stringstream into a string*/
	std::string drain(std::stringstream& ss)
	{
		std::string ret="",temp;
		while(true)
		{
			ss >> temp;
			ret+=" "+temp;
			if(ss.tellg()==-1) break;
		}
		return ret.substr(1,ret.size()-1);
	}
	std::string drain(std::istream& is)
	{
		std::string ret="",temp;
		while(true)
		{
			getline(is,temp);
			if(is.eof()) break;
			ret+=temp+"\n";
		}
		return ret;
	}
	/**@brief Find the index of a character in a string*/
	inline int indexOf(char c,const std::string& s) {return s.find(c);}
	/**@brief Find the last index of a character in a string*/
	inline int lastIndexOf(char c,const std::string& s) {return s.find_last_of(c);}
	/**@brief Find the index of a substring in a string*/
	inline int indexOf(const char* st,const std::string& s) {return s.find(st);}
	inline int indexOf(const std::string& st,const std::string& s) {return s.find(st);}
	/**@brief Find the index of any character from a list (given as a string) in another string. Priority is given to the earlier characters in the string.
	   @details Character priorities are decided by order in the list string. The earlier characters have higher priority.<br/>
		 The index of a character lower in priority will be returned <i>only</i> if no characters of higher priority exist in the string.
	*/
	inline int indexOfAnyByPriority(const std::string& list,const std::string& s)
	{
		int ind=-1,pri=list.length();
		for(int i=0;i<s.length();i++)
		{
			for(int j=0;j<list.length();j++)
			{
				if(s[i]==list[j])
				{
					if(!j) return i;
					if(j<pri) {pri=j; ind=i;}
					break;
				}
			}
		}
		return ind;
	}
	/**@brief Find the index of the last occurence of any character from a list (given as a string) in another string. Priority is given to the earlier characters in the string.
	   @details Character priorities are decided by order in the list string. The earlier characters have higher priority.<br/>
		 The index of a character lower in priority will be returned <i>only</i> if no characters of higher priority exist in the string.
	*/
	inline int lastIndexOfAnyByPriority(const std::string& list,const std::string& s)
	{
		int ind=-1,pri=list.length();
		for(int i=s.length()-1;i>=0;i--)
		{
			for(int j=0;j<list.length();j++)
			{
				if(s[i]==list[j])
				{
					if(!j) return i;
					if(j<pri) {pri=j; ind=i;}
					break;
				}
			}
		}
		return ind;
	}
	/**@brief split a string by a given delimiter and generate a std::vector of strings as output
		 @details This function also allows "protector characters" where a string between a start protector and end protector remains uncut
		<br/>For example ((1,2,3),(4,5,6)) when split over ',' will yield (1,2,3),(4,5,6) if '(' and ')' are used as protector charactes in the start and end positions respectively.
	*/
	std::vector<std::string> split(const std::string& s,char delim=' ',bool autotrim=false,char protB=(char)0,char protE=(char)0)
	{
		std::vector<std::string> ret;
		std::string temp="";
		int layer=0;
		for(int i=0;i<s.length();i++)
		{
			if(protB && s[i]==protB) layer++;
			if(protE && s[i]==protE) layer--;
			if(s[i]==delim && !layer)
			{
				if(!temp.length() && autotrim) continue;
				ret.push_back(temp); temp=""; continue;
			}
			temp+=s[i];
		}
		ret.push_back(temp);
		return ret;
	}
	/**@brief Replace a substring with another
	   @param[in] s: Original string
		 @param[in] val: Substring to replace
		 @param[in] ns: New string to replace with
	*/
	inline std::string& replace(std::string& s,std::string val,std::string ns)
	{
		int l=s.find(val);
		if(l==std::string::npos) return s;
		return s.replace(l,val.size(),ns);
	}
	/**@brief Trim (remove trailing spaces from) a string*/
	inline std::string trim(const std::string& aString) {if(aString.find_first_not_of(' ')==std::string::npos) return ""; return aString.substr(aString.find_first_not_of(' '), (aString.find_last_not_of(' ') - aString.find_first_not_of(' ')) + 1);}

	/**@brief Matches an identifier to an atom name. Used when loading rules and categories
	   @details This function also uses information of wild characters (like '*') when matching two strings.<br/>
		 Returns 'true' if there is a match and 'false' if not.
	*/
	bool matches(const std::string& seq,const std::string& str)
	{
		if(seq==str) return true;
		bool wild=false,pos=false;
		std::string cmp=seq;
		if(seq[seq.size()-1]=='*') {wild=true; pos=0; cmp=seq.substr(0,seq.size()-1);}
		if(seq[0]=='*') {wild=true; pos=1; cmp=seq.substr(1,seq.size()-1);}
		if(!wild) return false;
		if(!pos) return (indexOf(cmp,str)==0);
		else return (indexOf(cmp,str)==str.size()-cmp.size());
	}
}

namespace commons
{
	/*
	 * Basic functions including
	 * 1) abs
	 */
	static bool pointersAreArrays=true;
	int abs(int i) {return ((i<0)?-i:i);}
	//int bitSplice(long int i,)
	//Exception classes
	class IndexOutOfBoundsException : public std::exception
	{
		std::string msg;
	public:
		IndexOutOfBoundsException()
		{
			msg="Called for array with out of bounds index.";
		}
		IndexOutOfBoundsException(int i)
		{
			msg=(std::string("Called for array with out of bounds index: ")<<i);
		}
		IndexOutOfBoundsException(int i,int l)
		{
			msg=(std::string("Called array of length=")<<l<<std::string(" with out of bounds index: ")<<i);
		}
		virtual const char* what() throw() {return msg.c_str();}
	};
	class InvalidArrayIndexException : public std::exception
	{
		std::string msg;
	public:
		InvalidArrayIndexException()
		{
			msg="Called for array with bad index.";
		}
		InvalidArrayIndexException(int i)
		{
			msg=(std::string("Called for array with bad index: ")<<i);
		}
		InvalidArrayIndexException(int i,int l)
		{
			msg=(std::string("Called array of length=")<<l<<std::string(" with bad index: ")<<i);
		}
		virtual const char* what() throw() {return msg.c_str();}
	};
	class InvalidFileModeException : public std::exception
	{
	public:
		InvalidFileModeException() {std::cerr << "Sorry. File mode not set for requested operation"<<"\n";}
		InvalidFileModeException(char md)
		{
			std::cerr << "File mode '"<<md<<"' not suitable for requested operation" <<"\n";
		}
		InvalidFileModeException(char md,char req)
		{
			std::cerr << "File mode '"<<md<<"' not suitable for requested operation. Mode='"<<req <<"'\n";
		}
	};
	class NoSuchKeyException : public std::exception
	{
	public:
		NoSuchKeyException() {}
		virtual const char* what() throw() {return "Key not found in hashmap";}
	};
	class InvalidArgumentException : public std::exception
	{
	public:
		InvalidArgumentException() {}
		virtual const char* what() throw() {return "Invalid argument given to a function";}
	};
	long int readBytes(std::istream& is,int b)
	{
	  unsigned long int n=0;
	  char c;
	  for(int i=0;i<b;i++)
	  {
	    n<<=8;
	    //n*=SHIFT;
	    if(!is.get(c)) return 0;
	    if(((int)c)<0) n+=c+256;
	    else n+=c;
	    //cout << (int)c << "\n";
	    //n+=;
	  }
	  return n;
	}
	//Helpful structs
	union char_ptr
	{
		char* p;
		const char* ptr;
	};
	//Package classes
	class String;
	class X
	{

	};
	class String
	{
		std::string str="";
		int length=0;
		int i=0;
		public:
			//static const String* EMPTY=new String(new std::string(""));
			String() {str="";}
			String(char* s) : String(std::string(s)) {}
			String(char ch)
			{
				str="";
				str+=ch;
				length=1;
			}
			String(int n)
			{
				//std::cout << "Called int construct\n";
				length=n;
				str=std::string(" ");
				for(int i=1;i<n;i++)
					str+=" ";
			}
			String(std::string s) {str=s;length=str.length();}
			String(char* s,int l)
			{
				length=l;
				str="";
				for(int i=0;i<l;i++)
					str+=s[i];
				/*length=std::string(s).length();
				length=std::min(l,length);*/
			}
			String(const String& st)
			{
				length=st.length;
				str=st.str;
			}
			String(Object o);
			int indexOf(const char& ch) const
			{
				for(int i=0;i<length;i++)
				{
					if(ch==str[i])
						return i;
				}
				return -1;
			}
			int size() {return length;}
			int indexOf(const String& cmp) const
			{
				//std::cout << cmp.toString() << "\n";
				bool s=false;
				for(int i=0;i<length;i++)
				{
					if(i+cmp.getLength()>length)
						return -1;
					s=true;
					for(int j=i;j<i+cmp.getLength();j++)
					{
						//std::cout << str[i] << " "<< cmp[j-i]<<"\n";
						if(str[j]!=cmp[j-i])
						{
							s=false;
							break;
						}
						/*else
						 * 	continue;*/
					}
					if(s)
						return i;
				}
				return -1;
			}
			String substring(int begin,int end) const
			{
				/*if(begin<0)
					begin=length+begin;
				if(begin<0)
					throw IndexOutOfBoundsException(begin,length);*/
				if(begin<0 ||  begin>length)
					throw IndexOutOfBoundsException(begin,length);
				std::string ret="";
				for(int i=begin;i<end;i++)
					ret+=str[i];
				return String(ret);
			}
			inline String substring(int begin) const {return substring(begin,length);}
			inline char at(int i) const {return charAt(i);}
			bool startsWith(const String& string) const
			{
				if(string.getLength()>getLength())
					return false;
				for(int i=0;i<string.getLength();i++)
				{
					if(str[i]!=string.str[i])
						return false;
				}
				return true;
			}
			bool endsWith(const String& string) const
			{
				if(string.getLength()>getLength())
					return false;
				int K=string.getLength()-1;
				for(int i=length;i>=0;i--)
				{
					if(str[i]!=string.str[K--])
						return false;
				}
				return true;
			}
			String trim() const
			{
				if(length<1)
					return *this;
				int st=0,en=length-1;
				while(st<length && !isPrintable(str[st]))
					st++;
				while(en>st && !isPrintable(str[en]))
					en--;
				if(st>en)
					return "";
				return substring(st,en+1);
			}
			String replace(const String& str,const String& pat) const
			{
				String ret(*this);
				int ind;
				while((ind=ret.indexOf(str))!=-1)
					ret=ret.substring(0,ind)+pat+ret.substring(ind+str.getLength());
				return ret;
			}
			// Operator overloading
			char& operator[](int i) //Works like in python
			{
				if(i<0)
				{
					if(i<-length)
						throw IndexOutOfBoundsException(i,length);
					return str[length+i];
				}
				if(i>=length)
					throw IndexOutOfBoundsException(i,length);
				return str[i];
			}
			const char& operator[](int i) const //Works like in python
			{
				if(i<0)
				{
					if(i<-length)
						throw IndexOutOfBoundsException(i,length);
					return str[length+i];
				}
				if(i>=length)
					throw IndexOutOfBoundsException(i,length);
				return str[i];
			}
			bool operator==(const String& st) const
			{
				/*if(st.length!=length)
					return false;*/
				return str==st.str;
				/*for(int i=0;i<length;i++)
				{
					if(st.str[i]!=str[i])
						return false;
					if(str[i]=='\0')
						return true;
				}
				if(st.length==length)
					return true;
				else
					return (st.length>length && st[length]=='\0');*/
			}
			inline bool operator!=(const String& st)const {return !operator==(st);}
			//bool operator==(std::string st) {return strcmp(str,st.c_str())==0;}
			//inline bool operator!=(std::string st) {return !operator==(st);}
			bool operator==(char ch) const
			{
				if(length!=1)
					return false;
				return str[0]==ch;
			}
			inline bool operator!=(char ch) const {return !operator==(ch);}
			inline bool operator>(const String& st) const {return str>st.str;}
			inline bool operator>=(const String& st) const {return str>=st.str;}
			inline bool operator<(const String& st) const {return str<st.str;}
			inline bool operator<=(const String& st) const {return str<=st.str;}
			String operator+(const String& st) const
			{
				String ret(st.length+length);
				for(int i=0;i<length;i++)
					ret[i]=str[i];
				for(int i=length;i<length+st.length;i++)
					ret[i]=st[i-length];
				return ret;
			}
			String operator+(const char& ch) const
			{
				String ret(1+length);
				for(int i=0;i<length;i++)
					ret[i]=str[i];
				ret[length]=ch;
				return ret;
			}
			String operator*(int no) const
			{
				String ret="";
				for(int i=0;i<no;i++)
					ret=ret+*this;
				return ret;
			}
			/*String& operator=(const String& ext)
			{
				length=ext.length;
				str=new char[length];
				for(int i=0;i<length;i++)
					str[i]=ext.str[i];
				return *this;
			}*/
			//Casting
			operator std::string() {return str;}
			double toDouble() {return std::stod(toString());}

			//Simple data functions
			inline char charAt(int i) const {return str[i];}
			inline int getLength() const {return length;}
			char* toCharArray() const
			{
				char* ret=new char[length+1];
				ret[length]='\0';
				for(int i=0;i<length;i++)
					ret[i]=str[i];
				return ret;
			}
			inline std::string toString() const {return str;}

			int countChar(const char& ch) const
			{
				int c=0;
				for(int i=0;i<length;i++)
				{
					if(str[i]==ch)
						c++;
				}
				return c;
			}
			String reverse() const
			{
				String ret(length);
				for(int i=0;i<length;i++)
					ret.str[i]=str[length-1-i];
				return ret;
			}
			bool isEmpty() {return (trim().getLength()==0);}
			//static methods
			static String* cut(String str,char delim)
			{
				int len=str.countChar(delim)+1;
				String* ret=new String[len];
				int K=0,ind;
				while((ind=str.indexOf(delim))!=-1)
				{
					ret[K++]=str.substring(0,ind);
					str=str.substring(ind+1);
				}
				ret[K++]=str;
				return ret;
			}
			static inline String cut(String str,char delim,int fld) {return (cut(str,delim)[fld-1]);}
	};
	static std::istream& operator >>(std::istream& is,String& str);
	static String getString(std::string in);

	template<class T> class CircularArray
	{
		//Class denoted by 'T' is expected to have a functional assignment operator, i.e. operator=(const T& ext) {} in place
	protected:
		int size=0;
		int ori=0;
		T* array;
	private:
		int pos=0;

	public:
		CircularArray() : CircularArray(0) {}
		CircularArray(int s) {size=s;array=new T[s];}
		CircularArray(T* ptr,int s)// : CircularArray(s)
		{
			size=s;array=new T[s];
			for(int i=0;i<size;i++)
				array[i]=ptr[i];
		}
		CircularArray(const CircularArray<T>& arr) : CircularArray(arr.size)
		{
			for(int i=0;i<size;i++)
				array[i]=arr.array[i];
		}
		~CircularArray() {delete[] array;}


		void setOrigin(int i){ori+=i; ori%=size;}
		const T* begin() const {return &array[ori];}
		T* begin() {return &array[ori];}
		const T* end() const {return &operator[](ori-1);}
		T* end() {return &operator[](ori);}
		int getSize() const {return size;}
		bool contains(const T& el) const
		{
			for(int i=0;i<size;i++)
			{
				if(operator[](i)==el)
					return true;
			}
			return false;
		}
		T* toLinearArray()
		{
			T* ret=new T[size];
			for(int i=ori;i<ori+size;i++)
				ret[i-ori]=operator[](i);
			return ret;
		}
		T& push_back(T t)
		{
			size++;
			T* next=new T[size];
			for(int i=0;i<size-1;i++)
				next[i]=array[i];
			next[size-1]=t;
			delete[] array;
			array=next;
			return array[size-1];
		}

		void setAll(const T& t)
		{
			for(int i=ori;i+ori+size;i++)
				operator[](i)=t;
		}

		bool getTraversal(const T& e1,const T& e2)
		{
			int i1=-1,i2=-1;
			for(int i=0;i<size;i++)
			{
				if(array[i]==e1)
					i1=i;
				if(array[i]==e2)
					i2=i;
				if(i1!=-1 && i2!=-1)
					break;
			}
			if(i1==i2)
				throw InvalidArgumentException();
			if(abs(i2-i1)>size/2)
			{
				if(i1<i2)
					i1+=size;
				else
					i2+=size;
			}
			return (i1<i2);
		}
		std::vector<T> shortestPath(const T& e1,const T& e2) const
		{
			std::vector<T> ret;
			int i1=-1,i2=-1;
			for(int i=0;i<size;i++)
			{
				if(array[i]==e1)
					i1=i;
				if(array[i]==e2)
					i2=i;
				if(i1!=-1 && i2!=-1)
					break;
			}
			if(i1==i2)
			{
				ret.push_back(array[i1]);
				return ret;
			}
			if(abs(i2-i1)>size/2)
			{
				if(i1<i2)
					i1+=size;
				else
					i2+=size;
			}
			if(i1>i2)
			{
				for(int i=i1;i>=i2;i--)
					ret.push_back(operator[](i));
			}
			else
			{
				for(int i=i1;i<=i2;i++)
					ret.push_back(operator[](i));
			}
			return ret;
		}
		//operator overloading
		T& operator[](int i)
		{
			i%=size;
			if(i<0)
				i+=size;
			i+=ori;
			return array[i%size];
		}
		const T& operator[](int i) const
		{
			i%=size;
			if(i<0)
				i+=size;
			i+=ori;
			return array[i%size];
		}
		CircularArray<T>& operator=(const CircularArray<T>& ext)
		{
			delete[] array;
			array=new T[ext.size];
			for(int i=0;i<ext.size;i++)
				array[i]=ext.array[i];
			size=ext.size;
			ori=ext.ori;
			return *this;
		}

		//Array manipulations
		const T& operator*() const {return operator[](pos);}
		T& operator*() {return operator[](pos);}
		void operator++() {pos++;}

	};
	template<class T> static std::ostream& operator<<(std::ostream& os,const CircularArray<T>& ca)
	{
		for(int i=0;i<ca.getSize();i++)
			os << ca[i];
		return os;
	}
	class File
	{
		int lno=0,lcount=-1;
		std::ofstream file_out;
		std::ifstream file_in;
		String fnam;
		char mode='r';
	public:
		File() {}
		File(const String& str,char mod)
		{
			fnam=str;mode=mod;
			if(mode=='W')
				mode='w';
			setup();
		}
		//explicit File(std::string str,char mod) : File(*(new String(str)),mode) {}
		File(const File& other)
		{
			lno=0; lcount=other.lcount;
			fnam=other.fnam;
			mode=other.mode;
			setup();
		}
		~File() {this->close();}

	private:
		void setup()
		{
			if(mode=='w')
				file_out.open(fnam.toString());
			else
			{
				file_in.open(fnam.toString());
				if(!file_in && mode!='w')
					std::cerr << "WARN: file did not open: "<<fnam.toString()<<"\n";
			}
		}
	public:
		inline std::ostream& getOutputStream() {return file_out;}
		inline std::istream& getInputStream() {return file_in;}
		inline String getName() {return fnam;}
		void resetRead()
		{
			if(mode=='w')
				return;
			file_in.clear();                 // clear fail and eof bits
			file_in.seekg(0, std::ios::beg); //Go to beginning
		}
		void close()
		{
			if(mode=='w')
				file_out.close();
			else
				file_in.close();
		}
		int getLineCount()
		{
			std::string line;
			if(lcount!=-1)
				return lcount;
			lcount=0;
			resetRead();
			while(std::getline(file_in,line))
				lcount++;
			resetRead();
			for(int i=0;i<lno;i++)
				getline(file_in,line);
			return lcount;
		}
		String nextLine()
		{
			if(mode=='w')
				throw InvalidFileModeException(mode,'r');
			String in=*(new String());
			if(file_in >> in)
			{
				lno++;
				return in;
			}
			else
				return String((char)EOF);
		}
		String nextWord()
		{
			if(mode=='w')
				throw InvalidFileModeException(mode,'r');
			std::string in=*(new std::string(""));
			if(file_in >> in)
				return getString(in);
			else
				return *(new String((char)EOF));
		}
		void write(std::string str)
		{
			if(mode!='w')
				throw InvalidFileModeException(mode,'w');
			file_out << str;
		}
		void writeLine(std::string str)
		{
			write(str);
			file_out << "\n";
		}
		inline void write(String str) {write(str.toString());}
		inline void writeLine(String str) {writeLine(str.toString());}
		void resetFile() {setup();}
		//static functions
		static File openFile(String msg,char md)
		{
			String fn=*(new String());
			std::cout << msg.toString()<<": ";
			std::cin >> fn;
			std::cout << fn.toString() << "\n";
			return *(new File(fn,md));
		}
		inline static File openFileRead(String msg) {return openFile(msg,'r');}
		inline static File openFileRead() {return openFile("Enter filename: ",'r');}
		inline static File openFileWrite(String msg) {return openFile(msg,'w');}
		inline static File openFileWrite() {return openFile("Enter filename: ",'w');}
		inline static File openFile() {return openFileRead();}
	};
	//template<typename ...> class HashMap;
	template<typename K,typename V> class HashMap
	{
		std::vector<K> keys;
		std::vector<V> values;
	public:
		HashMap() {}

		void append(K key,V val)
		{
			if(::contains(keys,key))
				values[::find(keys,key)]=val;
			else
			{
				keys.push_back(key);
				values.push_back(val);
			}
		}
		int size() {return keys.size();}
		bool contains(const K& key) {return ::contains(keys,key);}
		V get(const K& key) {operator[](key);}
		const std::vector<K>& getKeys() const {return keys;}
		const std::vector<V>& getValues() const {return values;}
		std::vector<K>& getKeys()  {return keys;}
		std::vector<V>& getValues()  {return values;}
		bool containsKey(const K& key){return ::contains(keys,key);}
		void remove(K key)
		{
			int n=-1;
			for(int i=0;i<keys.size();i++)
			{
				if(keys[i]==key)
				{
					n=i;
					break;
				}
			}
			if(n==-1)
				return;
			keys.erase(keys.begin()+n);
			values.erase(values.begin()+n);
		}

		//operator overloading
		V& operator[](const K& key)
		{
			if(::contains(keys,key))
				return values[::find(keys,key)];
			throw NoSuchKeyException();
		}
	};
	template<class Sortable> class SortedList
	{
		std::vector<Sortable> elems=*(new std::vector<Sortable>());
	public:
		SortedList() {}
		SortedList(Sortable el) :SortedList() {elems.push_back(el);}
		SortedList(Sortable* el,int n)
		{
			for(int i=0;i<n;i++)
				elems.push_back(el[i]);
			elems=sort(elems);
		}
		SortedList(std::vector<Sortable> v) {elems=sort(v);}

		inline std::vector<Sortable>& getElements() {return elems;}
		void addElement(Sortable s)
		{
			elems.push_back(s);
			elems=sort(elems);
			//elems=bubbleSort(elems);
		}
		void removeDuplicates() {elems=set(elems);}
		int search(Sortable s)
		{
			int ul=elems.size()-1,ll=0;
			int m;
			while(ll<=ul)
			{
				m=(ul+ll)/2;
				if(elems[m]==s)
					return m;
				if(elems[m]<s)
					ll=m+1;
				else
					ul=m-1;
			}
			return -1;
		}

		//Operator overloading
		SortedList<Sortable> operator+(SortedList<Sortable> s2)
		{
			std::vector<Sortable> vec=elems;
			vec.insert(vec.end(),s2.elems.begin(),s2.elems.end());
			return *(new SortedList<Sortable>(vec));
			//SortedList ret=
		}
		SortedList<Sortable> operator+=(SortedList s2)
		{
			elems.insert(elems.end(),s2.elems.begin(),s2.elems.end());
			elems=sort(elems);
			return *this;
		}
	};
	/*class IndexOutOfBoundsException : public exception
	{
		public:
			IndexOutOfBoundsException()
			{
				std::cout << "Called for array with out of bounds index." <<"\n";
			}
			IndexOutOfBoundsException(int i)
			{
				std::cout <<"Called for array with out of bounds index: "<<i<<".\n";
			}
			IndexOutOfBoundsException(int i,int l)
			{
				std::cout << "Called array of length="<<l<<" with out of bounds index: "<<i<<".\n";
			}
	};*/
}
static commons::String commons::getString(std::string in)
{
	commons::String ret(in); //=*(new commons::String());
	return ret;
}
static std::ostream& operator<<(std::ostream& os,commons::String str)
{
	/*for(int i=0;i<str.getLength();i++)
		os << str[i];*/
	os << str.toString();
	return os;
}
static std::istream& commons::operator >>(std::istream& is,commons::String& str)
{
	//is.ignore();
	std::string in=*(new std::string(""));
	if(!std::getline(is,in))
	{
		str=*(new String((char)EOF));
		return is;
	}
	/*char* testStr=new char[in.length()];
	std::strcpy(testStr,in.c_str());
	str=*(new commons::String(testStr));*/
	str=getString(in);
	return is;
}

static commons::String format(const commons::String& str,int size,bool align) //align=true => left align,align=false=>right align, size is no of spaces to allocate
{
	commons::String ret="";
	if(size==str.getLength())
		return commons::String(str);
	if(size<str.getLength())
		return str.substring(0,size);
	ret=str;
	commons::String spaces=commons::String(" ")*(size-str.getLength());
	if(align)
		ret=ret+spaces;
	else
		ret=spaces+ret;
	return ret;
}
static commons::String format(double num,int acc,int size,bool align=false)
{
	std::ostringstream strs;
	strs.precision(acc);
	strs << num;
	std::string str = strs.str();
	return format(commons::String(str),size,align);
}

//Some template functions
template<class T> static std::ostream& operator<<(std::ostream& os,std::vector<T> v)
{
	os <<"{";
	for(T el : v)
		os << el<<",";
	os <<"}\n";
	return os;
}
template<class T,size_t N> constexpr size_t sizeOf(T(&)[N]) {return N;}
template<class T,size_t N> constexpr std::array<T,N+1> append(std::array<T,N> ta,T obj)
{
	std::array<T,N+1> ret;
	for(int i=0;i<N;i++)
		ret[i]=ta[i];
	ret[N]=obj;
	return ret;
}
template<class T> static inline commons::String toString(T obj)
{
	std::ostringstream strs;
	strs << obj;
	std::string str = strs.str();
	return commons::String(str);
}
template<class T> static inline commons::String getString(T obj) {return toString(obj);}
/*static std::istream& operator>>(std::istream& is,commons::String& st)
{
 //std::cout <<"Getline"<<"\n";
	std::getline(is,st.str);
	st.length=std::string(st.str).length();
	std::cout <<st.str<<"\n";
	return is;
}*/

//General functions that can later change
template<class T> T input(commons::String msg=commons::String("Input: "))
{
	std::cout << msg;
	T ret;
	std::cin >> ret;
	return ret;
}

class Object
{
public:
	Object() {}

	virtual commons::String toString() const {return commons::String("Object at ")+getString((long)this)+commons::String("\n");}
	//Object clone();

	bool operator==(const Object& o) {return this==&o;}
};
commons::String::String(Object o) : String(o.toString()) {}

static std::ostream& operator<<(std::ostream& os,const Object& o)
{
	os << o.toString();
	return os;
}
#endif
