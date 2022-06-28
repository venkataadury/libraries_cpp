#ifndef INCLUDED_FLEXIBILITY
#define INCLUDED_FLEXIBILITY 1
#include <algorithm>
#include <vector>
#include "commons/commons.h"

namespace flex
{
	template<class T> struct IteratorPos
	{
		T* el=nullptr;
		long int ind=0;
		IteratorPos() {el=nullptr; ind=0;}
		explicit IteratorPos(T& t,long int i=0) {el=&t; ind=i;}

		inline T& operator*() {return *el;}
		inline const T& operator*() const {return *el;}
		inline bool operator==(const IteratorPos<T>& i2) const {return (ind==i2.ind || el==i2.el);}
		inline bool operator!=(const IteratorPos<T>& i2) const {return !(ind==i2.ind || el==i2.el);}
	};
	template<class T> class Iterable
	{
	protected:
		IteratorPos<T> start,stop;
	public:
		explicit Iterable(IteratorPos<T> s,IteratorPos<T> e) {start=s; stop=e;}
		explicit Iterable(std::vector<IteratorPos<T>> se) {start=se[0];stop=se[1];}


	public:
		virtual void next(IteratorPos<T>& i) =0;
	public:
		class iterator: public std::iterator<
		std::input_iterator_tag,   // iterator_category
		IteratorPos<T>,                      // value_type
		long,                      // difference_type
		const T*,               // pointer
		T                       // reference
		>{
			IteratorPos<T> current;
			Iterable<T>* obj;
		public:
			explicit iterator(IteratorPos<T> t,Iterable<T>* o) : obj(o) {current=t;}
			iterator& operator++() {obj->next(current); return *this;}
			iterator operator++(int) {iterator retval = *this; ++(*this); return retval;}
			bool operator==(iterator other) const {return current == other.current;}
			bool operator!=(iterator other) const {return !(*this == other);}
			T& operator*() {return *current;}
		};
		virtual iterator begin() {return iterator(start,this);}
		virtual iterator end() {return iterator(stop,this);}
	};
	class Range : public Iterable<long>
	{
		long START,STOP;
	public:
		void next(IteratorPos<long>& cur) override
		{
			if(START>=STOP)
			{
				cur.ind--;
			}
			else
			{
				cur.ind++;
			}
		}
	public:
		Range(long st,long en) :START(st),STOP(en), Iterable(IteratorPos<long>(START,0),IteratorPos<long>(STOP,en-st)) {} //{std::cout << *start << "\t" << *stop << "\n"; std::cout << start.ind << "\t" << stop.ind << "\n";}
	};
	/*template <long START,long END> class Range
	{
		long num=START;
	public:
		typedef long iterator;
		typedef const long const_iterator;
		iterator begin() {return START;}
		iterator end() {return END>=START?END-1:END+1;}
		iterator& operator++() {if(END>=START) num=num+1; else num=num-1; return num;}
		bool operator==(iterator other) const {return num == other;}
		bool operator!=(iterator other) const {return !(num == other);}
		long operator*() const {return num;}
	};*/
}
#endif
