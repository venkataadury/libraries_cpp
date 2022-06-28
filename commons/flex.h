#ifndef INCLUDED_FLEXIBILITY
#define INCLUDED_FLEXIBILITY 1
#include <algorithm>
#include <vector>

namespace flex
{
	class Cloneable
	{
	public:
		virtual void* clone() const=0;
	};
	template<typename Base, typename T> inline bool instanceof(const T*) { return std::is_base_of<Base, T>::value;}
	template<class T> class Iterable
	{
		T start,stop;
	public:
		explicit Iterable(T s,T e) {start=s; stop=e;;}
		explicit Iterable(std::vector<T> se) {start=se[0];stop=se[1];}


	public:
		virtual void next(T& i) =0;
	public:
		class iterator: public std::iterator<
		std::input_iterator_tag,   // iterator_category
		T,                      // value_type
		long,                      // difference_type
		const T*,               // pointer
		T                       // reference
		>{
			T current;
			Iterable<T>* obj;
		public:
			explicit iterator(T t,Iterable<T>* o) : obj(o) {current=t;}
			iterator& operator++() {obj->next(current); return *this;}
			iterator operator++(int) {iterator retval = *this; ++(*this); return retval;}
			bool operator==(iterator other) const {return current == other.current;}
			bool operator!=(iterator other) const {return !(*this == other);}
			#ifdef RANDOMACCESS
				virtual int operator-(iterator other) const =0; //{return (**this)-(*other);}
			#endif
			const T& operator*() const {return current;}
		};
		virtual iterator begin() {return iterator(start,this);}
		virtual iterator end() {return iterator(stop,this);}
	};
	class Range : public Iterable<long>
	{
		long START,STOP;
	public:
		void next(long& cur) override
		{
			if(START>=STOP)
				cur--;
			else
				cur++;
		}
	public:
		Range(long st,long en) : Iterable(st,en) {START=st; STOP=en;}
	};
	template<class T,class P> class IndexIterator : public std::iterator<
	std::input_iterator_tag,   // iterator_category
	T,                      // value_type
	long,                      // difference_type
	const T*,               // pointer
	T                       // reference
	>
	{
		P* parent;
		int ind;
	public:
		explicit IndexIterator(int i,P* par) : parent(par) {ind=i;}
		IndexIterator& operator++() {ind++; return *this;}
		IndexIterator operator++(int) {IndexIterator retval = *this; ++(*this); return retval;}
		bool operator==(IndexIterator other) const {return ind == other.ind;}
		bool operator!=(IndexIterator other) const {return !(*this == other);}
		T& operator*() {return parent->operator[](ind);}
		const T& operator*() const {return parent->operator[](ind);}
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
