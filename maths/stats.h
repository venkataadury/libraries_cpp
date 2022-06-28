#ifndef STATS_INCLUDED
#define STATS_INCLUDED 1
#include "commons/commons.h"
#include "maths/maths.h"
using namespace std;
using namespace maths;
typedef std::vector<double> Sample;

//Basic Template Functions
namespace randgen
{
	double throwarandompoint_normal(double mean=0,double stdev=1.0) //PCMODIFIED: ADDITION (added mean and stdev parameters)
	{
		random_device rd; mt19937 eng(rd());
		normal_distribution<double> distr(mean,stdev);
	  return distr(eng);
	}
	double throwarandompoint(double min=0, double max=1){
		random_device rd; mt19937 eng(rd());
		uniform_real_distribution<> distr(min,max);
		return distr(eng);
	}
	inline unsigned long randInt(unsigned long int L) {return throwarandompoint(0,L);}
	template<class T> const T& choice(const std::vector<T>& v) {return v[(int)throwarandompoint(0,v.size())];}
	template<class T> T& choice(std::vector<T>& v) {return v[(int)throwarandompoint(0,v.size())];}
	int weightSelect(double* A,int n)
	{
		double maxV=sum(A,n);
		double rn=throwarandompoint(0,1)*maxV;
		int i;
		for(i=0;i<n;i++)
		{
			rn-=A[i];
			if(rn<0) break;
		}
		return i;
	}
	int weightSelect(const std::vector<double>& wts)
	{
		double v[wts.size()]; for(int i=0;i<wts.size();i++) v[i]=wts[i];
		return weightSelect(v,wts.size());
	}
}
namespace stats
{
	struct StatData {double mean,stddev;};
	struct StatRange {double b=0,e=1;};

	class SampleSizeMismatchException : public std::exception
	{
	public:
		virtual const char* what() throw() {return "Sample sizes did not match (where necessary) for some stats operation requested.";}
	};

	static Sample randomSample(int mn,int mx,int ss=1000)
	{
		cout << "Random int sample\n";
		Sample s;
		for(int i=0;i<ss;i++)
			s.push_back(randint(mn,mx));
		return s;
	}
	static Sample randomSample(double mn,double mx,int ss=1000)
	{
		Sample s;
		for(int i=0;i<ss;i++) s.push_back(randgen::throwarandompoint(mn,mx));
		return s;
	}
	std::vector<double> toNumList(char** str,int n)
	{
		std::vector<double> ret;
		for(int i=0;i<n;i++)
			ret.push_back(std::stod(str[i]));
		return ret;
	}
	std::vector<double> toNumList(commons::String* str,int n)
	{
		std::vector<double> ret;
		for(int i=0;i<n;i++)
			ret.push_back(std::stod(str[i]));
		return ret;
	}
	double sum(std::vector<double> nums)
	{
		double s=0;
		for(double d : nums) s+=d;
		return s;
	}
	double prod(std::vector<double> nums)
	{
		double s=1;
		for(double d : nums) s*=d;
		return s;
	}
	double mean(std::vector<double> nums)
	{
		double s=0;
		for(double d : nums)
			s+=d;
		return s/((double)nums.size());
	}
	inline double avg(Sample s) {return mean(s);}
	double stddev(std::vector<double> nums)
	{
		double mn=mean(nums);
		double dev=0;
		for(double d : nums)
			dev+=sqr(d-mn);
		dev/=(nums.size()-1);
		return sqrt(dev);
	}
	double chi2(const Sample& sam,const Sample& templ)
	{
		std::vector<double> res;
		if(sam.size()!=templ.size()) throw SampleSizeMismatchException();
		for(int i=0;i<sam.size();i++) res.push_back(::sqr(sam[i]-templ[i])/templ[i]);
		return sum(res);
	}

	Sample bootstrapSample(Sample s)
	{
		Sample bs;
		for(int i=0;i<s.size();i++)
			bs.push_back(s[(int)randgen::throwarandompoint(0,s.size())]);
		return bs;
	}
	StatRange bootstrap(Sample s,double pc,int bsamp=0)
	{
		if(bsamp==0)
			bsamp=100*pc+2;
		Sample means,bssamp;
		double meanV=0;
		for(int i=0;i<bsamp;i++)
		{
			bssamp=bootstrapSample(s);
			meanV=mean(bssamp);
			means.push_back(meanV);
		}
		means=sort(means);
		double r=(100-pc)/2.0,mult=bsamp/100.0;
		StatRange ret;
		ret.b=means[(int)(r*mult)]; ret.e=means[(int)((100-r)*mult)-1];
		return ret;
	}
	std::vector<maths::Point2D> generatePoints(double ll,double ul,double (*f)(double i),int ss=250)
	{
		Sample xs=randomSample(ll,ul,ss);
		//cout << xs << "\n";
		std::vector<Point2D> ret;
		for(int i=0;i<ss;i++) ret.push_back(maths::Point2D((double)xs[i],(double)f(xs[i])));
		return ret;
	}
	template<class T> std::vector<maths::Point2D> generatePoints(double ll,double ul,double (T::*f)(double i) const,const T& obj,int ss=250,double errmean=0,double errstddev=0)
	{
		Sample xs=randomSample(ll,ul,ss);
		//cout << xs << "\n";
		std::vector<Point2D> ret;
		for(int i=0;i<ss;i++) ret.push_back(maths::Point2D(xs[i],(obj.*f)(xs[i])+randgen::throwarandompoint_normal(errmean,errstddev)));
		return ret;
	}
}
namespace gametools
{
	class Coin
	{
		double fairness=0.5;
	public:
		Coin(double fn=0.5) {fairness=fn;}

		std::vector<bool> toss(int n)
		{
			std::vector<bool> tosses;
			for(int i=0;i<n;i++)
				tosses.push_back(maths::toss(fairness));
			return tosses;
		}
		bool toss() {return maths::toss(fairness);}
	};

	class GeneralCard
	{
	protected:
		int ID=0;
		const int& getID() const {return ID;}
	public:
		GeneralCard() {}
	};
	typedef GeneralCard Card;

	class PokerCard : public GeneralCard
	{
		short int category=0; //0-hearts, 1-spades, 2-diamonds, 3-flowers
		//0 - Joker/undefined 1-10 -> A->10; 11-13 -> J-K
	public:
		PokerCard(short int cat=maths::randint(4)) : PokerCard(cat,maths::randint(13)+1) {}
		PokerCard(short int cat,int id=maths::randint(13)) : GeneralCard()
		{
			category=cat;
			ID=id;
		}

		const short int& getCategory() const {return category;}
	};


	class CardSet
	{
		std::vector<Card> cards;
	public:
		CardSet(std::vector<Card> cs) {cards=cs;}
		CardSet() {}

		std::vector<Card> getCards() {return cards;}
		int getCardCount() {return cards.size();}
		const Card& getCard(int i) const {return cards[i];}
		Card& getCard(int i) {return cards[i];}
	};
	class Pack
	{
	protected:
		CardSet deck;
		CardSet discard;

	public:
		Pack(int n=1);
	};
	class StandardDeck : public Pack
	{
		int decks=1;
	public:
		StandardDeck(int n=1) {decks=1; initDecks();}

	private:
		void initDecks()
		{

		}
	};
}
#endif
