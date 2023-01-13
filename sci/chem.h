#ifndef INCLUDED_CHEM
#define INCLUDED_CHEM 1
#include "maths/maths.h"
#include <array>
#include <assert.h>
#include "commons/flex.h"
#include "commons/commons.h"
#include "commons/strfx.h"

/*
 * PLEASE READ USAGE especially for "Molecule" class
 */
using namespace commons;
namespace chemdata
{
	double getBondLength(String atom);
}
namespace molecules
{
	//Exceptions
	class MoleculeExceptions : public exception {};
	class NoAtomSelectedException : public MoleculeExceptions {virtual const char* what() {return "Trying to operate on molecule without selecting target atom.";}};
	class TooManyBondsException : public MoleculeExceptions {virtual const char* what() {return "Too many bonds for atom. Increase max bond count or check code.";}};
	class AtomNotinBondException : public MoleculeExceptions {virtual const char* what() {return "Requested atom is not paticipating in bond.";}};
	class NoDihedralAngleException : public AtomNotinBondException {virtual const char* what() {return "Requested atom(s) is/are not paticipating in dihedral angle.";}};
	class AtomNotInMoleculeException : public MoleculeExceptions {virtual const char* what() {return "Requested atom(s) is/are not the (same) molecule.";}};

	class Atom;
	class Molecule;

	//Methods declared at end of namespace
	static bool bonded(Atom* a1,Atom* a2);
	//typedef flex::Iterable<Atom*> Iterable<Atom*>;
	class Bond
	{
	public:
		Atom *a1,*a2;
		int order=1;
		bool aromatic=false,resonating=false;
		int gpflag=false,gpflag2=false; //general purpose flags (gps)
		int EZ=0; //(1) -> E ; (-1)->Z; 0->Unspecified

		Bond() {a1=a2=nullptr;}
		Bond(Atom* at1,Atom* at2)
		{
			a1=at1;a2=at2;
			if(a1==nullptr || a2==nullptr)
				cout << "Warn: Null atom in bond\n";
		}

		Atom& getSecond(Atom* a)
		{
			if(a1==a) return *a2;
			if(a2==a) return *a1;
			throw AtomNotinBondException();
		}
		Atom& getAtom1() {return *a1;}
		Atom& getAtom2() {return *a2;}
		void setEZ(int ez) {EZ=ez;}
		void setOrder(int o) {order=o;}
		int getOrder() {return order;}
		inline double getLength() const;
		inline Vector toVector() const;
		Vector toVector(const Atom* r) const;
		bool hasAtom(const Atom* ptr) {return (a1==ptr || a2==ptr);}
		void setAromatic(bool ar=true) {aromatic=ar;}
		bool isAromatic() {return aromatic;}
	};
	static Bond** getOrganicChain(int l);
	class Atom
	{
		Bond** bonds;
		int maxbonds;
		Molecule* myMolecule=nullptr;
		Vector position; //relative position from main atom of molecule
		int electrons=0;
		String symbol;
		int K=0;
		int charge=0;
		bool chirality=0,ischiral=false; // 0 is counter-clockwise, 1 is clockwise
	public:
		bool fillFlag=false,aromFlag=false;
		int gpflag=0;
		char chain;
	public:
		Atom() : Atom("H") {}
		Atom(char* str) : Atom(String(str)) {}
		Atom(String el,int bcount=1);
		explicit Atom(String el,int bcount,unsigned int elecount) : symbol(el)
		{
			position=Vector(0,0,0);
			electrons=(elecount==0)?bcount:elecount;
			try {bonds=new Bond*[bcount];} catch(std::bad_alloc ex) {cout<<el<<"\n"; throw ex;}
			for(int i=0;i<bcount;i++)
				bonds[i]=nullptr;
			maxbonds=bcount;
		}
		Atom(String atomdata,Atom* bondto,unsigned int elecount) : symbol(String::cut(atomdata,' ')[0])
		{
			//cout << "Enter\n";
			String* dp=String::cut(atomdata,' ');
			//cout  << "Process\n";
			double x=dp[1].toDouble(),y=dp[2].toDouble(),z=dp[3].toDouble();
			//cout << "Done\n";
			int nb=(int)dp[4].toDouble();
			//cout << x << " "<< y<< " " << z<<"\n";
			position=Vector(x,y,z);
			if(elecount=0)
				electrons=nb;
			else
				electrons=elecount;
			bonds=new Bond*[nb];
			for(int i=0;i<nb;i++)
				bonds[i]=nullptr;
			maxbonds=nb;
		}
		Atom(const Atom& atm)
		{
			maxbonds=atm.maxbonds;
			bonds=new Bond*[maxbonds];
			for(int i=0;i<maxbonds;i++)
				bonds[i]=nullptr;
			myMolecule=nullptr;
			position=Vector(0,0,0);
			electrons=atm.electrons;
			symbol=atm.symbol;
			fillFlag=atm.fillFlag;
			K=0;
		}
		~Atom()
		{
			delete[] bonds;
		}
		/*~Atom()
		{
			for(int i=0;i<maxbonds;i++)
			{
				if(bonds[i]==nullptr)
					continue;
				delete &bonds[i]->getSecond(this);
			}
		}*/


		//Positioning
		void setExactLocation(geom3D::Point3D p3); // Look at End
		void setPosition(geom3D::Point3D pt) {position=pt.toVector();}
		void setRelativePosition(Vector v,const Atom& a) {position=a.position+v;}
		Vector getLocation() const {return position;}
		geom3D::Point3D getExactLocation() const; //Look at end
		Vector getRelativePosition(const Atom& at)
		{
			if(at.getMolecule()==getMolecule())
				return (position-at.position);
			return (getExactLocation().toVector()-at.getExactLocation().toVector());
		}

		//Molecule related
		inline Molecule* getMolecule() const {return myMolecule;}
		inline void setMolecule(Molecule& m) {myMolecule=&m;}
		inline void setMolecule(Molecule* mptr) {myMolecule=mptr;}
		void recursiveSetMolecule(Molecule* mptr)
		{
			if(myMolecule!=mptr)
			{
				myMolecule=mptr;
				for(int i=0;i<maxbonds;i++)
				{
					if(bonds[i]==nullptr)
						continue;
					else
						bonds[i]->getSecond(this).recursiveSetMolecule(mptr);
				}
			}
		}
		inline int getMaxBondCount() const {return maxbonds;}
		void setMaxBondCount(int bc)
		{
			assert(bc>=0);
			Bond** rbonds=new Bond*[bc];
			for(int i=0;i<min(bc,maxbonds);i++)
				rbonds[i]=bonds[i];
			for(int i=min(bc,maxbonds);i<bc;i++)
				rbonds[i]=nullptr;
			maxbonds=bc;
			if(electrons<maxbonds)
				electrons=maxbonds;
			bonds=rbonds;
		}
		void unrestrain() {setMaxBondCount(6);}
		void addBond(Bond* bond)
		{
			if((getBondedElectronCount()+bond->getOrder())>maxbonds)
			{
				if(getBondedElectronCount()+bond->getOrder()<=electrons) // To be tested properly
				{
					setMaxBondCount(maxbonds+1);
					charge++;
					bonds[K++]=bond;
					return;
				}
				else
					throw TooManyBondsException();
			}
			bonds[K++]=bond;
		}
		void removeBond(int n)
		{
			for(int i=n+1;i<maxbonds;i++)
				bonds[i-1]=bonds[i];
			bonds[maxbonds-1]=nullptr;
			K--;
		}
		//More data
		inline void addCharge(int n)
		{
			charge+=n;electrons-=n;
			setMaxBondCount(maxbonds+maths::abs(n));
			//cout << this->getBondCount() << "\t"<<this->getElectronCount() << "\t"<<this->getBondedElectronCount() << ": "<<(long)this<<"\n";
		}
		inline int getCharge() {return charge;}
		inline void flag(bool f=true) {fillFlag=f;}
		inline bool getFlag() const {return fillFlag;}
		inline bool setChirality(bool c) {chirality=c; ischiral=true;}
		inline bool isChiral() const {return ischiral;}
		inline bool getChirality() const //WARN doesn't check if it is chiral
		{return chirality;}
		inline char getChain() {return chain;}
		inline void flipChirality() {ischiral=true; chirality=!chirality;}
		inline int getElectronCount() const {return electrons;}
		bool isAromatic() const {return aromFlag;}
		void setAromatic(bool v=true) {aromFlag=v;}
		//Working with bonds
		Bond* getBond(int ind)
		{
			if(ind>=maxbonds)
				throw IndexOutOfBoundsException(ind,maxbonds);
			return bonds[ind];
		}
		Bond* getFirstBond()
		{
			if(maxbonds==0 || bonds[0]==nullptr)
				throw IndexOutOfBoundsException(1,0);
			return bonds[0];
		}
		Bond* getLastBond()
		{
			if(maxbonds==0 || bonds[0]==nullptr)
				throw IndexOutOfBoundsException(1,0);
			for(int i=maxbonds-1;i>=0;i--)
			{
				if(bonds[i]==nullptr)
					continue;
				return bonds[i];
			}
			return nullptr;
		}
		Bond* nextBond(int X)
		{
			if(X>=maxbonds) return nullptr;
			for(int i=X;i<maxbonds;i++)
			{
				if(bonds[i]==nullptr)
					continue;
				return bonds[i];
			}
			return nullptr;
		}
		int findBond(Atom* at)
		{
			for(int i=0;i<maxbonds;i++)
			{
				if(bonds[i]=nullptr)
					continue;
				if(&(bonds[i]->getSecond(this))==at)
					return i;
			}
			return -1;
		}
		Bond* getBondWith(Atom* atm)
		{
			for(int i=0;i<maxbonds;i++)
			{
				if(bonds[i]==nullptr)
					continue;
				if((&bonds[i]->getSecond(this))==atm)
					return bonds[i];
			}
			return nullptr;
		}
		bool hasFreeBond() {return getBondedElectronCount()<maxbonds;}
		Bond* bondTo(Atom* a2,int order=1)
		{
			Bond* bnd=new Bond(this,a2);
			bnd->setOrder(order);
			this->addBond(bnd); a2->addBond(bnd);
			return bnd;
		}

		inline String toString() const {return symbol;}
		int getBondedElectronCount()
		{
			double s=0;
			for(int i=0;i<maxbonds;i++)
			{
				if(bonds[i]==nullptr)
					continue;
				s+=bonds[i]->getOrder();
				/*if(bonds[i]->isAromatic())
					s+=0.5;*/
			}
			return (int)s;
		}
		int getFreeElectronCount() {return electrons-charge-getBondedElectronCount();}
		int getBondCount()
		{
			int c=0;
			for(int i=0;i<maxbonds;i++)
			{
				if(bonds[i]==nullptr)
					continue;
				c++;
			}
			return c;
		}
		void organize(Bond* b,double le);
		//operator overloading
		std::array<Atom,2> operator|(Atom a2)
		{
			std::array<Atom,2> ret=*(new std::array<Atom,2>());
			ret[0]=*this;
			ret[1]=a2;
			return ret;
		}
		bool operator==(Atom a2) {return this==&a2;}
		Atom& operator=(const Atom& atm)
		{
			maxbonds=atm.maxbonds;
			bonds=new Bond*[maxbonds];
			for(int i=0;i<maxbonds;i++)
				bonds[i]=nullptr;
			myMolecule=nullptr;
			position=Vector(0,0,0);
			electrons=atm.electrons;
			symbol=atm.symbol;
			K=0;
			return *this;
		}
		Molecule* generateMolecule();

		//Friend functions
		friend double Bond::getLength() const;
		friend Vector Bond::toVector() const;
		friend Vector Bond::toVector(const Atom* at) const;
	};
	/*** Missing functions from Bond ***/
	Vector Bond::toVector(const Atom* resp) const
	{
		if(resp==a1) return (a2->position-a1->position);
		if(resp==a2) return (a1->position-a2->position);
		throw AtomNotinBondException();
	}
	inline Vector Bond::toVector() const {return Vector(a2->position-a1->position);}
	inline double Bond::getLength() const {return toVector().getMagnitude();}

	//Data classes
	class BondLengthData
	{
		String sym1,sym2; //atoms (eg. H,H or C,H or C,C)
		double val; //In angstroms (Bond length)
		int order=1;
	public:
		BondLengthData() {sym1=sym2="H";val=0.74;}
		BondLengthData(String s1,String s2,double len)
		{
			sym1=s1; sym2=s2;
			val=len;
		}
		inline double value() {return val;}
		inline int getOrder() {return order;}
		inline void setOrder(int o) {order=o;}
		inline bool matches(Atom* atoms,int ord=1) {return (order==ord && operator[](atoms));}
		inline bool matches(std::array<Atom,2> atoms,int ord=1) {return (order==ord && operator[](atoms));}
		//operator overloading
		bool operator[](Atom* atoms) //Assume atoms to have two elements
		{
			//cout << sym1 << " " << sym2 << "\t" << atoms[0].toString()<< " "<< atoms[1].toString() << "\n";
			if(sym1==atoms[0].toString())
				return (sym2==atoms[1].toString());
			if(sym2==atoms[0].toString())
				return (sym1==atoms[1].toString());
			return false;
		}
		inline bool operator[](std::array<Atom,2> atoms) {return operator[](atoms.data());}
	};
	class LengthTolerenceData
	{
		String sym1,sym2; //atoms (eg. H,H or C,H or C,C)
		double val; //In angstroms (Bond length)
		int order=1;
	public:
		LengthTolerenceData() {sym1=sym2="H";val=0.03;}
		LengthTolerenceData(String s1,String s2,double len)
		{
			sym1=s1; sym2=s2;
			val=len;
		}
		inline double value() {return val;}
		inline int getOrder() {return order;}
		inline void setOrder(int o) {order=o;}
		inline bool matches(Atom* atoms,int ord=1) {return (order==ord && operator[](atoms));}
		inline bool matches(std::array<Atom,2> atoms,int ord=1) {return (order==ord && operator[](atoms));}
		//operator overloading
		bool operator[](Atom* atoms) //Assume atoms to have two elements
		{
			//cout << sym1 << " " << sym2 << "\t" << atoms[0].toString()<< " "<< atoms[1].toString() << "\n";
			if(sym1==atoms[0].toString())
				return (sym2==atoms[1].toString());
			if(sym2==atoms[0].toString())
				return (sym1==atoms[1].toString());
			return false;
		}
		inline bool operator[](std::array<Atom,2> atoms) {return operator[](atoms.data());}
	};

	class BondAngleData
	{
		String sym1,sym2,sym3; //atoms (eg. H,O,H or O,C,H or C,C)
		int bo1=1,bo2=1;
		double val; //In maths units (see maths)
	public:
		BondAngleData() {sym1=sym3="H";sym2="O";val=toRadians(104);}
		BondAngleData(String s1,String s2,String s3,double ang)
		{
			sym1=s1; sym2=s2; sym3=s3;
			val=ang;
		}
		BondAngleData(String s1,String s2,double ang)
		{
			sym1=s1; sym2=s2; sym3=s1;
			val=ang;
		}

		double value()
		{
			if(isDegrees())
				return toDegrees(val);
			return val;
		}
		void setBondOrder1(int v) {bo1=v;}
		void setBondOrder2(int v) {bo2=v;}
		bool operator[](Atom* atoms)
		{
			if(sym1==atoms[0].toString())
				return (sym2==atoms[1].toString() && sym3==atoms[2].toString());
			if(sym3==atoms[0].toString())
				return (sym2==atoms[1].toString() && sym1==atoms[2].toString());
			return false;
		}
		bool matches(Atom* atoms,int o1=1,int o2=1) //Assume atoms to have three elements
		{
			if(sym1==atoms[0].toString())
				return (sym2==atoms[1].toString() && sym3==atoms[2].toString() && (bo1==o1 && bo2==o2));
			if(sym3==atoms[0].toString())
				return (sym2==atoms[1].toString() && sym1==atoms[2].toString() && (bo1==o2 && bo2==o1));
			return false;
		}
		inline bool operator[](std::array<Atom,3> atoms) {return operator[](atoms.data());}
	};
	class BondCountData
	{
		String atom;
		int count=1;
	public:
		BondCountData(String at,int c) {atom=at.trim();count=c;}
		BondCountData(Atom a,int c) : BondCountData(a.toString(),c) {}

		bool matches(const Atom& at) const {return (at.toString()==atom);}
		int get() const {return count;}
		//operator overloading
		inline bool operator[](const Atom& at) {return matches(at);}
	};
	class ElectronCountData
	{
		String atom;
		int count=1;
	public:
		ElectronCountData(String at,int c) {atom=at.trim();count=c;}
		ElectronCountData(Atom a,int c) : ElectronCountData(a.toString(),c) {}

		bool matches(const Atom& at) const {return (at.toString()==atom);}
		int get() const {return count;}

		//operator overloading
		inline bool operator[](const Atom& at) {return matches(at);}
	};
	class Data
	{
		static std::vector<BondLengthData> bond_lengths;
		static std::vector<LengthTolerenceData> bond_tolerences;
		static std::vector<BondAngleData> bond_angles;
		static std::vector<BondCountData> bond_counts;
		static std::vector<ElectronCountData> elec_counts;
	public:
		static constexpr double TETRAHEDRAL_ANGLE=109.4712206;
		static void addLength(Atom s1,Atom s2,double len,int ord=1)
		{
			BondLengthData added=*(new BondLengthData(s1.toString(),s2.toString(),len));
			added.setOrder(ord);
			bond_lengths.push_back(added);
		}
		static void addTolerence(Atom s1,Atom s2,double tol,int ord=1)
		{
			LengthTolerenceData added=*(new LengthTolerenceData(s1.toString(),s2.toString(),tol));
			added.setOrder(ord);
			bond_tolerences.push_back(added);
		}
		static void addAngle(Atom s1,Atom s2,Atom s3,double ang,int o1=1,int o2=1)
		{
			if(isDegrees())
				ang=toRadians(ang);
			BondAngleData nba=*(new BondAngleData(s1.toString(),s2.toString(),s3.toString(),ang));
			nba.setBondOrder1(o1);
			nba.setBondOrder2(o2);
			bond_angles.push_back(nba);
		}
		static void addBondCount(Atom a1,int n)
		{
			BondCountData bcd=*(new BondCountData(a1,n));
			bond_counts.push_back(bcd);
		}
		static void addElectronCount(Atom a,int n)
		{
			ElectronCountData ec=*(new ElectronCountData(a,n));
			elec_counts.push_back(ec);
		}
		//Bond length functions
		static bool hasLength(Atom a1,Atom a2,int ord=1) {return hasLength(a1|a2,ord);}
		static bool hasLength(Atom* ats,int ord=1) //assume two atoms in array
		{
			for(BondLengthData bld : bond_lengths)
			{
				if(bld[ats] && bld.getOrder()==ord)
					return true;
			}
			return false;
		}
		static bool hasLength(std::array<Atom,2> ats,int ord=1) {return hasLength(ats.data(),ord);}

		static double getLength(Atom* ats,int ord=1) //assume two atoms in array
		{
			for(BondLengthData bld : bond_lengths)
			{
				if(bld.matches(ats,ord))
					return bld.value();
			}
			return -1;
		}
		static inline double getLength(std::array<Atom,2> ats,int ord=1) {return getLength(ats.data(),ord);}
		inline static double getLength(Atom a1,Atom a2,int ord=1) {return getLength(a1|a2,ord);}

		//Bond length tolerence functions
		static bool hasLengthTolerence(Atom a1,Atom a2,int ord=1) {return hasLengthTolerence(a1|a2,ord);}
		static bool hasLengthTolerence(Atom* ats,int ord=1) //assume two atoms in array
		{
			for(LengthTolerenceData bld : bond_tolerences)
			{
				if(bld[ats] && bld.getOrder()==ord)
					return true;
			}
			return false;
		}
		static bool hasLengthTolerence(std::array<Atom,2> ats,int ord=1) {return hasLengthTolerence(ats.data(),ord);}

		static double getLengthTolerence(Atom* ats,int ord=1) //assume two atoms in array
		{
			for(LengthTolerenceData bld : bond_tolerences)
			{
				if(bld.matches(ats,ord))
					return bld.value();
			}
			return 0.026;
		}
		static inline double getLengthTolerence(std::array<Atom,2> ats,int ord=1) {return getLengthTolerence(ats.data(),ord);}
		inline static double getLengthTolerence(Atom a1,Atom a2,int ord=1) {return getLengthTolerence(a1|a2,ord);}


		//Bond angle functions
		inline static bool hasAngle(Atom a1,Atom a2,Atom a3,int o1=1,int o2=1) {return hasAngle(append((a1|a2),a3),o1,o2);}
		static bool hasAngle(Atom* ats,int o1=1,int o2=1) //assume three atoms in array
		{
			for(BondAngleData bad : bond_angles)
			{
				if(bad.matches(ats,o1,o2))
					return true;
			}
			return false;
		}
		static bool hasAngle(std::array<Atom,3> ats,int o1=1,int o2=1) {return hasAngle(ats.data(),o1,o2);}
		static double getAngle(Atom* ats,int o1=1,int o2=1) //assume three atoms in array
		{
			for(BondAngleData bad : bond_angles)
			{
				if(bad.matches(ats,o1,o2))
					return bad.value();
			}
			return -1;
		}
		inline static double getAngle(std::array<Atom,3> ats,int o1=1,int o2=1) {return getAngle(ats.data(),o1,o2);}
		inline static double getAngle(Atom a1,Atom a2,Atom a3,int o1=1,int o2=1) {return getAngle(append((a1|a2),a3),o1,o2);}

		//Bond count functions
		static bool hasBondCount(Atom a1)
		{
			for(BondCountData bcd : bond_counts)
			{
				if(bcd[a1])
				{
					//cout << bcd.atom << "\t"<<a1.toString() << "\n";
					return true;
				}
			}
			return false;
		}
		static int getBondCount(const Atom& a1)
		{
			for(BondCountData bcd : bond_counts )
			{
				if(bcd[a1])
					return bcd.get();
			}
			return -1;
		}

		//Electron count functions
		static bool hasElectronCount(Atom a1)
		{
			for(ElectronCountData ecd : elec_counts)
			{
				if(ecd[a1])
					return true;
			}
			return false;
		}
		static int getElectronCount(const Atom& a1)
		{
			for(ElectronCountData ecd : elec_counts )
			{
				if(ecd[a1])
					return ecd.get();
			}
			return -1;
		}
	};



	std::vector<BondLengthData> Data::bond_lengths;
	std::vector<LengthTolerenceData> Data::bond_tolerences;
	std::vector<BondAngleData> Data::bond_angles;
	std::vector<BondCountData> Data::bond_counts;
	std::vector<ElectronCountData> Data::elec_counts;


	class Molecule
	{
		/*
		 * Molecule class consists of a root atom, and then its bonds
		 * DO NOT assign a molecule from another. ONLY PASS AND USE POINTERS
		 * DO NOT use : Molecule m= *(new Molecule(Atom("H"))). Instead DO use: Molecule* m=new Molecule(Atom("H"))
		 */
		Atom& root;
		Atom* current=nullptr;
		Vector trans=Vector(0,0,0);
		bool staticroot=false;
		std::vector<CircularArray<Atom*>> rings;
	public:
		/*** Iterator ***/
		Atom* begin() {return &root;}
		Atom* end()
		{
			Atom* ret=&root;
			while(ret->getBondCount()>1)
				ret=&(ret->getLastBond()->getSecond(ret));
			return ret;
		}

		//Constructors
		//Molecule() #NOT ALLOWED#
		Molecule(Atom* at) : Molecule(*at) {}
		Molecule(Atom& r) : root(r) {current=&root;root.setMolecule(*this);recurse(current);}
		Molecule(const Molecule& mol,bool warn=true) : root(mol.root)
		{
			if(warn)
				cout << "WARN: Molecule Reconstruction. Please read instructions for using Molecule class safely\n";
			current=mol.current;
			setAllRoot(this);
		}

		void setAllRoot(Molecule* m);
	private:
		void recurse(Atom* at) {setAllRoot(this);}
	public:
		//Motion
		inline void translate(Vector v) {trans=trans+v;}
		inline Vector getTranslation() {return trans;}
		void rotateAbout(geom3D::Point3D pt,const geom3D::EulerAngle& EA);
		inline void rotateAbout(Atom* at,const geom3D::EulerAngle& EA) {rotateAbout(at->getExactLocation(),EA);}
		inline void rotate(const geom3D::EulerAngle& EA) {rotateAbout(&root,EA);}
		inline void fixRoot() {staticroot=true;}
		inline void detachRoot() {staticroot=false;}
		void setCurrent(Atom* at)
		{
			if(at->getMolecule()!=this)
				at->setMolecule(this);
			current=at;
		}
		void centralize();

		//Building of Molecule
		Atom* getRoot() {return &root;}
		Atom*& getCurrent() {return current;}
		void setTarget(Atom* at)
		{
			if(at->getMolecule()!=this)
			{
				cout << "WARN: Molecule::setTarget(Atom*): Forcing atom to change molecule.\n";
				at->setMolecule(this);
			}
			current=at;
		}
		void addAtom(Atom& a,int ord=1)
		{
			//Need to include bond angle and length calculations
			if(current==nullptr)
				throw NoAtomSelectedException();
			if(!Data::hasLength(a,*current,ord))
			{
				double blen=input<double>(String("Enter bond-length (")+a.toString()+String("-")+current->toString()+String("): "));
				Data::addLength(a,*current,blen,ord);
				//cout << "OK: "<<blen<<"A\n";
			}
			double bondlen=Data::getLength(a,*current); //Can be made more efficient
			Bond* b=new Bond(current,&a);
			b->setOrder(ord);
			current->addBond(b);
			a.addBond(b);
			a.setMolecule(*this);
			current->organize(b,bondlen);
		}
		inline void addAtom(Atom* ptr,int ord=1) {addAtom(*ptr,ord);}
		Molecule* addMolecule(Molecule* m,int ord=1)
		{
			Atom* t=m->current;
			recAdd(t,ord);
			if(!staticroot)
				current=&current->getLastBond()->getSecond(current);
			return this;
		}

		//Extensions (March2019+)
		int fillHydrogen();
		int partialFillHydrogen();

		//Extensions (Dec2019+)
		void defineRings(const std::vector<CircularArray<Atom*>>& rs) {rings=rs;}
		const std::vector<CircularArray<Atom*>>& getRings() const {return rings;}
	private:
		void recAdd(Atom* t,int ord=1)
		{
			static std::vector<Atom*> processed;
			Atom* at=new Atom(*t);
			addAtom(*at,ord);
			processed.push_back(t);
			Atom* tcur=current;
			current=at;
			//cout << t->toString() << " "<<t->getMaxBondCount() << "\n";
			for(int i=0;i<t->getMaxBondCount();i++)
			{
				if(t->getBond(i)==nullptr)
					continue;
				if(!contains(processed,&t->getBond(i)->getSecond(t)))
					recAdd(&t->getBond(i)->getSecond(t),t->getBond(i)->getOrder());
			}
			current=tcur;
		}
	public:
		void saveToPDB(String molName);
	};
	class MoleculeIterator : public flex::Iterable<Atom*>
	{
		Molecule* M=nullptr;
		std::vector<Atom*> visited;
		bool quiet=false;
	public:
		MoleculeIterator(Molecule& m,bool q=false) : flex::Iterable<Atom*>(m.begin(),nullptr) {M=&m;quiet=q;}

		void next(Atom*& atm) override
		{
			visited.push_back(atm);
			//if(atm->getBondCount()>1)
			//{
				for(int i=0;i<atm->getMaxBondCount();i++)
				{
					if(atm->getBond(i)==nullptr)
						continue;
					if(contains(visited,&(atm->getBond(i)->getSecond(atm))))
						continue;
					atm=&(atm->getBond(i)->getSecond(atm));
					return;
				}
				//cout << "WARN: Multiple bonds but only one/no atom(s) bonded to it.\n";
				Atom* base;
				try
				{
					int ind=find(visited,atm);
					if(ind==0)
					{
						atm=nullptr;
						if(!quiet)
							cout << "END\n";
						return;
					}
					base=visited[ind-1];
					//base=&(atm->getFirstBond()->getSecond(atm));
				}
				catch(IndexOutOfBoundsException ex)
				{
					atm=nullptr;
					if(!quiet)
						cout << "END\n";
					return;
				}
				/*if(atm==M->getRoot())
				{
					atm=nullptr;
					if(!quiet)
						cout << "END\n";
					return;
				}*/
				next(base);
				atm=base;

			//}

		}
	};
	class CommonAtoms
	{
	public:
		inline static Atom* H() {return new Atom("H");}
		inline static Atom* O() {return new Atom("O",2,6);}
		inline static Atom* C() {return new Atom("C",4);}
	} CA;
	class CommonGroups
	{
	public:
		static Molecule* Methyl()
		{
			Molecule* ret=new Molecule(*CA.C());
			for(int i=0;i<3;i++)
				ret->addAtom(*CA.H());
			return ret;
		}
		static Molecule* Alcohol()
		{
			Molecule* ret=new Molecule(*CA.O());
			ret->addAtom(*CA.H());
			return ret;
		}
		static Molecule* CH2()
		{
			Molecule* ret=new Molecule(*CA.C());
			ret->addAtom(*CA.H());
			ret->addAtom(*CA.H());
			return ret;
		}


	} CG;
	class CommonMolecules
	{
	public:
		static Molecule* Water()
		{
			Molecule* ret=new Molecule(*CA.O());
			ret->addAtom(*CA.H());
			ret->addAtom(*CA.H());
			return ret;
		}
		static Molecule* Methane()
		{
			Molecule* ret=CG.Methyl();
			ret->addAtom(*CA.H());
			return ret;
		}
	};

	//Other methods in molecules
	static bool bonded(Atom* a1,Atom* a2)
	{
		for(int i=0;i<a1->getMaxBondCount();i++)
		{
			if(a1->getBond(i)==nullptr)
				continue;
			if(&(a1->getBond(i)->getSecond(a1))==a2)
				return true;
		}
		return false;
	}
	static Bond* getBond(Atom* a1,Atom* a2)
	{
		for(int i=0;i<a1->getBondCount();i++)
		{
			if(a1->getBond(i)==nullptr)
				continue;
			if(&(a1->getBond(i)->getSecond(a1))==a2)
				return a1->getBond(i);
		}
		return nullptr;
	}
	static double dihedralAngle(Bond* bnd,Atom* a1,Atom* a2)
	{
		Atom *r1=&(bnd->getAtom1()), *r2=&(bnd->getAtom2());
		//cout << a1->toString() << "-" << r1->toString() << "-\t-" << r2->toString() << "-"<<a2->toString() <<"\n";
		if(bonded(r1,a1) && bonded(r2,a2)) {}
		else if(bonded(r2,a1) && bonded(r1,a2)) {return dihedralAngle(bnd,a2,a1);}
		else
			throw NoDihedralAngleException();
		Vector bondV=bnd->toVector(r1);
		Vector v1=r1->getLocation()-a1->getLocation();
		Vector v2=a2->getLocation()-r2->getLocation();
		Vector n1=v1.cross(bondV),n2=bondV.cross(v2);
		return -n1.angleWith(n2);
	}

	static void rotateBond(Bond* bond,const Bond* src,Atom* atom,const Atom* about,double angle,bool warn=true)
	{
		if(bond->getOrder()>1 && warn)
			cerr << "Warning: Rotating multiple bond.\n";
		//cout << "IN\n";
		Bond* bptr=nullptr;
		static Vector bdir=(src->toVector(about)); //Bond direction
		bool found=false;
		for(int i=0;i<atom->getMaxBondCount();i++)
		{
			bptr=atom->getBond(i);
			if(bptr==nullptr)
				continue;
			if(bptr==bond)
			{
				found=true;
				continue;
			}
			/*cout << "Bond: "<<i<<"\n";
			cout << &bptr->getAtom1() << "\n";
			cout << &bptr->getAtom2() << "\n";
			cout  <<"At: " << atom << "\n";*/
			//cout << (&bptr->getAtom1()==atom || &bptr->getAtom2()==atom) << "\n";
			rotateBond(bptr,src,&(bptr->getSecond(atom)),about,angle,warn);
		}
		if(!found && warn)
			cerr << "Rotated bond does not connect to the atom? Or one of the bonds passed is wrong. (chem/chem.h: 639 : in method rotateBond(Bond*, Bond*, Atom*, Atom*, double, bool))\n";
		cout << "Processed: "<<atom->toString()<<"\n";
		Vector pos=atom->getLocation()-about->getLocation(); // Relative position of atom
		pos=rotateVectorAbout(pos,bdir,angle); //Fails
		pos=pos+about->getLocation();
		atom->setPosition(pos);
	}
	static void rotateBond(Bond* bond,Atom* atom,double angle,bool warn=true) {rotateBond(bond,bond,atom,atom,angle,warn);}
	static void eclipseBond(Bond* bond,Atom* a1,Atom* a2) {rotateBond(bond,&(bond->getAtom2()),-dihedralAngle(bond,a1,a2));}
	static void setDihedral(Bond* bond,Atom* a1,Atom* a2,double dh) {eclipseBond(bond,a1,a2);rotateBond(bond,&(bond->getAtom2()),dh);}

	static bool validate(Molecule* m);

	static std::vector<Bond*> generateBonds(std::vector<Atom*> ats)
	{
		return std::vector<Bond*>();
	}
	static std::vector<Bond*> readFromPDB(File pdb)
	{
		//Molecule* ret;
		std::vector<Atom*> atoms;
		pdb.resetRead();
		String ln;
		String a;
		Vector pos;
		int bc;
		while((ln=pdb.nextLine())!=EOF)
		{
			if(!ln.startsWith("ATOM"))
				continue;
			a=ln.substring(77,79);
			Atom* atom=new Atom(a);
			bc=Data::getBondCount(*atom);
			if(bc<0)
			{
				bc=input<int>(String("Enter bond count for ")+a+String(": "));
				Data::addBondCount(*atom,bc);
			}
			atom->setMaxBondCount(bc);
			pos=Vector(std::stod(ln.substring(31,39)),std::stod(ln.substring(39,47)),std::stod(ln.substring(47,54)));
			//cout << "Set\n";
			atom->setExactLocation(pos);
			atoms.push_back(atom);
		}
		return generateBonds(atoms);
	}

}


//missing functions
molecules::Atom::Atom(String el,int bcount) : molecules::Atom::Atom(el,bcount,0)
{
	if(molecules::Data::hasElectronCount(*this))
		electrons=molecules::Data::getElectronCount(*this);
}
void molecules::Atom::organize(molecules::Bond* b,double le)
{
	if(K==1)
	{
		b->getSecond(this).setRelativePosition(Vector(le,0,0),*this);
		return;
	}
	EquationSystem ES(3);
	Vector temp;
	double tmag;
	bool trans=false;
	Matrix transM,revTransM;
	for(int i=0;i<K;i++)
	{
		if(bonds[i]==nullptr)
		{
			cout << "WARN: Null bond pointer in middle\n";
			continue;
		}
		if(bonds[i]==b)
			continue;
		temp=bonds[i]->toVector(this);
		tmag=temp.getMagnitude()*le;
		temp=temp/tmag;
		if(!trans)
		{
			transM=rotateXAxisTo(temp);
			revTransM=transM.getInverse();
			//cout << transM.determinant() << "\n";
			trans=true;
		}
		//cout <<"-----------------\n" << temp <<"-----------------\n"<<"\n\n";
		if(!Data::hasAngle(bonds[i]->getSecond(this),*this,b->getSecond(this),bonds[i]->getOrder(),b->getOrder()))
		{
			double ang=input<double>(String("Enter angle (")+bonds[i]->getSecond(this).toString()+String("-")+this->toString()+String("-")+b->getSecond(this).toString()+"): ");
			Data::addAngle(bonds[i]->getSecond(this),*this,b->getSecond(this),ang,bonds[i]->getOrder(),b->getOrder());
		}
		double angl=Data::getAngle(bonds[i]->getSecond(this),*this,b->getSecond(this),bonds[i]->getOrder(),b->getOrder());
		if(isDegrees())
			angl=toRadians(angl);
		ES.addEquation(Vector(transM*temp),cos(angl));
	}
	//Rotated all the vectors so that the first faces x-axis
	ES=ES.dropLastVariables();
	//cout << "\n\n\n" << ES.getMatrixForm() << "\n" << ES.getResults() <<"\n" << ES.getMatrixForm().getInverse() << "\n\n\n------------\n";
	Vector soln=ES.solve();
	//cout << soln << "\n";
	if(soln.getSize()>3)
	{
		cout << "WARN: Too many equations, and only 3 coordinates.\n";
		soln=soln.subVector(3);
	}
	Vector loc(3);
	for(int i=0;i<soln.getSize();i++)
		loc.at(i)=soln.at(i);
	//cout << loc << "\n";
	for(int i=0;i<3;i++)
	{
		if(loc.at(i)==0)
		{
			//If this assertion fails, please check if you use maths::setDegrees(true/false) properly to choose unit of angle.
			assert(le>=loc.getMagnitude());
			loc.at(i)=sqrt(sqr(le)-sqr(loc.getMagnitude()));
			break;
		}
	}
	b->getSecond(this).setRelativePosition(Vector(revTransM*loc),*this);
}
molecules::Molecule* molecules::Atom::generateMolecule() {return new Molecule(*this);}
static molecules::Bond** molecules::getOrganicChain(int l) //Bond* array of all bonds in main skeleten of chain. Molecule* can be extracted as Bond**[0]->getAtom1().getMolecule()
{
	if(l<=0) return new Bond*[0];
	Bond** ret=new Bond*[l-1];
	Atom** atoms=new Atom*[l+2];
	for(int i=0;i<l;i++)
	{
		ret[i]=nullptr;
		atoms[i]=nullptr;
	}
	atoms[l]=nullptr; atoms[l+1]=nullptr;
	Molecule* chain=molecules::CommonGroups::Methyl();
	atoms[0]=&(chain->getRoot()->getBond(0)->getSecond(chain->getRoot()));
	for(int i=1;i<l;i++)
	{
		atoms[i]=chain->getCurrent();
		chain->addMolecule(molecules::CommonGroups::CH2());
	}
	atoms[l]=chain->getCurrent();
	atoms[l+1]=&(atoms[l]->getBond(1)->getSecond(atoms[l]));
	Bond* bnd;
	double ANG=angleDegrees(180);
	for(int i=1;i<l;i++)
	{
		/*if(atoms[i]==nullptr || atoms[i-1]==nullptr || atoms[i+2]==nullptr || atoms[i+1]==nullptr)
			cout << "Null\n";*/
		bnd=molecules::getBond(atoms[i],atoms[i+1]);
		if(bnd==nullptr)
			cout << "Err. Null\n";
		setDihedral(bnd,atoms[i-1],atoms[i+2],ANG);
		ret[i-1]=bnd;
		cout << "set: "<<i<<"\n";
	}
	return ret;
}
inline geom3D::Point3D molecules::Atom::getExactLocation() const {/*cout << myMolecule->getTranslation()<< "\n";*/if(myMolecule==nullptr) return position; else return geom3D::Point3D(myMolecule->getTranslation()+position);}
inline void molecules::Atom::setExactLocation(geom3D::Point3D pt) {if(myMolecule==nullptr) position=pt.toVector(); else position=Vector((myMolecule->getTranslation())-pt.toVector()); position=position*-1;}
	//template function
template<size_t N> static std::array<molecules::Atom,N+1> operator|(std::array<molecules::Atom,N> in,molecules::Atom a)
{
	std::array<molecules::Atom,N> ret;
	for(int i=0;i<N;i++)
		ret[i]=in[i];
	return append(ret,a);
}

//Extra functions
inline molecules::MoleculeIterator atomsOf(molecules::Molecule& m,bool q=false) {return molecules::MoleculeIterator(m,q);}
inline molecules::MoleculeIterator atomsOf(molecules::Molecule*& m,bool q=false) {return molecules::MoleculeIterator(*m,q);}
void molecules::Molecule::centralize()
{
	geom3D::Point3D loc=root.getExactLocation();
	for(Atom* a : atomsOf(*this,false))
		a->setExactLocation(a->getExactLocation()-loc);
}
void molecules::Molecule::setAllRoot(molecules::Molecule* m)
{
	for(Atom* a : atomsOf(*this,true))
		a->setMolecule(m);
}
void molecules::Molecule::rotateAbout(geom3D::Point3D pt,const geom3D::EulerAngle& EA)
{
	//cout << pt << "\n";
	for(Atom* atm : atomsOf(*this,true))
		atm->setExactLocation(rotateVector(atm->getExactLocation().toVector(),EA,pt));
}
int molecules::Molecule::partialFillHydrogen()
{
	int hC=0;
	Atom* temp;
	//cout << "Fill\n";
	for(Atom* a : atomsOf(*this,true))
	{
		if(a->getFlag() || a->isAromatic())
			continue;
		while(a->getBondedElectronCount()<a->getMaxBondCount())
		{
			temp=new Atom("H");
			temp->setMolecule(this);
			a->bondTo(temp);
			hC++;
		}
	}
	return hC;
}
int molecules::Molecule::fillHydrogen()
{
	int hC=0;
	Atom* temp;
	//cout << "Fill\n";
	for(Atom* a : atomsOf(*this,true))
	{
		if(a->getFlag())
			continue;
		while(a->getBondedElectronCount()<a->getMaxBondCount())
		{
			temp=new Atom("H");
			temp->setMolecule(this);
			a->bondTo(temp);
			hC++;
		}
	}
	return hC;
}
static std::ostream& operator<<(std::ostream& os,molecules::Molecule* mol)
{
	cout << "Translation: "<< geom3D::Point3D(mol->getTranslation()) << "\n";;
	for(molecules::Atom* at : atomsOf(mol))
	{
		os << at->toString()<<"\t";
		os << at->getExactLocation() << "\n";
	}
	return os;
}
static void detailedPrint(molecules::Molecule* mol)
{
	cout << "Translation: "<< geom3D::Point3D(mol->getTranslation()) << "\n";;
	for(molecules::Atom* at : atomsOf(mol))
	{
		cout << at->toString();
		if(at->isAromatic())
			cout << "*";
		cout <<"\t"<< at->getExactLocation();
		cout << "\t bonds to\t";
		for(int i=0;i<at->getMaxBondCount();i++)
		{
			if(at->getBond(i)==nullptr)
				continue;
			cout << at->getBond(i)->getSecond(at).toString() << at->getBond(i)->getSecond(at).getExactLocation()<<","<<at->getBond(i)->getOrder()<< " ";
		}
		cout << "\n";
	}
}
static std::ostream& operator<<(std::ostream& os,molecules::Atom* mol)
{
	if(mol==nullptr)
	{
		cout << "NULL\n";
		return os;
	}
	geom3D::Point3D pt=mol->getExactLocation();
	cout << "Atom: "<< mol->toString() <<" at " << pt.x << "," << pt.y << "," << pt.z << "\tpointing at: "<<((long)mol);
	return os;
}
void molecules::Molecule::saveToPDB(String molName)
{
	File fl=File::openFileWrite();
	String temp;
	int K=0;
	geom3D::Point3D pt;
	for(Atom* at : atomsOf(*this,true))
	{
		pt=at->getExactLocation();
		temp="ATOM  ";
		temp=temp+format(toString(K),5,false)+format(at->toString()+toString(K),5,false)+String(" "); //Column 17 last occupied
		//temp=temp+String(" ")+format(at->toString(),4,true)+String(" "); //Column 17 last occupied
		temp=temp+format(molName,3,false); //Column 20 last occupied
		temp=temp+" A   1    "; //Column 30 last occupied
		temp=temp+format((pt.x),3,8,false);
		temp=temp+format((pt.y),3,8,false);
		temp=temp+format((pt.z),3,8,false); //Column 54 last occupied
		temp=temp+"  1.0";//Column 50 last occupied
		temp=temp+String(" ")*26; //Column 76 last occupied
		temp=temp+format(at->toString(),2,false);
		cout << temp <<"\n";
		fl.writeLine(temp);
		K++;
	}
	fl.close();
	cout << "Successfully written.\n";
}

/*static auto operator|(molecules::Atom in[2],molecules::Atom pl[2])
{
	int T=2;
	auto ret=new molecules::Atom[T+T];
	for(int i=0;i<T;i++)
		ret[i]=in[i];
	for(int i=T;i<T+T;i++)
		ret[i]=pl[T-i];
	return ret;
}*/
//namespace functions
static bool molecules::validate(Molecule* m)
{
	for(Atom* at : atomsOf(m,true))
	{
		if(at->getMolecule()!=m)
			return false;
		for(int i=0;i<at->getMaxBondCount();i++)
		{
			if(at->getBond(i)==nullptr) continue;
			if(at->getBond(i)->hasAtom(at)) continue;
			return false;
		}
	}
	return true;
}
//other functions
using molecules::Atom;
static int getBondCountData(String a)
{
	Atom* t=new Atom(a);
	if(molecules::Data::hasBondCount(*t))
		return molecules::Data::getBondCount(*t);
	else
	{
		int no;
		cout << "Enter bond count for '"<<a<<"' atom: ";
		cin >> no;
		molecules::Data::addBondCount(*t,no);
		return no;
	}
}
inline static int getBondCountData(char ch) {return getBondCountData(String(ch));}

/*static std::vector<CircularArray<Atom*>> getRings(molecules::Molecule* m)
{
	std::vector<CircularArray<molecules::Atom*>> ret;
	std::vector<molecules::Atom*> processed;
	CircularArray<Atom*> *p=nullptr;
	Atom* last=nullptr;
	molecules::Bond* bndp;
	for(molecules::Atom* a : atomsOf(m,true))
	{
		if(!(a->isAromatic()))
			continue;
		if(contains(processed,a))
			continue;
		processed.push_back(a);
		if(p==nullptr)
		{
			p=new CircularArray<Atom*>;
			p.push_back(a);
			last=a;
			continue;
		}
		if(last==nullptr) {cout << "Last atom is null\n";} //Still to check
		bndp=a->getBond(a->findBond(last));
		if(!bndp->isAromatic())
			continue;
	}
	return ret;
}*/
//OPERATOR OVERLOADED for SORT
static bool operator<(CircularArray<Atom*> a1,CircularArray<Atom*> a2) {return (a1.getSize()<a2.getSize());}
static bool operator>(CircularArray<Atom*> a1,CircularArray<Atom*> a2) {return (a1.getSize()>a2.getSize());}
static bool operator==(CircularArray<Atom*> a1,CircularArray<Atom*> a2)
{
	if(a1.getSize()!=a2.getSize())
		return false;
	for(int i=0;i<a1.getSize();i++)
	{
		if(!a2.contains(a1[i]))
			return false;
	}
	return true;
}
static bool hasBondOfOrder(Atom* at,int o)
{
	molecules::Bond* b;
	for(int i=0;i<at->getMaxBondCount();i++)
	{
		b=at->getBond(i);
		if(b==nullptr)
			continue;
		if(b->getOrder()==o)
			return true;
	}
	return false;
}
inline static int sumHybridizationElectronsPure(molecules::Atom* atm)
{
	//cout << atm << "\n";
	//cout <<"Atom: "<<atm->toString()<<" "<< atm->getBondCount() << "\t"<<atm->getElectronCount() << "\t"<<atm->getBondedElectronCount() << "\n";
	return (atm->getBondCount()+(int)ceil((atm->getElectronCount()-atm->getBondedElectronCount())/2.0));
	/*else
	 * {
	 *	if(aromAttachedCheck(atm))
	 *		cout <<atm <<"\t"<< atm->getBondCount()<<" "<<atm->getElectronCount()<<" "<<atm->getBondedElectronCount() <<"-"<<hc<<"\n";
}*/
	//cout << hc <<"\n-----===------\n";
	//cout << "\t"<< hc<<"\n";
}
static int sumHybridizationElectrons(molecules::Atom* atm)
{
	int hc=sumHybridizationElectronsPure(atm);
	if(atm->isAromatic() && !hasBondOfOrder(atm,2) && !hasBondOfOrder(atm,3))
	{
		hc--;
		if(atm->toString()=="O")
			cout << atm << "\n";
	}
	return hc;
}
static std::vector<CircularArray<Atom*>> subRingSolver(std::vector<CircularArray<Atom*>> rings)
{
	//cout << "SRS: solver\n";
	std::vector<CircularArray<Atom*>> ret;
	std::vector<Atom*> passed;
	CircularArray<Atom*> r1,r2;
	Atom *ta,*ta2,*tmp;
	molecules::Bond* tb;
	bool suc,suc2;
	int i1;
	for(int i=0;i<rings.size();i++)
	{
		r1=rings[i];
		for(int j=0;j<rings.size();j++)
		{
			if(i==j)
				continue;
			r2=rings[j];
			suc=true;
			for(int p=0;p<r2.getSize();p++)
			{
				if(!(r1.contains(r2[p])))
				{
					suc=false;
					break;
				}
			}
			if(!suc)
				continue;
			for(int p=0;p<r2.getSize();p++)
			{
				ta=r2[p];
				if(contains(passed,ta))
					continue;
				i1=-1;
				CircularArray<Atom*> temp;
				for(int b=0;b<ta->getMaxBondCount();b++)
				{
					if(ta->getBond(b)==nullptr)
						continue;
					tb=ta->getBond(b);
					ta2=&(tb->getSecond(ta));
					if(r1.contains(ta2) && !r2.contains(ta2))
					{
						i1=b;
						break;
					}
				}
				if(i1==-1)
					continue;
				bool trav=r1.getTraversal(ta,ta2);
				i1=-1;
				for(int l=0;l<r1.getSize();l++)
				{
					if(r1[l]==ta)
					{
						i1=l;
						break;
					}
				}
				assert(i1!=-1); //If this fails, it means that an atom of sub-ring is not part of super-ring
				if(trav)
					i1++;
				else
					i1--;
				temp.push_back(ta);
				passed.push_back(ta);
				for(int t=i1;;(trav)?t++:t--)
				{
					tmp=r1[t];
					temp.push_back(tmp);
					if(r2.contains(tmp))
						break;
				}
				std::vector<Atom*> sp=r2.shortestPath(ta,tmp);
				passed.push_back(tmp);
				for(int i=sp.size()-2;i>0;i--)
					temp.push_back(sp[i]);
				/*cout << "\t";
				for(int p1=0;p1<r1.getSize();p1++)
					cout << r1[p1]->toString();
				cout << "\n";
				cout << "\t";
				for(int p1=0;p1<r2.getSize();p1++)
					cout << r2[p1]->toString();
				cout << "\t\t";
				for(int p1=0;p1<temp.getSize();p1++)
					cout << temp[p1]->toString();
				cout << "\n\n";*/
				if(!contains(ret,temp) && !contains(rings,temp))
					ret.push_back(temp);
			}
		}
	}
	for(auto r : rings)
		ret.push_back(r);
	if(ret.size()!=rings.size())
		ret=subRingSolver(ret);
	return ret;
}
static bool ringIsSolved(CircularArray<Atom*> ring,int sP=0)
{
	Atom* ra;
	bool solved=true;
	for(int i=sP;i<ring.getSize()+sP;i++)
	{
		ra=ring[i];
		//cout<< ra <<","<<sumHybridizationElectronsPure(ra)<<" "<<hasBondOfOrder(ra,2)<<"\n";
		if(hasBondOfOrder(ra,2))
			continue;
		if(ra->getFreeElectronCount()>=2 && sumHybridizationElectronsPure(ra)>3)
		{
			if(hasBondOfOrder(ring[i-1],2) || hasBondOfOrder(ring[i+1],2))
				continue;
		}
		//cout << ":BREAK\n";
		solved=false;
		break;
	}
	cout << "\n";
	return solved;
}
static void solveAromaticRing(std::vector<CircularArray<Atom*>> rings,int ind,int sP=0)
{
	CircularArray<Atom*> ring=rings[ind];
	if(ring.getSize()<=3) //=3?
		return;
	if(sP>ring.getSize())
	{
		cout <<"Warning: probably unsolved aromatic ring\n";
		return;
	}
	static std::vector<molecules::Bond*> chb;
	if(ringIsSolved(ring))
	{
		for(Atom* a : ring)
			a->setAromatic();
		chb=std::vector<molecules::Bond*>();
		return;
	}
	/*cout << ring.getSize() << "\t";
	for(Atom* a : ring)
		cout << a->toString();
	cout <<":" << sP<<"\n";*/
	for(molecules::Bond* b : chb)
		b->setOrder(1);
	chb=std::vector<molecules::Bond*>();
	Atom *temp1=ring[sP],*temp2=ring[sP+1];
	for(int i=sP;i<ring.getSize()+sP+1;i++)
	{
		temp1=ring[i]; temp2=ring[i+1];
		if(!hasBondOfOrder(temp1,2) && !hasBondOfOrder(temp2,2) && temp1->hasFreeBond() && temp2->hasFreeBond())
		{
			temp1->getBondWith(temp2)->setOrder(2);
			chb.push_back(temp1->getBondWith(temp2));
			i++;
		}
	}
	solveAromaticRing(rings,ind,sP+1);
}
static molecules::Molecule* solveAromatics(molecules::Molecule* m,std::vector<CircularArray<Atom*>> rings)
{
	Atom *temp1,*temp2,*temp;
	molecules::Bond* bnd,*chk;
	std::vector<Atom*> passed;
	std::vector<int> passedR;
	rings=sort(subRingSolver(rings));
	cout << rings.size() << "\n";
	/*for(auto r : rings)
	{
		for(int i=0;i<r.getSize();i++)
			cout <<r[i]->toString();
		cout << "\n";
	}*/
	int ext=1;
	bool overr=false,st=true,ar=true;
	for(int I=0;I<rings.size();I++)
	{
		auto r=rings[I];
		ar=true;
		temp1=r[0];
		overr=false;
		for(int i=0;i<r.getSize();i++)
		{
			if(!(r[i]->isAromatic()))
			{
				ar=false;
				break;
			}
		}
		if(!ar)
			continue;
		/*cout << "TEST: ";
		for(int i=0;i<r.getSize();i++)
			cout <<r[i]->toString();
		cout << "\n";*/
		solveAromaticRing(rings,I,0);
	}
	cout << "End\n";
	return m;
}
static String preprocess(const String& str)
{
	/*std::vector<String> elements;
	elements.push_back("Cl");
	elements.push_back("Br");*/
	String res=str;
	return res;
}

//Helper functions begin
String nextAtom(const String& str,int& i,bool noF=false)
{
	if(isLower(str[i]))
		return str[i];
	if(!isAlphabet(str[i]))
	{
		if(str[i]=='/' || str[i]=='\\')
		{
			i++;
			return nextAtom(str,i,noF);
		}
		cout << "Warning: not an alphabet\n";
		return str[i];
	}
	//Definitely upper-case alphabet
	if(i+1>=str.getLength() || !isLower(str[i+1]))
		return str[i];
	if(noF) //nofill>0
		return String(str[i-1])+String(str[++i]);
	else
	{
		String ret=String(str[i])+String(str[i+1]);
		if(ret=="Cl" || ret=="Br") //May change list
		{
			i++;
			return ret;
		}
		return ret[0];
	}
}
Atom* nextNewAtom(const String& str,int& i,bool noF)
{
	String el=nextAtom(str,i,noF);
	bool ar=false;
	if(isLower(el[0]))
	{
		el[0]=toUpper(el[0]);
		ar=true;
	}
	Atom* t=new Atom(el,getBondCountData(el));
	if(noF)
		t->flag();
	if(ar)
		t->setAromatic();
	return t;
}
//Helper functions end

static HashMap<int , Atom*> numrefs;
static std::vector<CircularArray<Atom*>> rings;
struct LoadedMoleculeData
{
	molecules::Molecule* m;
	std::vector<CircularArray<Atom*>> rings;
	std::vector<int> closedR;
};
static LoadedMoleculeData internalLoadMolecule(const String& smiles,HashMap<int,CircularArray<Atom*>> trs=HashMap<int,CircularArray<Atom*>>())
{
	std::vector<int> closed;
	HashMap<int,CircularArray<Atom*>> temprings=trs;

	cout <<"Smiles: "<<smiles<<"\n";
	char ch;
	String elm;
	Atom *temp,*root;
	elm=String(smiles[0]);
	int start=0;
	bool fas=false; double bo=1;
	int nofill=0,ori=0;
	if(smiles[start]=='[')
	{
		nofill=1;
		start++;
		fas=true;
	}
	root=nextNewAtom(smiles,start,nofill>0);
	start++;
	molecules::Molecule* mol=new molecules::Molecule(root);
	for(CircularArray<Atom*>& l : temprings.getValues())
		l.push_back(mol->getRoot());
	molecules::Bond* bnd;
	for(int i=start;i<smiles.getLength();i++)
	{
		ch=smiles[i];
		/*for(auto r : temprings.getValues())
		{
			for(int a=0;a<r.getSize();a++)
				cout << r[a]->toString();
			cout << "\n";
		}
		cout << "\n";*/
		if(ch=='=')
		{
			bo=2;
			continue;
		}
		else if(ch=='#')
		{
			bo=3;
			continue;
		}
		else if(ch=='-' && nofill<=0)
			continue;
		else if(ch=='[')
		{
			nofill++;
			continue;
		}
		else if(ch==']')
		{
			nofill--;
			fas=false;
			continue;
		}
		else if(ch=='@')
		{
			mol->getCurrent()->flipChirality();
			continue;
		}
		else if(ch=='+' || ch=='-')
		{
			cout << "Hit"<<ch<<"\n";
			std::string str="";
			while(i+1<smiles.getLength() && isNumeric(smiles[++i]))
				str.append(1u,smiles[i]);
			i--;
			str.append("\0");
			if(str=="" || str=="\0")
				str="1";
			int n=std::stoi(str);
			cout << "\t\t\t"<<n<<"\n";
			if(ch=='+')
				mol->getCurrent()->addCharge(n);
			else
				mol->getCurrent()->addCharge(-n);
			continue;
		}
		else if((ch=='\\'))
		{
			/*cout << "SLS\n";
			if(ori==0)
			{
				ori=1;
				cout <<"NXT\n";
				continue;
			}
			else
			{
				if(ori==1)
					mol->getCurrent()->getLastBond()->setEZ(-1);
				else
					mol->getCurrent()->getLastBond()->setEZ(1);
				cout << "NXT\n";
				continue;
			}*/
			continue;
		}
		else if(ch=='/')
		{
			/*if(ori==0)
				ori=-1;
			else
			{
				if(ori==-1)
					mol->getCurrent()->getLastBond()->setEZ(-1);
				else
					mol->getCurrent()->getLastBond()->setEZ(1);
			}*/
			continue;
		}
		else if(ch=='(')
		{
			int eI=getMatchingBracket(smiles,'(',')',i);
			String subq=smiles.substring(i+1,eI);
			i=eI;
			if(subq.getLength()<=0)
				continue;
			molecules::Molecule* smol;
			char ch=subq[0];
			if(!isAlphabet(ch) && !isBracket(ch))
				subq=subq.substring(1);
			double bod=1;
			if(ch=='=')
				bod=2;
			else if(ch=='#')
				bod=3;
			LoadedMoleculeData tlmd=internalLoadMolecule(subq,temprings);
			smol=tlmd.m;
			for(int n : tlmd.closedR)
				temprings.remove(n);
			Atom* curr=mol->getCurrent();
			//cout << "Root start\n";
			smol->setAllRoot(mol);
			//cout << "Root end\n";
			molecules::Bond* nbon=nullptr;
			nbon=curr->bondTo(smol->getRoot(),bod);
			mol->setCurrent(curr);
			continue;
		}
		else if(ch==')')
		{
			cout << "WARN: extra closing bracket ')'\n";
			continue;
		}
		else if(isNumeric(ch) || ch=='%')
		{
			cout << "Enc\n";
			String str;
			if(ch=='%')
			{
				str=String(smiles[i+1])+String(smiles[i+2]);
				i+=2;
			}
			else
				str=ch;
			str=str+'\0';
			int nov=std::stoi(str);
			if(nofill>0)
			{
				//eg. CC[CH2-]
				Atom* atype=&(mol->getCurrent()->getLastBond()->getSecond(mol->getCurrent()));
				Atom* base=mol->getCurrent();
				for(int i=1;i<nov;i++)
					base->bondTo(new Atom(atype->toString(),atype->getMaxBondCount()));
				mol->setCurrent(base);
				continue;
			}
			if(!numrefs.contains(nov))
			{
				cout << nov << "\n";
				numrefs.append(nov,mol->getCurrent());
				CircularArray<Atom*> ar;
				ar.push_back(mol->getCurrent());
				temprings.append(nov,ar);
			}
			else
			{
				//cout << "Matched: "<<nov<<"\n";
				bnd=mol->getCurrent()->bondTo(numrefs[nov],bo);
				//flag=false;
				bo=1;
				numrefs.remove(nov);
				//cout << "\tAfterEnc\n";
				rings.push_back(temprings[nov]); //Read problem (Next loop)
				temprings.remove(nov); // Unreadable temprings?
				closed.push_back(nov);
				//cout <<"Ring "<<ch<<" completed\n";
			}
			cout << "\tOK\n";
			continue;
		}
		//cout << "Atom\t";
		temp=nextNewAtom(smiles,i,nofill>0);
		temp->setMolecule(mol);
		bnd=mol->getCurrent()->bondTo(temp,bo);
		//cout << temp << "\t"<< (i+1)<<"\n";
		if(nofill>0 && fas)
			mol->setCurrent(&(temp->getFirstBond()->getSecond(temp)));
		else
		{
			mol->setCurrent(temp);
			for(CircularArray<Atom*>& l : temprings.getValues())
				l.push_back(temp);
		}
		if(bo==1)
			ori=0;
		//bond order reset
		bo=1;
		if(nofill>0)
			fas=true;
	}
	cout << "Done. "<< rings.size() << " ring(s).\n\n";
	LoadedMoleculeData lmd;
	lmd.m=mol;
	lmd.rings=rings;
	lmd.closedR=closed;
	for(auto r : rings)
	{
		cout << "Ring size = "<<r.getSize() << "\n";
		for(int i=0;i<r.getSize();i++)
			cout << r[i]->toString();
		cout << "\n\n";
	}
	/*numrefs=HashMap<int,Atom*>();
	rings=std::vector<CircularArray<Atom*>>();*/
	return lmd;
}


namespace chemanalysis
{
	static molecules::Atom* CARBON=new molecules::Atom("C");
	static molecules::Atom* HYDROGEN=new molecules::Atom("H");
	static molecules::Atom* OXYGEN=new molecules::Atom("O");
	static molecules::Atom* NITROGEN=new molecules::Atom("N");
	static molecules::Atom* SULPHUR=new molecules::Atom("S");
	static molecules::Atom* PHOSPHORUS=new molecules::Atom("P");
	bool ElementMatch(molecules::Atom* src,molecules::Atom* comptemp) {return (src->toString()==comptemp->toString());}
	class AtomType
	{
	protected:
		molecules::Atom* comp;
		String name,identifier;
		bool (*checkf)(molecules::Atom*)=nullptr;

	public:
		//AtomType()  (Disallowed)
		AtomType(molecules::Atom* atm) {comp=atm; checkf=nullptr;}
		AtomType(molecules::Atom* atm,bool (*cf)(molecules::Atom*),String nam="Atom Type",String iden="AT") {comp=atm; checkf=cf;name=nam; identifier=iden;}

		//String toString() {return comp->toString();}
		String toString() const {if(identifier!="AT") return identifier+"\t"+name; else return name;}
		virtual bool satisfies(molecules::Atom* atm) const
		{
			if(checkf==nullptr)
				return ElementMatch(atm,comp);
			else
				return (ElementMatch(atm,comp) && checkf(atm));
		}
	};

	//Basic checking functions
	bool alwaysTrue(molecules::Atom* src) {return true;}
	bool isBondedTo(molecules::Atom* atm,molecules::Atom* atm2,double bo=1,bool arom=false) //requires arom?
	{
		molecules::Bond* bnd;
		for(int i=0;i<atm->getMaxBondCount();i++)
		{
			bnd=atm->getBond(i);
			if(bnd==nullptr)
				continue;
			if(ElementMatch(&(bnd->getSecond(atm)),atm2) && (bo==0 || bnd->getOrder()==bo || (bnd->isAromatic() && arom)))
				return true;
		}
		return false;
	}
	inline bool isBondedToAnyOrder(molecules::Atom* atm,molecules::Atom* atm2) {return isBondedTo(atm,atm2,0);}
	bool isOnlyBondedTo(molecules::Atom* atm,molecules::Atom* atm2,double bo=1)
	{
		molecules::Bond* bnd;
		bool ex=false;
		for(int i=0;i<atm->getMaxBondCount();i++)
		{
			bnd=atm->getBond(i);
			if(bnd==nullptr)
				continue;
			ex=true; //Has at-least one bond
			if(!(ElementMatch(&(bnd->getSecond(atm)),atm2) && (bo==0 || bnd->getOrder()==bo)))
				return false; // (true && false)->false
		}
		return ex; //(ex && true) -> ex
	}
	inline bool isOnlyBondedToAnyOrder(molecules::Atom* atm,molecules::Atom* atm2) {isOnlyBondedTo(atm,atm2,0);}
	bool isOnlyBondedTo(molecules::Atom* atm,molecules::Atom* atm2,bool (*targfx)(molecules::Atom*),double bo=1)
	{
		molecules::Bond* bnd;
		bool ex=false;
		for(int i=0;i<atm->getMaxBondCount();i++)
		{
			bnd=atm->getBond(i);
			if(bnd==nullptr)
				continue;
			ex=true; //Has at-least one bond
			if(!(ElementMatch(&(bnd->getSecond(atm)),atm2) && (bo==0 || bnd->getOrder()==bo)))
				return false; // (true && false)->false
			else
			{
				if(!targfx(&(bnd->getSecond(atm))))
					return false;
			}
		}
		return ex; //(ex && true) -> ex
	}
	inline bool isOnlyBondedToAnyOrder(molecules::Atom* atm,molecules::Atom* atm2,bool (*targfx)(molecules::Atom*)) {isOnlyBondedTo(atm,atm2,targfx,0);}
	bool isBondedTo(molecules::Atom* atm,molecules::Atom* atm2,bool (*targfx)(molecules::Atom*),double bo=1,bool arom=false) //requires arom?
	{
		molecules::Bond* bnd;
		for(int i=0;i<atm->getMaxBondCount();i++)
		{
			bnd=atm->getBond(i);
			if(bnd==nullptr)
				continue;
			if(ElementMatch(&(bnd->getSecond(atm)),atm2) && (bo==0 || bnd->getOrder()==bo || (bnd->isAromatic() && arom)))
			{
				if(targfx(&(bnd->getSecond(atm))))
					return true;
				else
					continue;
			}
		}
		return false;
	}
	inline bool isBondedToAnyOrder(molecules::Atom* atm,molecules::Atom* atm2,bool (*targfx)(molecules::Atom*)) {return isBondedTo(atm,atm2,targfx,0);}
	bool isAromaticAtom(molecules::Atom* atm) {return atm->isAromatic();}
	bool aromAttachedCheck(molecules::Atom* src)
	{
		molecules::Bond* temp;
		for(int i=0;i<src->getMaxBondCount();i++)
		{
			temp=src->getBond(i);
			if(temp==nullptr)
				continue;
			if(isAromaticAtom(&(temp->getSecond(src))))
				return true;
		}
		return false;
	}
	int countAtomsIn(molecules::Atom* atm,molecules::Atom* typ)
	{
		molecules::Bond* temp;
		int cnt=0;
		for(int i=0;i<atm->getMaxBondCount();i++)
		{
			temp=atm->getBond(i);
			if(temp==nullptr)
				continue;
			if(ElementMatch(&(temp->getSecond(atm)),typ))
				cnt++;
		}
		return cnt;
	}
	inline int countHydrogenAtoms(molecules::Atom* atm) {countAtomsIn(atm,HYDROGEN);}
	//Added Dec2019+
	bool conjCheck(molecules::Atom* src)
	{
		if(!hasBondOfOrder(src,2)) return false;
		Atom *sdb=nullptr,*tatm;
		molecules::Bond* bnd;
		for(int i=0;i<src->getMaxBondCount();i++)
		{
			bnd=src->getBond(i);
			if(bnd==nullptr) continue;
			if(bnd->getOrder()==2)
			{
				//if(sdb) return false;
				sdb=&(bnd->getSecond(src)); //Second atom of double bond
				continue;
			}
			tatm=&(bnd->getSecond(src));
			if(hasBondOfOrder(tatm,2) || tatm->getFreeElectronCount()>=2 || tatm->getCharge()!=0) return true;
		}
		for(int i=0;i<sdb->getMaxBondCount();i++)
		{
			bnd=sdb->getBond(i);
			if(bnd==nullptr) continue;
			if(&(bnd->getSecond(sdb))==src) continue;
			//if(bnd->getOrder()==2) return false;
			tatm=&(bnd->getSecond(sdb));
			//cout << tatm->toString() << " "<<hasBondOfOrder(tatm,2)<<" "<< tatm->getFreeElectronCount() <<" "<<tatm->getCharge()<<"\n";
			if(hasBondOfOrder(tatm,2) || tatm->getFreeElectronCount()>=2 || tatm->getCharge()!=0) return true;
		}
		return false;
	}
	bool isInRingCheck(molecules::Atom* atm,int ringsize=3) //Default size for ease of code (Cyclopropyl carbon definititon uses it)
	{
		std::vector<CircularArray<Atom*>> rs=atm->getMolecule()->getRings();
		for(CircularArray<Atom*> r : rs)
		{
			if(r.getSize()!=ringsize) continue;
			if(r.contains(atm)) return true;
		}
		return false;
	}

	/*Common Atom types*/
	bool sp3Check(molecules::Atom* src) {return sumHybridizationElectrons(src)==4;}
	bool sp2Check(molecules::Atom* src) {return sumHybridizationElectrons(src)==3;}
	bool spCheck(molecules::Atom* src) {return sumHybridizationElectrons(src)==2;}
	bool chiralCheck(molecules::Atom* src)  {return src->isChiral();}
	inline bool isCarbonBonded(molecules::Atom* src) {return isBondedTo(src,CARBON);}
	inline bool isOxygenBonded(molecules::Atom* src) {return isBondedTo(src,OXYGEN);}
	inline bool isOxygenBonded2(molecules::Atom* src) {return isBondedTo(src,OXYGEN,2);}
	inline bool isNitrogenBonded(molecules::Atom* src) {return isBondedTo(src,NITROGEN);}
	inline bool isOnlyNitrogenBonded(molecules::Atom* src) {return isOnlyBondedTo(src,NITROGEN);}
	inline bool isOnlyNitrogenBondedAnyOrder(molecules::Atom* src) {return isOnlyBondedToAnyOrder(src,NITROGEN);}

	namespace gaff
	{
		//Carbon
		bool CarbonylCheck(molecules::Atom* src) {return (isBondedTo(src,OXYGEN,2) || isBondedTo(src,SULPHUR,2));}

		class CarbonylCarbon : public AtomType
		{
		public:
			CarbonylCarbon() : AtomType(CARBON,&CarbonylCheck,"Carbonyl-group Carbon") {}
		} CARBOCARB;
		class AromaticCarbon : public AtomType
		{
		public:
			AromaticCarbon() : AtomType(CARBON,&isAromaticAtom,"Aromatic (SP2) Carbon") {}
		} AROMCARB;
		class SP2Carbon : public AtomType // All (including the above two)
		{
		public:
			SP2Carbon() : AtomType(CARBON,&sp2Check,"SP2 Carbon") {}
		} SP2CARB;
		class SPCarbon : public AtomType
		{
		public:
			SPCarbon() : AtomType(CARBON,&spCheck,"SP Carbon") {}
		} SPCARB;
		class SP3Carbon : public AtomType
		{
		public:
			SP3Carbon() : AtomType(CARBON,&sp3Check,"SP3 Carbon") {}
		} SP3CARB;
		class ChiralCarbon : public AtomType
		{
		public:
			ChiralCarbon() : AtomType(CARBON,&chiralCheck,"Chirality Specified (SP3) Carbon") {}
		} CHICARB;

		//Nitrogen
		bool amideNCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&CarbonylCheck);}
		bool nitroNCheck(molecules::Atom* src) {return (isBondedTo(src,OXYGEN,2) && isBondedTo(src,OXYGEN,1));}
		bool amineNCheck(molecules::Atom* src) {return !(amideNCheck(src) && nitroNCheck(src));}
		bool sp2NCheck(molecules::Atom* src) {return (sp2Check(src) && src->getBondedElectronCount()<=3);}
		bool sp3NCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondedElectronCount()<=3);}
		bool sp3NRCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondedElectronCount()==4);}
		bool sp2NRCheck(molecules::Atom* src) {return (sp2Check(src) && src->getBondedElectronCount()==4);}
		bool aromAttachedNCheck(molecules::Atom* src) {return aromAttachedCheck(src);}
		class AmideNitrogen : public AtomType
		{
		public:
			AmideNitrogen() : AtomType(NITROGEN,&amideNCheck,"Amide Nitrogen") {}
		} AMIDEN;
		class NitroNitrogen : public AtomType
		{
		public:
			NitroNitrogen() : AtomType(NITROGEN,&nitroNCheck,"Nitro-group Nitrogen") {}
		} NITRON;
		class SP1Nitrogen : public AtomType
		{
		public:
			SP1Nitrogen() : AtomType(NITROGEN,&spCheck,"SP Nitrogen") {}
		} SP1N;
		class SP2Nitrogen : public AtomType
		{
		public:
			SP2Nitrogen() : AtomType(NITROGEN,&sp2NCheck,"SP2 Nitrogen") {}
		} SP2N;
		class SP3Nitrogen : public AtomType
		{
		public:
			SP3Nitrogen() : AtomType(NITROGEN,&sp3NCheck,"SP3 Nitrogen") {}
		} SP3N;
		class SP3NitrogenHypervalent : public AtomType
		{
		public:
			SP3NitrogenHypervalent() : AtomType(NITROGEN,&sp3NRCheck,"SP3 (Hypervalent) Nitrogen") {}
		} SP3TETRAN;
		class SP2NitrogenHypervalent : public AtomType
		{
		public:
			SP2NitrogenHypervalent() : AtomType(NITROGEN,&sp2NRCheck,"SP2 (Hypervalent) Nitrogen") {}
		} SP2TRIN;
		class AromaticAttachedNitrogen : public AtomType
		{
		public:
			AromaticAttachedNitrogen() : AtomType(NITROGEN,&aromAttachedNCheck,"Aromatic-Attached Nitrogen") {}
		} AROMATTN;
		class AromaticNitrogen : public AtomType
		{
		public:
			AromaticNitrogen() : AtomType(NITROGEN,&isAromaticAtom,"Aromatic Nitrogen") {}
		} AROMN;

		//Oxygen
		bool hydroxyOCheck(molecules::Atom* src) {return (sp3Check(src) && isBondedTo(src,HYDROGEN));}
		class SP2Oxygen : public AtomType
		{
		public:
			SP2Oxygen() : AtomType(OXYGEN,&sp2Check,"SP2 Oxygen") {}
		} SP2OXY;
		class HydroylOxygen : public AtomType
		{
		public:
			HydroylOxygen() : AtomType(OXYGEN,&hydroxyOCheck,"Hydroxyl Oxygen") {}
		} OHOXY;
		class SP3Oxygen : public AtomType
		{
		public:
			SP3Oxygen() : AtomType(OXYGEN,&sp3Check,"SP3 Oxygen") {}
		} SP3OXY;


		//Sulphur
		bool thiolSCheck(molecules::Atom* src) {return (sp3Check(src) && isBondedTo(src,HYDROGEN));}
		bool sp3SCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondedElectronCount()==2);}
		bool threebondSCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondedElectronCount()==3);}
		bool fourbondSCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondedElectronCount()==4);}
		class SP2Sulphur : public AtomType
		{
		public:
			SP2Sulphur() : AtomType(SULPHUR,&sp2Check,"SP2 Sulphur") {}
		} SP2SUL;
		class ThiolSulphur : public AtomType
		{
		public:
			ThiolSulphur() : AtomType(SULPHUR,&thiolSCheck,"Thiol Sulphur") {}
		} SHSUL;
		class SP3Sulphur : public AtomType
		{
		public:
			SP3Sulphur() : AtomType(SULPHUR,&sp3SCheck,"SP3 Sulphur") {}
		} SP3SUL;
		class SP3SulphurHypervalent : public AtomType
		{
		public:
			SP3SulphurHypervalent() : AtomType(SULPHUR,&threebondSCheck,"Hypervalent (3 bond) Sulphur") {}
		} HYPERSUL;
		class SP3SulphurSupervalent : public AtomType
		{
		public:
			SP3SulphurSupervalent() : AtomType(SULPHUR,&fourbondSCheck,"Supervalent (4 bond) Sulphur") {}
		} SUPERSUL;

		//Phosphorus
		bool threebondPCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondCount()==3);}
		bool fourbondPCheck(molecules::Atom* src) {return (sp3Check(src) && src->getBondCount()==4);}
		class SP2Phosphorus : public AtomType
		{
		public:
			SP2Phosphorus() : AtomType(PHOSPHORUS,&sp2Check,"SP2 Phosphorus") {}
		} SP2PHOS;
		class SP3Phosphorus : public AtomType
		{
		public:
			SP3Phosphorus() : AtomType(PHOSPHORUS,&sp3Check,"SP3 Phosphorus") {}
		} SP3PHOS;
		class Hyper3Phosphorus : public AtomType
		{
		public:
			Hyper3Phosphorus() : AtomType(PHOSPHORUS,&threebondPCheck,"Hypervalent (3 bond) Phosphorus") {}
		} HYPER3PHOS;
		class Hyper4Phosphorus : public AtomType
		{
		public:
			Hyper4Phosphorus() : AtomType(PHOSPHORUS,&fourbondPCheck,"Hypervalent (4 bond) Phosphorus") {}
		} HYPER4PHOS;


		//Hydrogen
		bool normcarbonhydCheck(molecules::Atom* src) {return isBondedTo(src,CARBON);}
		bool aromcarbonhydCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&isAromaticAtom);}
		bool nithydCheck(molecules::Atom* src) {return isBondedTo(src,NITROGEN);}
		bool oxyhydCheck(molecules::Atom* src) {return isBondedTo(src,OXYGEN);}
		bool sulhydCheck(molecules::Atom* src) {return isBondedTo(src,SULPHUR);}
		bool phoshydCheck(molecules::Atom* src) {return isBondedTo(src,PHOSPHORUS);}
		class NormalCarbonHydrogen : public AtomType
		{
		public:
			NormalCarbonHydrogen() : AtomType(HYDROGEN,&normcarbonhydCheck,"Aliphatic-carbon attached Hydrogen") {}
		} NORMHYD;
		class AromaticCarbonHydrogen : public AtomType
		{
		public:
			AromaticCarbonHydrogen() : AtomType(HYDROGEN,&aromcarbonhydCheck,"Aromatic-carbon attached Hydrogen") {}
		} AROMHYD;
		class NitrogenHydrogen : public AtomType
		{
		public:
			NitrogenHydrogen() : AtomType(HYDROGEN,&nithydCheck,"Nitrogen attached Hydrogen") {}
		} NHYD;
		class OxygenHydrogen : public AtomType
		{
		public:
			OxygenHydrogen() : AtomType(HYDROGEN,&oxyhydCheck,"Oxygen attached Hydrogen") {}
		} OXYHYD;
		class SulphurHydrogen : public AtomType
		{
		public:
			SulphurHydrogen() : AtomType(HYDROGEN,&sulhydCheck,"Sulphur attached Hydrogen") {}
		} SULHYD;
		class PhosphorusHydrogen : public AtomType
		{
		public:
			PhosphorusHydrogen() : AtomType(HYDROGEN,&phoshydCheck,"Phosphorus attached Hydrogen") {}
		} PHOSHYD;


		//Halogens
		class FluorineAtom : public AtomType
		{
		public:
			FluorineAtom() : AtomType(new molecules::Atom("F"),&alwaysTrue,"Fluorine Atom") {}
		} FATM;
		class ChlorineAtom : public AtomType
		{
		public:
			ChlorineAtom() : AtomType(new molecules::Atom("Cl"),&alwaysTrue,"Chlorine Atom") {}
		} CLATM;
		class BromineAtom : public AtomType
		{
		public:
			BromineAtom() : AtomType(new molecules::Atom("Br"),&alwaysTrue,"Bromine Atom") {}
		} BRATM;
		class IodineAtom : public AtomType
		{
		public:
			IodineAtom() : AtomType(new molecules::Atom("I"),&alwaysTrue,"Iodine Atom") {}
		} IATM;
		static std::vector<AtomType> getStandardAtomTypes() //Read note inside
		{
			/*** Please note
			 * Counts 1 SP2 Oxygen and 1 SP3 Oxygen for each nitro group.
			 * May count aromatic nitrogen as Hypervalent (SP2) nitrogen is aromaticity is due to lone pair on the nitrogen (giving it + charge in resonance)
			 ***/
			std::vector<AtomType> ret;
			ret.push_back(CARBOCARB);
			ret.push_back(AROMCARB);
			ret.push_back(SP2CARB);
			ret.push_back(SPCARB);
			ret.push_back(CHICARB);
			ret.push_back(SP3CARB);

			ret.push_back(AMIDEN);
			ret.push_back(NITRON);
			ret.push_back(SP3TETRAN);
			ret.push_back(SP2TRIN);
			ret.push_back(AROMN);
			ret.push_back(AROMATTN);
			ret.push_back(SP1N);
			ret.push_back(SP2N);
			ret.push_back(SP3N);

			ret.push_back(OHOXY);
			ret.push_back(SP2OXY);
			ret.push_back(SP3OXY);

			ret.push_back(HYPERSUL);
			ret.push_back(SUPERSUL);
			ret.push_back(SHSUL);
			ret.push_back(SP2SUL);
			ret.push_back(SP3SUL);

			ret.push_back(HYPER3PHOS);
			ret.push_back(HYPER4PHOS);
			ret.push_back(SP2PHOS);
			ret.push_back(SP3PHOS);

			ret.push_back(AROMHYD);
			ret.push_back(NORMHYD);
			ret.push_back(NHYD);
			ret.push_back(OXYHYD);
			ret.push_back(SULHYD);
			ret.push_back(PHOSHYD);

			ret.push_back(FATM);
			ret.push_back(CLATM);
			ret.push_back(BRATM);
			ret.push_back(IATM);

			return ret;
		}
	}

	namespace charmm //cgenff
	{
		static molecules::Atom* FLUORINE=new molecules::Atom("F");
		static molecules::Atom* CHLORINE=new molecules::Atom("Cl");


		//Sulphur checks:
		bool elemSCheck(molecules::Atom* src) {return isOnlyBondedToAnyOrder(src,SULPHUR);}
		bool disulphSCheck(molecules::Atom* src) {return (isBondedTo(src,SULPHUR,&isCarbonBonded,1) && isBondedTo(src,CARBON,1));}
		bool sulphateSCheck(molecules::Atom* src) {return isBondedTo(src,OXYGEN,2);}

		bool CarbonylCheck(molecules::Atom* src) {return isBondedTo(src,OXYGEN,2);}
		bool quatSP3CCheck(molecules::Atom* src) {return (countHydrogenAtoms(src)==0 && sp3Check(src));}
		bool terSP3CCheck(molecules::Atom* src) {return (countHydrogenAtoms(src)==1 && sp3Check(src));}
		bool secSP3CCheck(molecules::Atom* src) {return (countHydrogenAtoms(src)==2 && sp3Check(src));}
		bool priSP3CCheck(molecules::Atom* src) {return (countHydrogenAtoms(src)>=3 && sp3Check(src));}
		bool isOnlyWithPriCarbon(molecules::Atom* src) {return (isOnlyBondedTo(src,CARBON,&priSP3CCheck));}
		bool tetmetammCCheck(molecules::Atom* src) {return (countHydrogenAtoms(src)>=3 && isBondedTo(src,NITROGEN,&isOnlyWithPriCarbon));}
		bool thiolateGroupSulphurCheck(molecules::Atom* src) //Assuming src to be Sulphur atom
		{
			if(src->getCharge()<0 || src->getBondedElectronCount()<2) return true;
			else {return false;}
		}
		bool thiolateCCheck(molecules::Atom* src) {return isBondedTo(src,SULPHUR,&thiolateGroupSulphurCheck);}
		bool CE1SP2Check(molecules::Atom* src) {return (isBondedTo(src,CARBON,2) && countHydrogenAtoms(src)<=1);}
		bool CE2SP2Check(molecules::Atom* src) {return (isBondedTo(src,CARBON,2) && countHydrogenAtoms(src)>=2);}
		bool CO2CCheck(molecules::Atom* src) {return (isOnlyBondedTo(src,OXYGEN,2));}
		bool cyanCCheck(molecules::Atom* src) {return isBondedTo(src,NITROGEN,3);}
		bool monofluoroCCheck(molecules::Atom* src) {return (countAtomsIn(src,FLUORINE)==1);} //Assuming sp3 carbon
		bool difluoroCCheck(molecules::Atom* src) {return (countAtomsIn(src,FLUORINE)==2);} //Assuming sp3 carbon
		bool trifluoroCCheck(molecules::Atom* src) {return (countAtomsIn(src,FLUORINE)==3);} //Assuming sp3 carbon
		bool cyclopropCCheck(molecules::Atom* src) {return isInRingCheck(src,3);} //Assuming sp3 carbon
		bool connectsCarboxyl(molecules::Atom* src)
		{
			molecules::Bond* bnd;
			Atom* temp;
			for(int i=0;i<src->getMaxBondCount();i++)
			{
				bnd=src->getBond(i);
				if(bnd==nullptr) continue;
				temp=&(bnd->getSecond(src));
				if(CarbonylCheck(temp)) return true;
			}
			return false;
		}

		class CarbonylCarbon : public AtomType
		{
		public:
			CarbonylCarbon() : AtomType(CARBON,&CarbonylCheck,"Carbonyl-group Carbon","C") {}
		} CARBOCARB;
		class AromaticCarbon : public AtomType
		{
		public:
			AromaticCarbon() : AtomType(CARBON,&isAromaticAtom,"Aromatic (SP2) Carbon","CA") {}
		} AROMCARB;
		class SP3QUATCarbon : public AtomType //SP3 Carbon with no hydrogens
		{
		public:
			SP3QUATCarbon() : AtomType(CARBON,&quatSP3CCheck,"Quaternary (SP3) Carbon","CT") {}
		} SP3QC;
		class SP3TERCarbon : public AtomType //SP3 Carbon with 1 hydrogens
		{
		public:
			SP3TERCarbon() : AtomType(CARBON,&terSP3CCheck,"Tertiary (SP3) Carbon","CT1") {}
		} SP3TC;
		class SP3SECCarbon : public AtomType //SP3 Carbon with 2 hydrogens
		{
		public:
			SP3SECCarbon() : AtomType(CARBON,&secSP3CCheck,"Secondary (SP3) Carbon","CT2") {}
		} SP3SC;
		class SP3PRICarbon : public AtomType //SP3 Carbon with 3+ hydrogens
		{
		public:
			SP3PRICarbon() : AtomType(CARBON,&priSP3CCheck,"Primary (SP3) Carbon","CT3") {}
		} SP3PC;
		class SP3WITHQUATCarbon : public AtomType //SP3 Carbon with 3 hydrogens attached to a Quaternary nitrogen
		{
		public:
			SP3WITHQUATCarbon() : AtomType(CARBON,&tetmetammCCheck,"Quaternary Nitrogen attached Primary (SP3) Carbon","CTL5") {}
		} SP3QNC;
		class ThiolateCarbon : public AtomType //Carbon attached to thiolate group
		{
		public:
			ThiolateCarbon() : AtomType(CARBON,&thiolateCCheck,"Thiolate Group attached Carbon","CS") {}
		} TLSC;
		class SP2CE1Carbon : public AtomType //SP2 Carbon with <=1 Hydrogen atoms
		{
		public:
			SP2CE1Carbon() : AtomType(CARBON,&CE1SP2Check,"(SP2) Carbon","CE1") {}
		} SP2RC;
		class SP2CE2Carbon : public AtomType //SP2 Carbon with 2 Hydrogen atoms
		{
		public:
			SP2CE2Carbon() : AtomType(CARBON,&CE2SP2Check,"Terminal (SP2) Carbon","CE2") {}
		} SP2HC;
		class ConjENECarbon : public AtomType //SP2 Carbon with conjugation
		{
		public:
			ConjENECarbon() : AtomType(CARBON,&conjCheck,"Conjugated Carbon","CC1A") {}
		} SP2CC;
		class DioxideCarbon : public AtomType //SP2 Carbon with conjugation
		{
		public:
			DioxideCarbon() : AtomType(CARBON,&CO2CCheck,"Carbon-Dioxide Carbon","CST") {}
		} CO2C;
		class CyanideCarbon : public AtomType //Cyanide Carbon
		{
		public:
			CyanideCarbon() : AtomType(CARBON,&cyanCCheck,"Cyanide Carbon","CN") {}
		} CNC;
		class MonoFluoroCarbon : public AtomType //Monofluoro Carbon
		{
		public:
			MonoFluoroCarbon() : AtomType(CARBON,&monofluoroCCheck,"Monofluoro Carbon","CF1") {}
		} CF1C;
		class DiFluoroCarbon : public AtomType //Difluoro Carbon
		{
		public:
			DiFluoroCarbon() : AtomType(CARBON,&difluoroCCheck,"Difluoro Carbon","CF2") {}
		} CF2C;
		class TriFluoroCarbon : public AtomType //Trifluoro Carbon
		{
		public:
			TriFluoroCarbon() : AtomType(CARBON,&trifluoroCCheck,"Trifluoro Carbon","CF3") {}
		} CF3C;
		class CyclopropaneCarbon : public AtomType //Cyclopropane Carbon
		{
		public:
			CyclopropaneCarbon() : AtomType(CARBON,&cyclopropCCheck,"Cyclopropane Carbon","C3") {}
		} CYPC;

		//Nitrogen
		bool amideNCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&CarbonylCheck);}
		bool amineNCheck(molecules::Atom* src) {return !(amideNCheck(src));}
		bool cyanideNCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,3);}
		bool pepNCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,connectsCarboxyl);}
		bool guanNCheck(molecules::Atom* src) {return isBondedToAnyOrder(src,CARBON,&isOnlyNitrogenBondedAnyOrder);}

		class AmideNitrogen : public AtomType //Amide Nitrogen
		{
		public:
			AmideNitrogen() : AtomType(NITROGEN,&amideNCheck,"Amide Nitrogen","NH2") {}
		} AMDN;
		class AmineNitrogen : public AtomType //Amine Nitrogen
		{
		public:
			AmineNitrogen() : AtomType(NITROGEN,&amineNCheck,"Amine Nitrogen","NH3") {}
		} ANN;
		class CyanideNitrogen : public AtomType //Cyanide Nitrogen
		{
		public:
			CyanideNitrogen() : AtomType(NITROGEN,&cyanideNCheck,"Cyanide Nitrogen","NC") {}
		} CNN;
		class PeptideNitrogen : public AtomType //Peptide Nitrogen
		{
		public:
			PeptideNitrogen() : AtomType(NITROGEN,&pepNCheck,"Peptide Nitrogen","NH1") {}
		} PEPN;
		class GuanidineNitrogen : public AtomType //Guanidine Nitrogen
		{
		public:
			GuanidineNitrogen() : AtomType(NITROGEN,&guanNCheck,"Guanidine Nitrogen","NH3L") {}
		} GNDN;

		//Oxygen
		bool hydroxyOCheck(molecules::Atom* src) {return (sp3Check(src) && isBondedTo(src,HYDROGEN));}
		bool carboxOCheck(molecules::Atom* src) {return (isBondedTo(src,CARBON,&isOxygenBonded,2) || (src->getCharge()<0 && isBondedTo(src,CARBON,&isOxygenBonded2,1)));} //Includes =O in carboxylic acid
		bool esterOCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&CarbonylCheck,1);} //The single bonded one
		bool phosEsterOCheck(molecules::Atom* src) {return isBondedTo(src,PHOSPHORUS,&CarbonylCheck,1);} //The single bonded one
		bool carbOCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,2);} //Includes =O in the ester
		bool dioxideOCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&CO2CCheck,2);}
		bool phossulOCheck(molecules::Atom* src) {return (isBondedTo(src,SULPHUR,2) || isBondedTo(src,PHOSPHORUS,2) || (src->getCharge()<0 && isBondedTo(src,SULPHUR,&sulphateSCheck,1)));}
		bool phosohOCheck(molecules::Atom* src) {return (isBondedTo(src,HYDROGEN) && isBondedTo(src,PHOSPHORUS));}
		bool acidOCheck(molecules::Atom* src) {return (isBondedTo(src,CARBON,&isOxygenBonded2,1) && isBondedTo(src,HYDROGEN));}

		class CarbonylOxygen : public AtomType //Carbonyl Oxygen
		{
		public:
			CarbonylOxygen() : AtomType(OXYGEN,&carbOCheck,"Carbonyl Oxygen","O") {}
		} CARBOO;
		class CarboxylateOxygen : public AtomType //Carboxylate Oxygen
		{
		public:
			CarboxylateOxygen() : AtomType(OXYGEN,&carboxOCheck,"Carboxylate Oxygen","OC") {}
		} CARBOXO;
		class AcidOxygen : public AtomType //Carboxylic acid Oxygen
		{
		public:
			AcidOxygen() : AtomType(OXYGEN,&acidOCheck,"Carboxylic Acid Oxygen","OB") {}
		} ACO;
		class HydroxylOxygen : public AtomType //Hydroyl Oxygen
		{
		public:
			HydroxylOxygen() : AtomType(OXYGEN,&hydroxyOCheck,"Hydroxyl Oxygen","OH1") {}
		} HYDOXYO;
		class EsterOxygen : public AtomType //Ester Oxygen
		{
		public:
			EsterOxygen() : AtomType(OXYGEN,&esterOCheck,"Ester Oxygen","OS") {}
		} ESTO;
		class DioxideOxygen : public AtomType //Carbon Dioxide's Oxygens
		{
		public:
			DioxideOxygen() : AtomType(OXYGEN,&dioxideOCheck,"Dioxide Oxygen","OST") {}
		} DIOXO;
		class PSOxide : public AtomType //Phosphate or Sulphate Oxygens
		{
		public:
			PSOxide() : AtomType(OXYGEN,&phossulOCheck,"Phosphate/Sulphate Oxygen","O2L") {}
		} PSO;
		class PHydroOxygen : public AtomType //Phosphorus hydroxyl Oxygens
		{
		public:
			PHydroOxygen() : AtomType(OXYGEN,&phosohOCheck,"Phosphorus Hydroxyl Oxygen","OHL") {}
		} POHO;
		class PEsterOxygen : public AtomType //Phosphorus hydroxyl Oxygens
		{
		public:
			PEsterOxygen() : AtomType(OXYGEN,&phosEsterOCheck,"Phosphorus Ester Oxygen","OSL") {}
		} PEL;

		//Sulphur check functions at the beginning
		//Sulphur
		class ElementalSulphur : public AtomType //Elemental Sulphur
		{
		public:
			ElementalSulphur() : AtomType(SULPHUR,&elemSCheck,"Elemental Sulphur","S") {}
		} ELS;
		class DisulphideSulphur : public AtomType //C-S-S-C type Sulphur
		{
		public:
			DisulphideSulphur() : AtomType(SULPHUR,&disulphSCheck,"Disulphide Sulphur","SM") {}
		} DSS;
		class ThiolateSulphur : public AtomType //Sulphur attached to thiolate group
		{
		public:
			ThiolateSulphur() : AtomType(SULPHUR,&thiolateGroupSulphurCheck,"Thiolate Sulphur","SS") {}
		} TLS;
		class SulphateSulphur : public AtomType //Sulphur in Sulphate
		{
		public:
			SulphateSulphur() : AtomType(SULPHUR,&sulphateSCheck,"Sulphate Sulphur","SL") {}
		} SULPHS;

		class PhosphorusAtom : public AtomType //Sulphur in Sulphate
		{
		public:
			PhosphorusAtom() : AtomType(PHOSPHORUS,&alwaysTrue,"Phosphorus","PL") {}
		} PATM;

		//Halogens
		bool monofluoroFCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&monofluoroCCheck,1);}
		bool difluoroFCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&difluoroCCheck,1);}
		bool trifluoroFCheck(molecules::Atom* src) {return isBondedTo(src,CARBON,&trifluoroCCheck,1);}

		class ChlorineAtom : public AtomType
		{
		public:
			ChlorineAtom() : AtomType(new molecules::Atom("Cl"),&alwaysTrue,"Chlorine Atom","CLAL") {}
		} CLATM;
		class HydrogenAtom : public AtomType
		{
		public:
			HydrogenAtom() : AtomType(HYDROGEN,&alwaysTrue,"Hydrogen Atom","HALL") {}
		} HALL;
			//Fluorine
		class AromFluorine : public AtomType //Fluorine on aromatic system
		{
		public:
			AromFluorine() : AtomType(FLUORINE,&aromAttachedCheck,"Aromatic-Attached Fluorine","FA") {}
		} ARF;
		class MonofluoroFluorine : public AtomType //Monofluoro Fluorine
		{
		public:
			MonofluoroFluorine() : AtomType(FLUORINE,&monofluoroFCheck,"Monofluoro Fluorine","F1") {}
		} MFF;
		class DifluoroFluorine : public AtomType //Difluoro Fluorine
		{
		public:
			DifluoroFluorine() : AtomType(FLUORINE,&difluoroFCheck,"Difluoro Fluorine","F2") {}
		} DFF;
		class TrifluoroFluorine : public AtomType //Difluoro Fluorine
		{
		public:
			TrifluoroFluorine() : AtomType(FLUORINE,&trifluoroFCheck,"Trifluoro Fluorine","F3") {}
		} TFF;



		static std::vector<AtomType> getStandardAtomTypes() //Read note inside
		{
			/*** Please note
			 * Requires S to have negative charge (or only one bond) to count as thiolate
			 * May count aromatic nitrogen as Hypervalent (SP2) nitrogen is aromaticity is due to lone pair on the nitrogen (giving it + charge in resonance)
			 ***/
			std::vector<AtomType> ret;
			ret.push_back(CO2C);
			ret.push_back(CARBOCARB);
			ret.push_back(AROMCARB);
			ret.push_back(CNC);
			ret.push_back(TLSC);
			ret.push_back(SP2CC);
			ret.push_back(SP2RC);
			ret.push_back(SP2HC);
			//add more carbon types above this line
			//Below this: SP3 confirmed
			ret.push_back(SP3QNC);
			ret.push_back(CF3C);
			ret.push_back(CF2C);
			ret.push_back(CF1C);
			ret.push_back(CYPC); //Even considers carbon atoms in heterocyclic 3 membered rings. eg. C1CO1
			ret.push_back(SP3QC);
			ret.push_back(SP3TC);
			ret.push_back(SP3SC);
			ret.push_back(SP3PC);

			//Nitrogen
			ret.push_back(PEPN);
			ret.push_back(GNDN);
			ret.push_back(AMDN);
			ret.push_back(CNN);
			ret.push_back(ANN);

			//Oxygen
			//For some specific cases: read above comments (where check functions are declared)
			ret.push_back(PSO);
			ret.push_back(DIOXO);
			ret.push_back(POHO);
			ret.push_back(ACO);
			ret.push_back(HYDOXYO);
			ret.push_back(PEL);
			ret.push_back(CARBOXO);
			ret.push_back(ESTO);
			ret.push_back(CARBOO);

			//Sulphur
			ret.push_back(TLS);
			ret.push_back(SULPHS);
			ret.push_back(ELS);
			ret.push_back(DSS);

			//Phosphorus
			ret.push_back(PATM);
			ret.push_back(CLATM);
			ret.push_back(ARF);
			ret.push_back(MFF);
			ret.push_back(DFF);
			ret.push_back(TFF);

			//All hydrogens have been grouped here
			ret.push_back(HALL);


			return ret;
		}
	}
}

static molecules::Molecule* loadMolecule(const String& smiles)
{
	cout << smiles << "\n";
	numrefs=HashMap<int , Atom*>();
	rings=std::vector<CircularArray<Atom*>>();
	String smilesP=preprocess(smiles);
	LoadedMoleculeData m = internalLoadMolecule(smilesP);
	molecules::Molecule* ret=solveAromatics(m.m,m.rings);
	ret->defineRings(rings);
	return ret;
}
static std::vector<molecules::Molecule*> loadMolecules(const String& smiles)
{
	std::vector<molecules::Molecule*> mols;
	int c=smiles.countChar('.')+1;
	for(int i=1;i<=c;i++)
		mols.push_back(loadMolecule(String::cut(smiles,'.',i)));
	return mols;
}
/*
 * else if(isNumeric(ch) || ch=='%')
		{
			String str;
			if(ch=='%')
			{
				str=String(smiles[i+1])+smiles[i+2];
				i+=2;
			}
			else
				str=ch;
			str=str+'\0';
			int nov=std::stoi(str);
			if(nofill>0)
			{
				Atom* atype=&(mol->getCurrent()->getLastBond()->getSecond(mol->getCurrent()));
				Atom* base=mol->getCurrent();
				for(int i=1;i<nov;i++)
					base->bondTo(new Atom(atype->toString(),atype->getMaxBondCount()));
				mol->setCurrent(base);
				continue;
			}
			if(!numrefs.contains(nov))
			{
				numrefs.append(nov,mol->getCurrent());
				CircularArray<Atom*> ar;
				ar.push_back(mol->getCurrent());
				temprings.append(nov,ar);
			}
			else
			{
				//cout << "Matched: "<<bo<<"\n";
				bnd=mol->getCurrent()->bondTo(numrefs[nov],bo);
				//flag=false;
				if((int)bo!=bo)
					bnd->setAromatic();
				bo=1;
				numrefs.remove(nov);
				rings.push_back(temprings[nov]);
				temprings.remove(nov);
				//cout <<"Ring "<<ch<<" completed\n";
			}
			continue;
		}
		else if(ch=='+' || ch=='-')
		{
			cout << "Hit"<<ch<<"\n";
			String str="";
			while(i+1<smiles.getLength() && isNumeric(smiles[++i]))
				str=str+smiles[i];
			i--;
			str=str+'\0';
			if(str=="" || str=='\0')
				str="1";
			int n=std::stoi(str);
			if(ch=='+')
				mol->getCurrent()->addCharge(n);
			else
				mol->getCurrent()->addCharge(-n);
			continue;
		}
		else if(ch=='(')
		{
			int eI=getMatchingBracket(smiles,'(',')',i);
			String subq=smiles.substring(i+1,eI);
			i=eI;
			if(subq.getLength()<=0)
				continue;
			molecules::Molecule* smol;
			char ch=subq[0];
			if(!isAlphabet(ch) && !isBracket(ch))
				subq=subq.substring(1);
			double bod=1;
			if(ch=='=')
				bod=2;
			else if(ch=='#')
				bod=3;
			else if(ch=='*')// To check
			{
				bod=1.5;
				cout << "aromatic side chain\n";
			}
			//cout << subq << "\n";
			smol=internalLoadMolecule(subq);
			Atom* curr=mol->getCurrent();
			//cout << "Root start\n";
			smol->setAllRoot(mol);
			//cout << "Root end\n";
			molecules::Bond* nbon=nullptr;
			nbon=curr->bondTo(smol->getRoot(),bod);
			if((int)bod!=bod)
			{
				nbon->setAromatic();
				for(auto l : temprings.getValues())
					l.push_back(numrefs[nov]); //branched ring eg. c1cc(ccc1)
			}
			mol->setCurrent(curr);
			//flag=false;
			continue;
		}
		else if(ch==')')
		{
			cout << "WARN: Extra ')' found.\n";
			continue;
		}

		elm=String(ch);
		if(i<smiles.getLength()-1 && isLower(smiles[i+1]))
			elm=elm+smiles[++i];
		//cout << elm << "\n";
		temp=new Atom(elm,getBondCountData(elm));
		if(nofill>=1)
		{
			temp->flag();
			if(!fas)
				temp->unrestrain();
		}
		temp->setMolecule(mol);
		bnd=mol->getCurrent()->bondTo(temp,bo);

		if((int)bo!=bo) // To be tested
			bnd->setAromatic();
		if(nofill>0 && fas)
			mol->setCurrent(&(temp->getFirstBond()->getSecond(temp)));
		else
			mol->setCurrent(temp);
 */
#endif
