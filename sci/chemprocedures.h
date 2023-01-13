#ifndef INCLUDED_CHEMPROC
#define INCLUDED_CHEMPROC
#include "commons/commons.h"
#include "sci/chem.h"
using namespace molecules;
//using namespace gaff;
//NOTE: Built for gaff
namespace chemprocs
{
	class AtomsNotBondedException : public MoleculeExceptions {virtual const char* what() {return "Requested atom(s) is/are not bonded.";}};
	static std::vector<Molecule*> breakBond(Atom* a1,Atom* a2,Molecule* mol=nullptr)
	{
		if(mol==nullptr)
			mol=a1->getMolecule();
		if(a1->getMolecule()!=mol || a2->getMolecule()!=mol)
			throw AtomNotInMoleculeException();
		Bond* temp;
		int n1=-1,n2=-1;
		for(int i=0;i<a1->getMaxBondCount();i++)
		{
			if(&(a1->getBond(i)->getSecond(a1)) == a2)
			{
				temp=a1->getBond(i);
				n1=i;
				break;
			}
		}
		if(n1==-1)
			throw AtomsNotBondedException();
		for(int i=0;i<a2->getMaxBondCount();i++)
		{
			if(a2->getBond(i)==temp)
			{
				n2=i;
				break;
			}
		}
		if(n2==-1)
		{
			cout << "WARN: Seemingly corrupt structure: Bond connecting two atoms must be present in both the atom objects (but was in only one).\n";
			throw AtomsNotBondedException();
		}
		a1->removeBond(n1); a2->removeBond(n2);
		bool a1f=0,a2f=0;
		for(Atom* a : atomsOf(mol,true))
		{
			if(a==a1)
				a1f=1;
			if(a==a2)
				a2f=1;
		}
		std::vector<Molecule*> ret;
		ret.push_back(mol);
		if(a1f && a2f)
			return ret;
		if(a1f)
			ret.push_back(new Molecule(a2));
		if(a2f)
			ret.push_back(new Molecule(a1));
		delete temp;
		return ret;
	}

	static Atom* findAtomByType(Molecule* m,const chemanalysis::AtomType& tp)
	{
		for(Atom* a : atomsOf(m,true))
		{
			if(tp.satisfies(a))
				return a;
		}
		return nullptr;
	}

	class GenericChemicalOperation
	{
	protected:
		Atom* (*second)(Atom* f)=nullptr;
		const chemanalysis::AtomType t1;
		std::vector<Molecule*> (*operation)(Atom* at1,Atom* at2)=nullptr;

		void setOperation(std::vector<Molecule*> (*oper)(Atom* at1,Atom* at2)) {operation=oper;}
	public:
		GenericChemicalOperation(const chemanalysis::AtomType& tp,Atom* (*sec)(Atom* f)=nullptr)  : t1(tp) {second=sec;}

		std::vector<Molecule*> operate(molecules::Molecule* mol)
		{
			static std::vector<Molecule*> ret;
			Atom* a1=findAtomByType(mol,t1);
			if(a1==nullptr)
				return ret;
			Atom* a2=nullptr;
			if(second)
				a2=second(a1);
			cout << "Operate\n";
			std::vector<Molecule*> r=operation(a1,a2);
			for(Molecule* mp : r)
			{
				if(!contains(ret,mp))
					ret.push_back(mp);
			}
			for(Molecule* mp : ret)
				operate(mp);
			return ret;
		}
	};
}
using namespace chemprocs;
namespace chemanalysis::gaff
{
	//Amides
	namespace specificprocs
	{
		static Atom* amideCarbonyl(Atom* a)
		{
			Bond* bnd=nullptr;
			for(int i=0;i<a->getMaxBondCount();i++)
			{
				bnd=a->getBond(i);
				if(bnd==nullptr)
					continue;
				if(gaff::CARBOCARB.satisfies(&(bnd->getSecond(a))))
					return &(bnd->getSecond(a));
			}
			return nullptr;
		}
		static std::vector<Molecule*> amideHydrolysis(Atom* C,Atom* N)
		{
			if(chemanalysis::ElementMatch(N,CARBON)) //CARBON from GAFF
				return amideHydrolysis(N,C);
			Bond* b=nullptr;
			for(int i=0;i<C->getMaxBondCount();i++)
			{
				b=C->getBond(i);
				if(b==nullptr)
					continue;
				if(&(b->getSecond(C))==N)
					break;
			}
			assert(b!=nullptr);

			std::vector<Molecule*> ret=breakBond(C,N);
			Molecule* tm=C->getMolecule();
			Atom* oxy=new Atom("O");
			oxy->bondTo(new Atom("H"));
			C->bondTo(oxy);
			N->bondTo(new Atom("H"));
			return ret;
		}
		class AmideCleavage : public GenericChemicalOperation
		{
		public:
			AmideCleavage() : GenericChemicalOperation(AMIDEN,&amideCarbonyl) {setOperation(&amideHydrolysis);} //AMIDEN from GAFF
		} AMIDECLEAVAGE;
	}
}
#endif
