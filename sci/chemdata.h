#ifndef INCLUDED_CHEMDATA
#define INCLUDED_CHEMDATA 1
#include "sci/chem.h"
#include "commons/commons.h"
#include "commons/strfx.h"
#include "maths/maths.h"

using namespace commons;

class NotEnoughBondsException : public exception {virtual const char* what() {return "Not enough bonds to calculate angle.";}};

static void loadCommonValences()
{
	molecules::Data::addBondCount("H",1);
	molecules::Data::addBondCount("D",1);
	molecules::Data::addBondCount("Br",1);
	molecules::Data::addBondCount("Cl",1);
	molecules::Data::addBondCount("F",1);
	molecules::Data::addBondCount("I",1);
	molecules::Data::addBondCount("C",4);
	molecules::Data::addBondCount("O",2);
	molecules::Data::addBondCount("S",2);
	molecules::Data::addBondCount("N",3);
	molecules::Data::addBondCount("P",5);
	molecules::Data::addBondCount("As",5);
	molecules::Data::addBondCount("Se",2);
	molecules::Data::addBondCount("Te",2);
}
static void loadCommonElectronCounts() // outermost orbital only
{
	molecules::Data::addElectronCount("H",1);
	molecules::Data::addElectronCount("D",1);
	molecules::Data::addElectronCount("N",5);
	molecules::Data::addElectronCount("P",5);
	molecules::Data::addElectronCount("O",6);
	molecules::Data::addElectronCount("S",6);
	molecules::Data::addElectronCount("Se",6);
	molecules::Data::addElectronCount("Te",6);
	molecules::Data::addElectronCount("Cl",7);
	molecules::Data::addElectronCount("F",7);
	molecules::Data::addElectronCount("Br",7);
	molecules::Data::addElectronCount("I",7);
}

static void loadExtraIonicValencies()
{
	molecules::Data::addBondCount("Li",1);
	molecules::Data::addBondCount("Na",1);
	molecules::Data::addBondCount("K",1);
	molecules::Data::addBondCount("Cs",1);
	molecules::Data::addBondCount("Zn",2);
	molecules::Data::addBondCount("Cu",2);
	molecules::Data::addBondCount("Co",3);
	molecules::Data::addBondCount("Ni",2);
	molecules::Data::addBondCount("Fe",4);
	molecules::Data::addBondCount("B",3);
	molecules::Data::addBondCount("Al",3);
	molecules::Data::addBondCount("Mg",2);
	molecules::Data::addBondCount("Ca",2);
	molecules::Data::addBondCount("Sr",2);
	molecules::Data::addBondCount("Ba",2);
}
static void loadExtraIonicElectronCounts()
{
	molecules::Data::addElectronCount("As",5);
	molecules::Data::addElectronCount("Ag",9);
}
static void loadRareIonicValencies()
{
	molecules::Data::addBondCount("Sn",4);
	molecules::Data::addBondCount("Sc",3);
	molecules::Data::addBondCount("Te",2);
	molecules::Data::addBondCount("V",5);
	molecules::Data::addBondCount("Pt",4);
	molecules::Data::addBondCount("Pd",4);
	molecules::Data::addBondCount("Ga",3);
	molecules::Data::addBondCount("Mn",7);
	molecules::Data::addBondCount("Au",3);
	molecules::Data::addBondCount("Ag",2);
	molecules::Data::addBondCount("Hg",2);
	molecules::Data::addBondCount("Be",4);
	molecules::Data::addBondCount("Cd",4);
	molecules::Data::addBondCount("Cr",6);
	molecules::Data::addBondCount("Pb",2);
	molecules::Data::addBondCount("Ir",3);
	molecules::Data::addBondCount("Os",4);
	molecules::Data::addBondCount("Ru",6);
	molecules::Data::addBondCount("Mo",6);
	molecules::Data::addBondCount("Tb",6);
}

static void loadCommonBondLengths()
{
	molecules::Data::addLength("C","H",1.09);
	molecules::Data::addLength("C","N",1.48); //1.48
		molecules::Data::addTolerence("C","N",0.06);
	molecules::Data::addLength("C","N",1.28,2);
		//molecules::Data::addTolerence("C","N",0.03,2);
	molecules::Data::addLength("C","N",1.14,3);
	molecules::Data::addLength("C","C",1.52);
		molecules::Data::addTolerence("C","C",0.07);
	molecules::Data::addLength("C","C",1.34,2);
		molecules::Data::addTolerence("C","C",0.012,2);
	molecules::Data::addLength("C","C",1.2,3);
	molecules::Data::addLength("C","O",1.43); //1.42
		molecules::Data::addTolerence("C","O",0.08);
		//molecules::Data::addTolerence("C","C",0.036);
	molecules::Data::addLength("C","O",1.235,2);
		molecules::Data::addTolerence("C","O",0.029,2);
	molecules::Data::addLength("N","O",1.4); //1.4
	molecules::Data::addLength("N","O",1.21,2);
	molecules::Data::addLength("O","H",0.98);
	molecules::Data::addLength("O","O",1.48);
	molecules::Data::addLength("C","P",1.84);
		molecules::Data::addTolerence("C","P",0.05);
	molecules::Data::addLength("C","S",1.81);
		molecules::Data::addTolerence("C","S",0.05);
	molecules::Data::addLength("C","S",1.6,2);
	molecules::Data::addLength("S","O",1.62);
		molecules::Data::addTolerence("S","O",0.05);
	molecules::Data::addLength("S","O",1.45,2);
	molecules::Data::addLength("S","S",2.05);
	molecules::Data::addLength("S","H",1.34);
	molecules::Data::addLength("N","H",1.01);
	molecules::Data::addLength("C","Se",2);
	molecules::Data::addLength("C","Te",2.18);

	molecules::Data::addLength("P","O",1.60); //1.42
		molecules::Data::addTolerence("P","O",0.05);
	molecules::Data::addLength("P","O",1.45,2); //1.42
		molecules::Data::addTolerence("P","O",0.024,2);

	molecules::Data::addLength("C","F",1.35);
	molecules::Data::addLength("C","Cl",1.77);
	molecules::Data::addLength("C","Br",1.94);
	molecules::Data::addLength("C","I",2.14);

}


static int countAtoms(molecules::Molecule* mol)
{
	int c=0;
	for(Atom* a : atomsOf(mol,true))
		c++;
	return c;
}
static std::vector<int> countAtomTypes(molecules::Molecule* mol,const std::vector<chemanalysis::AtomType>& types)
{
	std::vector<int> count(types.size());
	for(int i=0;i<types.size();i++)
		count[i]=0;
	bool found=false;
	for(Atom* a : atomsOf(mol,true))
	{
		//cout << a <<"\t"<<chemanalysis::sumHybridizationElectrons(a) <<"\n";
		found=false;
		for(int i=0;i<types.size();i++)
		{
			if(types[i].satisfies(a))
			{
				/*if(types[i].toString()=="Aromatic-carbon attached Hydrogen")
					cout <<"ArC Hs: "<< a <<"\t"<< a->chain<<" "<<a->gpflag<<"\n";*/
				//cout << types[i].toString()<<"\t"<<a<<"\n";
				found=true;
				count[i]++;
				break;
			}
		}
		if(!found)
			cout << a << "\n";
	}
	if(sum(count)<countAtoms(mol))
		cout << "WARN: Some atoms have been left unaccounted. (See above this line)\n";
	return count;
}
static std::vector<int> countAtomTypes(std::vector<molecules::Molecule*> mols,const std::vector<chemanalysis::AtomType>& types)
{
	std::vector<int> count(types.size());
	for(auto m : mols)
	{
		std::vector<int> tI=countAtomTypes(m,types);
		for(int i=0;i<tI.size();i++)
			count[i]+=tI[i];
	}
	return count;
}

using molecules::Atom;
using molecules::Bond;
using maths::Vector;
using molecules::Data;
using namespace chemanalysis;

namespace chemfiledata
{
	struct AtomStatus
	{
		char chain;
		int residue=-1;
	};
	inline static bool operator==(const AtomStatus& a1,const AtomStatus& a2) {return (a1.chain==a2.chain && a1.residue==a2.residue);}
	inline static bool operator!=(const AtomStatus& a1,const AtomStatus& a2) {return !(operator==(a1,a2));}
	struct AtomDistanceMeasure
	{
		Atom* atom=nullptr;
		double dist=0;
	};
	bool satisfies(const AtomStatus& as,Atom* a)
	{
		return (a->chain==as.chain && as.residue==a->gpflag);
	}
	static double getBondCountByAngle(Atom* at)
	{
		Vector v1,v2;
		double ang=0;
		int K=0;
		for(int i=0;i<at->getMaxBondCount();i++)
		{
			if(at->getBond(i)==nullptr)
				continue;
			if(K==0)
				v1=at->getBond(i)->getSecond(at).getRelativePosition(*at);
			else
			{
				v2=at->getBond(i)->getSecond(at).getRelativePosition(*at);
				ang+=v1.angleWith(v2);
			}
			K++;
		}
		double aa;
		if(K<2)
			throw NotEnoughBondsException();
		else
			aa= ang/(K-1);
		//cout << "The Atom: "<<at << "\t";
		//cout << v1 << "\n" << v2 << "\n";
		//cout <<"Angle: "<< aa << " "<<maths::abs(aa-120.0)<<"\t";
		if(maths::abs(aa)>=116.4 && maths::abs(aa)<=122.4) // originally abs(120-aa)<=3.6 or 5
		{
			//cout << 3 << "\n";
			return 3;
		}
		else if(maths::abs(aa-Data::TETRAHEDRAL_ANGLE)<=5.5)
		{
			//cout << 4 << "\n";
			return 4;
		}
		else if(maths::abs(aa-180)<=10)
		{
			//cout << 2 << "\n";
			return 2;
		}
		else
			return 5;
	}

	class MoleculeModel
	{
		HashMap<Atom*,std::vector<AtomDistanceMeasure>> neighbours;
		String emptystring;
		double dM;
	public:
		std::vector<Atom*> atoms;
		std::vector<Bond*> bonds;
		const double DEFAULT_LENGTH=10,CUTOFF_LENGTH=3,TOL=0.026,TOLD=0.045;
		int ngen=0;

		MoleculeModel() {}

		inline std::vector<Atom*> getAtoms() {return atoms;}
		inline std::vector<Bond*> getGeneratedBonds() {return bonds;}

		MoleculeModel& operator=(const MoleculeModel& mo)
		{
			atoms=mo.atoms;
			bonds=mo.bonds;
			emptystring=mo.emptystring;
		}


		void fix()
		{
			//Filter mirrored atoms
		}
		void addAtom(Atom* at)
		{
			const double REMOVAL_DIST=1;
			bool add=true;
			AtomDistanceMeasure adm;
			std::vector<AtomDistanceMeasure> nbs;
			for(int i=0;i<atoms.size();i++)
			{
				dM=atoms[i]->getRelativePosition(*at).getMagnitude();
				if(dM<REMOVAL_DIST && (atoms[i]->toString()!=at->toString()))
				{
					add=false;
					break;
				}
				if(dM<CUTOFF_LENGTH)
				{
					adm.atom=at;
					adm.dist=dM;
					neighbours[atoms[i]].push_back(adm);
					adm.atom=atoms[i];
					nbs.push_back(adm);
				}
			}
			if(add)
			{
				atoms.push_back(at);
				neighbours.append(at,nbs);
			}
		}
		bool startLoop(Atom* a) //Atom* sr=nullptr)
		{
			static std::vector<Bond*> prcbnd;
			//cout <<"Arom: "<< a << "\n";
			Atom *sel,*chc=nullptr;
			int i=0;
			loopb: sel=nullptr;
			for(;i<a->getMaxBondCount();i++)
			{
				if(a->getBond(i)==nullptr)
					continue;
				if(hasBondOfOrder(&(a->getBond(i)->getSecond(a)),2) || !(a->getBond(i)->getSecond(a).isAromatic()))
					continue;
				if(!(a->getBond(i)->getSecond(a).hasFreeBond()))
				{
					chc=&(a->getBond(i)->getSecond(a));
					continue;
				}
				sel=&(a->getBond(i)->getSecond(a));
				i++;
				break;
				//cout << "\t" << sel << "\n";
			}
			if(sel==nullptr)
			{
				if(chc==nullptr)
					return false;
				sel=chc;
				chc=nullptr;
			}
			//cout <<"Made double bond: "<< a<< " "<<sel << "\n";
			a->getBondWith(sel)->setOrder(2);
			prcbnd.push_back(a->getBondWith(sel));
			//prcbnd[prcbnd.size()-1]->setOrder(2);
			//cout << "Success: " << hasBondOfOrder(a,2) <<"\t" << sel << "\n";
			for(int j=0;j<sel->getMaxBondCount();j++)
			{
				if(sel->getBond(j)==nullptr)
					continue;
				if(!(sel->getBond(j)->getSecond(sel).isAromatic()) || hasBondOfOrder(&(sel->getBond(j)->getSecond(sel)),2) || sel->getBond(j)->getSecond(sel).toString()!="C")
					continue;
				//cout << "Tried: " <<a<< "\t" << sel  << "\n";
				if(!startLoop(&(sel->getBond(j)->getSecond(sel))))
				{
					//cout << "Secondary: " <<a<< "\t" << &(sel->getBond(j)->getSecond(sel)) << "\n";
					int ind=find(prcbnd,a->getBondWith(sel));
					for(int I=ind;I<prcbnd.size();I++)
						prcbnd[I]->setOrder(1);
					prcbnd.erase(prcbnd.begin()+ind,prcbnd.end());
					//a->getBondWith(sel)->setOrder(1);
					goto loopb;
				}
			}
			return true;
		}
		void solveAromaticBonds()
		{
			for(Atom* a : atoms)
			{
				if(!a->isAromatic())
					continue;
				for(int i=0;i<a->getMaxBondCount();i++)
				{
					if(a->getBond(i)==nullptr)
						continue;
					if(a->getBond(i)->getSecond(a).isAromatic())
						a->getBond(i)->setOrder(1);
				}
			}
			for(int i=0;i<atoms.size();i++)
			{
				if(atoms[i]->isAromatic() && !hasBondOfOrder(atoms[i],2) && atoms[i]->toString()=="C")
					startLoop(atoms[i]);
			}
		}
		void generateHBonds()
		{
			cout <<"Started genHBonds\n";
			//bonds=std::vector<Bond*>();
			Atom* temp;
			double rd,td;
			for(int i=0;i<atoms.size();i++)
			{
				//cout <<"'"<< atoms[i]->toString() << "'\n";
				if(atoms[i]->toString()!="H")
					continue;
				temp=nullptr;
				rd=1000;
				for(int j=0;j<atoms.size();j++)
				{
					if(atoms[j]->toString()=="H")
						continue;
					td=atoms[i]->getRelativePosition(*atoms[j]).getMagnitude();
					if(td<rd)
					{
						temp=atoms[j];
						rd=td;
					}
				}
				assert(temp!=nullptr); //If this assertion fails, it means limiting value of rd=1000 was not sufficient
				bonds.push_back(atoms[i]->bondTo(temp));
			}
			cout << "Generated "<<bonds.size()<<" H bonds\n";
		}
		molecules::Molecule** generateBonds(double SCA=1)
		{
			bonds=std::vector<Bond*>();
			std::vector<Atom*> passedAts;
			//static const double TOL=0.026;
			double TOLE;
			int genB=0;
			double d,dm=1000,sdm=1000;
			Atom *cl=nullptr,*scl=nullptr;
			Bond* bp;
			int cnt=0;

			//cout << SCA << "\n";
			//cin >> d;
			for(int i=0;i<atoms.size();i++)
			{
				passedAts.push_back(atoms[i]);
				for(AtomDistanceMeasure adm : neighbours[atoms[i]])
				{
					if(!atoms[i]->hasFreeBond() || !adm.atom->hasFreeBond() || contains(passedAts,adm.atom))
						continue;
					try
					{
					d=adm.dist;
					dm=Data::getLength(*atoms[i],*(adm.atom));
					TOLE=Data::getLengthTolerence(*atoms[i],*(adm.atom));
					if(d<dm+TOLE)
					{
						bp=atoms[i]->bondTo(adm.atom);
						bonds.push_back(bp);
					}
					if(d<dm-TOLD)
					{
						dm=Data::getLength(*atoms[i],*(adm.atom),2);
						TOLE=Data::getLengthTolerence(*atoms[i],*(adm.atom),2);
						if(dm==-1)
							cout << "WARNING: Resonance where double bond data is not available\n";
						if(d<dm+TOLE)
						{
							//cout << "Forged double bond: "<<atoms[i] << " "<<(adm.atom) << "\n";
							bp->setOrder(2);
							if(d<dm-TOLD)
							{
								dm=Data::getLength(*atoms[i],*(adm.atom),3);
								if(d<dm+TOL)
									bp->setOrder(3);
							}
						}
						else
							bp->resonating=true; //bp resonates between 1 and 2
					}
					}
					catch(molecules::TooManyBondsException ex)
					{
						cout << i << "\t"<<atoms[i] << " "<<(adm.atom) << "\n";
						for(int I=0;I<atoms[i]->getMaxBondCount();I++)
						{
							if(atoms[i]->getBond(I)==nullptr)
								continue;
							cout << &(atoms[i]->getBond(I)->getSecond(atoms[i])) << "\t"<<atoms[i]->getBond(I)->getOrder()<<" "<<atoms[i]->getRelativePosition(atoms[i]->getBond(I)->getSecond(atoms[i])).getMagnitude()<<"\n";
						}
						cout << "--\n";
						for(int I=0;I<(adm.atom)->getMaxBondCount();I++)
						{
							if((adm.atom)->getBond(I)==nullptr)
								continue;
							cout << &((adm.atom)->getBond(I)->getSecond(adm.atom)) << "\t"<<(adm.atom)->getBond(I)->getOrder()<<" "<<(adm.atom)->getRelativePosition((adm.atom)->getBond(I)->getSecond(adm.atom)).getMagnitude()<<"\n";
						}
						//return nullptr;
					}
				}
			}

			//Basic bonds generated
			bool exbf=false;
			for(Atom* a : atoms)
			{
				cnt=0;
				exbf=false;
				for(int i=0;i<a->getMaxBondCount();i++)
				{
					if(a->getBond(i)==nullptr)
						continue;
					if(a->getBond(i)->resonating)
						cnt++;
					if(a->getBond(i)->getOrder()>1)
					{
						exbf=true;
						if(a->getBond(i)->getSecond(a).isAromatic())
						{
							a->setAromatic();
							break;
						}
					}
					if(cnt>=2 || (exbf && cnt>=1))
					{
						a->setAromatic();
						break;
					}
				}
			}
			bool dbf=false;
			int set=1;
			while(set>=1)
			{
				set=0;
			for(Atom* a : atoms)
			{
				cnt=0;
				dbf=false;
				for(int i=0;i<a->getMaxBondCount();i++)
				{
					if(a->getBond(i)==nullptr)
						continue;
					if(a->getBond(i)->getOrder()==2)
					{
						if(!dbf)
							dbf=true;
						else
							a->getBond(i)->setOrder(1);
					}
					if(a->getBond(i)->getSecond(a).isAromatic())
						cnt++;
				}

				if(a->isAromatic())
				{
					if(cnt<2)
					{
						set++;
						//cout << "###$ "<<a << " "<<dbf<<"\n";
						a->setAromatic(false);
						if(!dbf)
						{
							double p=0;
							Bond* sel=nullptr;
							for(int i=0;i<a->getMaxBondCount();i++)
							{
								if(a->getBond(i)==nullptr)
									continue;
								if(hasBondOfOrder(&(a->getBond(i)->getSecond(a)),2) || !a->getBond(i)->getSecond(a).hasFreeBond())
									continue;
								if(a->isAromatic() && !(a->getBond(i)->getSecond(a).isAromatic()))
									continue;
								d=a->getBond(i)->getLength()*SCA;
								dm=Data::getLength(*a,a->getBond(i)->getSecond(a),2);
								d=dm/d;
								if(d>p)
								{
									p=d;
									sel=a->getBond(i);
								}
							}
							if(sel!=nullptr)
							{
								//cout <<"Bond Angle: "<<getBondCountByAngle(a)<<"\t"<<a<<"\n";
								if(getBondCountByAngle(a)-a->getMaxBondCount()<0)
								{
									sel->setOrder(2);
									//cout <<"***$\t"<< a << "\t";
								}
							}
							else {}
						}
					}
					}

				/*else
					cout << a << "\n";*/
			}
			}

			//cin >> d;
			int nm=0;
			molecules::Molecule** mols=new molecules::Molecule*[0];
			std::vector<AtomStatus> astats;
			bool bel=false,sta=false;
			for(int i=0;i<atoms.size();i++)
			{
				bel=false;
				sta=false;
				for(AtomStatus& as : astats)
				{
					if(satisfies(as,atoms[i]))
					{
						sta=true;
						break;
					}
				}
				for(int j=0;j<nm;j++)
				{
					if(atoms[i]->getMolecule()==mols[j])
					{
						bel=true;
						break;
					}
				}
				if(!sta)
				{
					AtomStatus ns;
					ns.chain=atoms[i]->chain;
					ns.residue=atoms[i]->gpflag;
					astats.push_back(ns);
				}
				if(bel)
					continue;


				//cout << atoms[i] << "\n";
				molecules::Molecule* M = new molecules::Molecule(*atoms[i]);
				atoms[i]->recursiveSetMolecule(M);
				if(sta)
				{
					//cout << atoms[i] << "\n";
					Atom *ca1=nullptr,*ca2=nullptr;
					double tdist=CUTOFF_LENGTH;
					for(Atom* an : atomsOf(M,true))
					{
						if(!an->hasFreeBond())
							continue;
						for(AtomDistanceMeasure adm : neighbours[an])
						{
							if((an->chain!=(adm.atom)->chain) || (an->gpflag!=(adm.atom)->gpflag) || an->getMolecule()==(adm.atom)->getMolecule() || !(adm.atom)->hasFreeBond())
								continue;
							if(adm.dist>=tdist)
								continue;
							tdist=adm.dist;
							ca1=an;
							ca2=(adm.atom);
						}
					}
					if(ca1!=nullptr && ca2!=nullptr)
					{
						ca1->recursiveSetMolecule(ca2->getMolecule());
						//cout << ca1 << "\t" << ca2 << "\n";
						ca1->bondTo(ca2);
						continue;
					}
					else
						cout << "WARN: Broken residue: "<< M->getRoot()->chain<<" "<<M->getRoot()->gpflag << "\n";
				}
				mols=append<molecules::Molecule*>(mols,M,nm);
				//detailedPrint(M);
				nm++;
			}
			solveAromaticBonds();
			ngen=nm;
			Bond *obnd, *dbnd;
			for(Atom* a : atoms)
			{
				obnd=nullptr;
				dbnd=nullptr;
				if(a->isAromatic() || !hasBondOfOrder(a,2))
					continue;
				for(int i=0;i<a->getMaxBondCount();i++)
				{
					if(a->getBond(i)==nullptr)
						continue;
					if(a->getBond(i)->getOrder()>=2)
						dbnd=a->getBond(i);
					if(a->getBond(i)->getSecond(a).toString()=="O")
					{
						if(a->getBond(i)->getSecond(a).hasFreeBond())
							obnd=a->getBond(i);
					}
				}
				if(obnd!=nullptr && obnd!=dbnd)
				{
					//Tautomerize
					dbnd->setOrder(1);
					obnd->setOrder(2);
				}
			}
			return mols;
		}
		void printDetails()
		{
			for(Atom* a : atoms)
				cout << a << "\n";
		}
	};

	class MoleculeFile
	{
		std::vector<MoleculeModel> models;
		File file;
		String format;
		bool pushed=false;
		MoleculeModel tmdl;
	public:
		MoleculeFile(String fname,bool hatm=false,String form="pdb") : file(File(fname,'r'))
		{
			format=form;
			std::vector<String> FILELINES;
			String temp;
			if(format=="pdb")
			{
				while((temp=file.nextLine())!=EOF)
					FILELINES.push_back(temp);
				loadPDBFromStrings(FILELINES,hatm);
				if(!pushed)
					models.push_back(tmdl);
			}
		}
		MoleculeFile(const std::vector<String>& lines,bool hatm=false,String form="pdb")
		{
			if(form=="pdb")
			{
				loadPDBFromStrings(lines,hatm);
				if(!pushed)
					models.push_back(tmdl);
			}
		}
		~MoleculeFile()
		{
			for(MoleculeModel& mdl : models)
			{
				for(Atom* a : mdl.atoms)
					delete a;
				for(Bond* b : mdl.bonds)
					delete b;
			}
		}
		std::vector<MoleculeModel> getModels() {return models;}

		void loadPDBFromStrings(const std::vector<String>& ts,bool hatm)
		{
			String temp,chg;
			String aname;
			geom3D::Point3D loc;
			Atom* ta;
			float x,y,z;
			int charge=0;
			char nc,cc;
			for(const String& temp : ts)
			{
				if(temp.startsWith("MODEL"))
					{
						tmdl=MoleculeModel();  //ID required?
						pushed=false;
					}
					if(temp.startsWith("ENDMDL"))
					{
						models.push_back(tmdl);
						pushed=true;
					}
					if(!temp.startsWith("ATOM"))// && !temp.startsWith("HETATM"))
					{
						if(!hatm || (hatm && !temp.startsWith("HETATM")))
							continue;
					}
					if(temp.charAt(16)>='B')
						continue;
					charge=0;
					aname=temp.substring(76,78).trim();
					if(aname=="H") //Ignore Hydrogen
						continue;
					if(aname.size()>1)
						aname[1]=toLower(aname[1]);
					chg=temp.substring(78).trim();
					if(chg.size()>0)
					{
						if(chg[0]=='-' || chg[0]=='+')
						{
							nc=chg[1];
							cc=chg[0];
						}
						else
						{
							nc=chg[0];
							cc=chg[1];
						}
						std::string chgstr=""; chgstr.append(1u,nc);
						charge=std::stoi(chgstr);
						if(cc=='-')
							charge*=-1;
					}
					x=std::stof(temp.substring(31,39).toString());
					y=std::stof(temp.substring(39,47).toString());
					z=std::stof(temp.substring(47,55).toString());
					loc=geom3D::Point3D(x,y,z);
					ta=new Atom(aname,molecules::Data::getBondCount(aname));
					if(charge!=0)
						ta->addCharge(charge);
					ta->setPosition(loc);
					ta->gpflag=std::stoi(temp.substring(23,26).toString());
					//cout << temp << "\n";
					ta->chain=temp.charAt(21);
					tmdl.addAtom(ta);
			}
		}

	};
}
#endif
