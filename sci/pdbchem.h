#ifndef PDBCHEM_INCLUDED
#define PDBCHEM_INCLUDED 1
#include  "sci/chemdata.h"
#include "sci/chem.h"

using namespace molecules;
using namespace chemfiledata;
//using namespace chemdata;

std::vector<AtomType> stdtypes=chemanalysis::gaff::getStandardAtomTypes();
namespace pdbchem
{
	static std::vector<int>& addV(std::vector<int>& v1,std::vector<int>& v2)
	{
		if(v1.size()!=v2.size())
			cout << "WARNING: Size mismatch in vector addition: "<<v1.size()<<" vs "<<v2.size()<<"\n";
		if(v1.size()>=v2.size())
		{
			for(int j=0;j<v2.size();j++)
				v1[j]+=v2[j];
			return v1;
		}
		else
		{
			for(int j=0;j<v1.size();j++)
				v2[j]+=v1[j];
			return v2;
		}
	}
	static HashMap<String,std::vector<int>> STDCOUNTS;
	static std::vector<String> AMINO_ACIDS,DNA;
	std::vector<int> loadMF(MoleculeFile&,double ts=1);
	std::vector<int> processLines(std::vector<String> lines,bool hatm=false,double ts=1)
	{
		MoleculeFile mf(lines,hatm);
		return loadMF(mf,ts);
	}
	std::vector<int> process(String fn,double ts=1)
	{
		maths::setDegrees(true);
		//String fln;

		chemfiledata::MoleculeFile fl(fn);
		cout << "Done\n";
		return loadMF(fl);
	}
	std::vector<int> loadMF(MoleculeFile& fl,double ts)
	{
		//fl.printDetails();
		int ngen=0;
		std::vector<MoleculeModel> models=fl.getModels();
		molecules::Molecule** res=new molecules::Molecule*[0];
		for(MoleculeModel& mod : models)
		{
			molecules::Molecule** ext=mod.generateBonds(ts);
			res=merge(res,ext,ngen,mod.ngen);
			ngen+=mod.ngen;
			cout <<"\t"<< mod.getGeneratedBonds().size() << " bonds calculated.\n";
		}
		//fl.printDetails();

		int ha;
		for(int i=0;i<ngen;i++)
		{
			//detailedPrint(res[i]);
			ha = res[i]->fillHydrogen();
			//cout <<  ha << " hydrogen atoms added.\n";
		}
		std::vector<molecules::Molecule*> mol;
		for(int i=0;i<ngen;i++)
			mol.push_back(res[i]);
		std::vector<int> counts=countAtomTypes(mol,stdtypes);
		/*for(int i=0;i<stdtypes.size();i++)
		 *	cout << stdtypes[i].toString()<<((stdtypes[i].toString().endsWith("s"))?"es: ":"s: ")<<counts[i]<<"\n";*/
		std::vector<Atom*> delats;
		molecules::Molecule* ptr;
		/*for(int i=0;i<ngen;i++)
		{
			ptr=res[i];
			for(Atom* a : atomsOf(ptr))
				delats.push_back(a);
			delete ptr;
		}
		for(Atom* a : delats)
			delete a;*/
		return counts;
	}
	static std::vector<int> loadCountFromFile(const String& fl,String ex)
	{
		ex=ex+fl+".pdb";
		return process(ex);
	}
	static std::vector<int> getCountFromSMILES(const String& smiles)
	{
		std::vector<Molecule*> mol=loadMolecules(smiles);
		cout << "Loaded "<<mol.size()<<" molecules\n";
		int ha=0;
		for(auto m : mol)
			ha+=m->fillHydrogen();
		std::vector<int> counts=countAtomTypes(mol,stdtypes);
		return counts;
	}

	static void loadStdCounts()
	{
		loadCommonElectronCounts();
		loadCommonValences();
		loadCommonBondLengths();
		AMINO_ACIDS.push_back("GLY");
		AMINO_ACIDS.push_back("ALA");
		AMINO_ACIDS.push_back("VAL");
		AMINO_ACIDS.push_back("LEU");
		AMINO_ACIDS.push_back("ILE");
		AMINO_ACIDS.push_back("MET");
		AMINO_ACIDS.push_back("PHE");
		AMINO_ACIDS.push_back("TYR");
		AMINO_ACIDS.push_back("TRP");
		AMINO_ACIDS.push_back("SER");
		AMINO_ACIDS.push_back("THR");
		AMINO_ACIDS.push_back("ASN");
		AMINO_ACIDS.push_back("GLN");
		AMINO_ACIDS.push_back("CYS");
		AMINO_ACIDS.push_back("PRO");
		AMINO_ACIDS.push_back("ARG");
		AMINO_ACIDS.push_back("HIS");
		AMINO_ACIDS.push_back("LYS");
		AMINO_ACIDS.push_back("ASP");
		AMINO_ACIDS.push_back("GLU");
		AMINO_ACIDS.push_back("SEC");
		AMINO_ACIDS.push_back("TYS");
		for(String str : AMINO_ACIDS)
			STDCOUNTS.append(str,loadCountFromFile(str,"/home/venkata/cpp/practise/pdbs/"));
	}

	static void loadDNACounts()
	{
		File dnafile("/home/venkata/cpp/practise/pdbs/DNA.smiles",'r');
		String temp,ID;
		while((temp=dnafile.nextLine())!=EOF)
		{
			ID=String::cut(temp,',',1).trim();
			temp=String::cut(temp,',',2).trim();
			DNA.push_back(ID);
			STDCOUNTS.append(ID,getCountFromSMILES(temp));
		}
	}

	static void modify(std::vector<int>& counts,const String& v,int amt)
	{
		if(amt==0) return;
		int i;
		for(i=0;i<stdtypes.size();i++)
		{
			if(stdtypes[i].toString()==v)
				break;
		}
		if(i>=stdtypes.size())
			return;
		counts[i]+=amt;
	}
	inline static std::vector<int> getCounts(const String& ama) {try {return STDCOUNTS[ama];} catch(NoSuchKeyException ex) {cout <<"No such key! '"<<ama<<"'\n"; throw exception();}}
	static std::vector<int> loadProteinSequence(const std::vector<String>& v,bool fcorr=false)
	{
		if(v.size()==0)
			return std::vector<int>();
		std::vector<int> counts,tcount;
		int i=0,I;
		while(i<v.size() && !(STDCOUNTS.containsKey(v[i])))
			i++;
		I=i;
		if(i==v.size())
			return std::vector<int>();
		if(STDCOUNTS.containsKey(v[i]))
			counts=getCounts(v[i++]);
		else
			return std::vector<int>();
		for(;i<v.size();i++)
		{
			if(STDCOUNTS.containsKey(v[i]))
			{
				tcount=getCounts(v[i]);
				counts=addV(counts,tcount);
			}
			else
				cout <<"No such key! '"<<v[i]<<"'\n";
		}
		int CR=v.size()-I;
		modify(counts,chemanalysis::gaff::OHOXY.toString(),(fcorr)?-(CR):-(CR-1)); //Hydroxyl Oxygen
		modify(counts,chemanalysis::gaff::OXYHYD.toString(),(fcorr)?-CR:-(CR-1)); //Hydroxyl Hydrogen
		modify(counts,chemanalysis::gaff::SP3N.toString(),(fcorr)?-CR:-(CR-1)); //SP3 Nitrogen (Amines)
		modify(counts,chemanalysis::gaff::NHYD.toString(),(fcorr)?-CR:-(CR-1)); //SP3 Nitrogen attached hydrogens(Amines)
		modify(counts,chemanalysis::gaff::AMIDEN.toString(),(fcorr)?CR:(CR-1)); //Hydroxyl Oxygen
		return counts;
	}
	static std::vector<int> loadDNASequence(const std::vector<String>& v,bool fcorr=false)
	{
		if(v.size()==0)
			return std::vector<int>();
		std::vector<int> counts,tcount;
		int i=0,I;
		while(i<v.size() && !(STDCOUNTS.containsKey(v[i])))
			i++;
		if(i>=v.size())
		{
			cout << "All empty\n";
			return std::vector<int>();
		}
		I=i;
		if(STDCOUNTS.containsKey(v[i]))
			counts=getCounts(v[i++]);
		else
			return std::vector<int>();
		for(;i<v.size();i++)
		{
			if(STDCOUNTS.containsKey(v[i]))
			{
				tcount=getCounts(v[i]);
				counts=addV(counts,tcount);
			}
			else
				cout <<"No such key! (N)'"<<v[i]<<"'\n";
		}
		int CR=v.size()-I;
		modify(counts,chemanalysis::gaff::OHOXY.toString(),2*((fcorr)?-(CR):-(CR-1))); //Hydroxyl Oxygen
		modify(counts,chemanalysis::gaff::OXYHYD.toString(),2*((fcorr)?-CR:-(CR-1))); //Hydroxyl Hydrogen
		modify(counts,chemanalysis::gaff::SP3OXY.toString(),(fcorr)?CR:(CR-1)); //Hydroxyl Oxygen
		return counts;
	}

	/*static Molecule* getSideChain(String code)
	{
		if(code=="GLY" || code=="G") {return new Molecule(new Atom("H"));}
		else if(code=="ALA" || code=="A") {return new Molecule(new Atom("C"));}
		else if(code=="VAL" || code=="V")
		{
			Molecule* N=new Molecule*(new Atom("C"));
			Atom* t=new Atom("C"); t->setMolecule(N);
			N->getRoot()->bondTo(t);
			Atom* t=new Atom("C"); t->setMolecule(N);
			N->getRoot()->bondTo(t);
			return N;
		}
	}
	static Molecule* getAminoAcid(String code)
	{
		code=code.trim()
		Molecule* ret=new Molecule(new Atom("N"));
		Atom* rt=ret->getRoot();
		Atom* C1=new Atom("C"); C1->setMolecule(ret);
		ret->getRoot()->bondTo(C1);
		Atom* C2=new Atom("C"); C2->setMolecule(ret);
		C1->bondTo(C2);
		//ret->setCurrent(C1);
		Molecule* N=getSideChain(code);
		for(Atom* a : atomsOf(N))
			a->setMolecule(M);
		C1->bondTo(N->getRoot());
		M->setCurrent(rt);
		return M;
	}*/
}
#endif
