#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

namespace MD
{
	class FrameHeaderUnavailableException : public std::exception {};
	namespace md_stringfx
	{
		inline std::string& string_replace(std::string& s,std::string val,std::string ns)
		{
			int l=s.find(val);
			if(l==std::string::npos) return s;
			return s.replace(l,val.size(),ns);
		}
		inline std::string string_trim(const std::string& aString) {if(aString.find_first_not_of(' ')==std::string::npos) return ""; return aString.substr(aString.find_first_not_of(' '), (aString.find_last_not_of(' ') - aString.find_first_not_of(' ')) + 1);}
	}

	class Atom
	{
		float* coords;
		std::string aname="";
	public:
		Atom(float x,float y,float z)
		{
			coords=new float[3];
			coords[0]=x;
			coords[1]=y;
			coords[2]=z;
		}
		Atom(float* r)
		{
			coords=new float[3];
			coords[0]=r[0];
			coords[1]=r[1];
			coords[2]=r[2];
		}
		~Atom() {delete[] coords;}

		inline void setName(const std::string& n) {aname=n;}
		inline const std::string& getName() const {return aname;}

		inline float getX() const {return coords[0];}
		inline float getY() const {return coords[1];}
		inline float getZ() const {return coords[2];}
	};

	class Frame
	{
		std::vector<Atom*> atoms;

		void loadFromFile(const std::string& xyzfile)
		{
			for(Atom* at : atoms) delete at;
			atoms=std::vector<Atom*>();
			std::ifstream df; df.open(xyzfile);
			loadFromStream(df);
			df.close();
		}
		void loadFromStream(std::istream& df)
		{
			int ncrd=0;
			std::string line;
			getline(df,line);
			line=md_stringfx::string_trim(line);
			while(!line.length())
			{
				if(df.eof()) throw FrameHeaderUnavailableException();
				getline(df,line);
				line=md_stringfx::string_trim(line);
				//line=line.replace(" ","").replace("\t","");
			}
			ncrd=std::stoi(line);
			std::string at;
			float x,y,z;
			while(ncrd>0)
			{
				getline(df,line);
				if(df.eof()) break;
				std::stringstream ss(line);
				ss >>at>> x >> y >> z;
				atoms.push_back(new Atom(x,y,z));
				atoms[atoms.size()-1]->setName(at);
				ncrd--;
			}
		}

	public:
		Frame(const std::string& xyzfile) {loadFromFile(xyzfile);}
		Frame(std::istream& df) {loadFromStream(df);}

		void dumpToStreamAsXYZ(std::ostream& out)
		{
			out <<" "<< atoms.size()<<"\n";
			for(int i=0;i<atoms.size();i++) out << atoms[i]->getName() <<" "<<atoms[i]->getX()<<" "<<atoms[i]->getY()<<" "<<atoms[i]->getZ()<<"\n";
		}
		void dumpToFileAsXYZ(const std::string& fname)
		{
			std::ofstream ofl; ofl.open(fname);
			dumpToStreamAsXYZ(ofl);
			ofl.close();
		}
	};

	class Trajectory
	{
		std::vector<Frame> frames;
	public:
		Trajectory(const std::string& fname)
		{
			std::ifstream infile; infile.open(fname);
			while(!infile.eof())
			{
				try{frames.push_back(Frame(infile));}
				catch(FrameHeaderUnavailableException& ex) {break;}
			}
			infile.close();
		}

		inline int getNumFrames() const {return frames.size();}
		inline Frame& getFrame(int idx) {if(idx>=0) return frames[idx]; else return frames[frames.size()+idx];}
		inline const Frame& getFrame(int idx) const {if(idx>=0) return frames[idx]; else return frames[frames.size()+idx];}
	};
}
