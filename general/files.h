#ifndef INCLUDED_FILES
#define INCLUDED_FILES 1
#include <fstream>
#include <vector>
namespace files
{
  static inline double convertString(const std::string& s,double d) {return std::stod(s);}
  static inline int convertString(const std::string& s,int i) {return std::stoi(s);}
  static inline std::string convertString(const std::string& s,std::string str) {return s;}
  //Add more to increase support

  /**@brief Load CSV files into templated vectors
    @details Loads CSV files line by line. Each line is converted into an std::vector<T> list. The type T must be specified at function call.<br/>
    A few of the starting lines (such as headers) can be skipped using skiplines option. If autoskip is chosen, lines that cannot be read (i.e. conversion to chosen data type is not possible) are skipped
  */
  template<class T> static std::vector<std::vector<T>> load_csv(const std::string& file,int skiplines=0,const std::vector<int>& cols=std::vector<int>(),bool autoskip=false)
  {
    std::ifstream loadedf; loadedf.open(file);
    std::vector<std::vector<T>> ret;
    int lno=0; std::string line;
    while(!loadedf.eof())
    {
      getline(loadedf,line);
      if(loadedf.eof()) break;
      if(lno++<skiplines) continue;
      cout << line << "\n";
      std::vector<std::string> opts=stringfx::split(line,',');
      std::vector<T> add;
      for(const std::string& n : opts) add.push_back(convertString(n,T()));
      cout <<"\n";
      ret.push_back(add);
    }
    loadedf.close();
    return ret;
  }
}
#endif
