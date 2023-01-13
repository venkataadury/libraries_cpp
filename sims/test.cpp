//#include "maths/maths.h"
#include <iostream>
#include <vector>
using namespace std;

int& testF(std::vector<int>& in)
{
  int& tmp=in[0];
  cout << tmp << "\n";
  in.erase(in.begin()/*+0*/);
  return tmp;
}
int main()
{
  std::vector<int> v;
  v.push_back(2); v.push_back(3); v.push_back(4);
  int tmp=testF(v);
  cout <<tmp<<"\n";
}
