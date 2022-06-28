#include "maths/graphs.h"

using namespace graphs;
int main()
{
  Vertex<std::string,int> v1("N1"),v2("N2"),v3("N3"),v4("N4"),v5("N5");
  v1.joinToVertex(v2);
  v2.joinToVertex(v3);
  v3.joinToVertex(v4);
  v4.joinToVertex(v1);
  //v4.joinToVertex(v2);
  //Edge<int,std::string> e(v1,v2);
  cout << v2 << "\n";
  Graph<std::string,int> g(v1);
  cout << g << "\n";
  auto path=algo::eulerPath(g,v1);
}
