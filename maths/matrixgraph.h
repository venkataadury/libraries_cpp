#ifndef INCLUDED_GRAPHMAT
#define INCLUDED_GRAPHMAT 1
#include "maths/graphs.h"
#include "general/general.h"

using graphs::Graph,graphs::Vertex,graphs::Edge;
namespace matrixgraph
{
  template<class T,class U=T> Graph<T,U> getGraphFromMatrix(const general::matrices::Matrix<T>& mymat)
  {
    int w,h=mymat.getRowCount();
    if(h==0) return Graph<T,T>();
    else w=mymat[0].getSize();
    general::matrices::Matrix<Vertex<T,U>*> newmat(w,h);
    Graph<T,U> ret;
    for(int i=0;i<h;i++)
    {
      for(int j=0;j<w;j++)
      {
        graphs::Vertex<T,U>* nver=new graphs::Vertex<T,U>(mymat[i][j]);
        newmat[i][j]=nver;
        if(i>0) nver->joinToVertex(*(newmat[i-1][j]));
        if(j>0) nver->joinToVertex(*(newmat[i][j-1]));
        ret.growGraph(*nver);
      }
    }
    return ret;
  }
}
#endif
