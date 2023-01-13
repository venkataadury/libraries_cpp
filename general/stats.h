#ifndef INCLUDED_GENSTATS
#define INCLUDED_GENSTATS 1
#include "general/general.h"

namespace general
{
  namespace stats
  {
    inline static double matrixAverage(const NumericMatrix<double>& m) {return m.average();}
    inline static double matrixMax(const NumericMatrix<double>& m)
    {
      if(m.getRowCount()==0 || m.getColumnCount()==0) return 0;
      double M=m[0][0];
      for(int i=0;i<m.getRowCount();i++)
      {
        for(int j=0;j<m.getColumnCount();j++) {if(m[i][j]>M) M=m[i][j];}
      }
      return M;
    }
    inline static double matrixMin(const NumericMatrix<double>& m)
    {
      if(m.getRowCount()==0 || m.getColumnCount()==0) return 0;
      double M=m[0][0];
      for(int i=0;i<m.getRowCount();i++)
      {
        for(int j=0;j<m.getColumnCount();j++) {if(m[i][j]<M) M=m[i][j];}
      }
      return M;
    }
  }
}
#endif
