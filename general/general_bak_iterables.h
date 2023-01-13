#include "maths/maths.h"
//#include "commons/flex.h"
#include <cmath>

namespace general
{
  class NonUniformMatrixException : public exception
  {
  public:
    NonUniformMatrixException()
    {
      cout <<"ERROR: General Matrix construction requested with uneven rows (all rows don't have same column count)"<<"\n";
    }
    NonUniformMatrixException(int s)
    {
      cout <<"ERROR: General Matrix construction requested with uneven rows (all rows don't have same column count: "<<s<<")"<<"\n";
    }
    NonUniformMatrixException(int i,int s)
    {
      cout <<"ERROR: General Matrix construction requested with uneven rows (all rows don't have same column count! Found "<<i<<" and "<<s<<" both as potential column counts)"<<"\n";
    }
  };
  namespace matrices
  {
    template<class T> class MatrixRow : public flex::Iterable<T>
    {
      std::vector<T> elems;
    public:
      MatrixRow() : flex::Iterable<T>(flex::IteratorPos<T>(),flex::IteratorPos<T>()) {}
      template<class U>
      MatrixRow(int s,U creator) : elems(s), flex::Iterable<T>(flex::IteratorPos<T>(elems[0],0),flex::IteratorPos<T>(elems[s],s))
      {
        //elems=std::vector<T>(s);
        for(int i=0;i<s;i++) elems[i]=T(creator);
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
      }
      MatrixRow(int s) : MatrixRow(s,T()) {}
      MatrixRow(T* const arr,int s) : elems(s), flex::Iterable<T>(flex::IteratorPos<T>(elems[0],0),flex::IteratorPos<T>(elems[s],s))
      {
        //elems=std::vector<T>(s);
        for(int i=0;i<s;i++) elems[i]=arr[i];
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
      }
      template<class U>
      MatrixRow(const std::vector<U>& arr) : elems(arr.size()), flex::Iterable<T>(flex::IteratorPos<T>(elems[0],0),flex::IteratorPos<T>(elems[arr.size()-1],arr.size()))
      {
        cout << "MatrixRow<T> constructor called\n";
        for(int i=0;i<arr.size();i++) elems[i]=T(arr[i]);
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
      }
      template<class U>
      MatrixRow(const MatrixRow<U>& r) : MatrixRow(r.elems)
      {
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
      }

      inline int size() const {return elems.size();}
      inline int getSize() const {return elems.size();}

      void next(flex::IteratorPos<T>& p) override
      {
        if(p.ind<=elems.size()) p.el=&elems[++p.ind];
      }

      void addElement(const T& t)
      {
        elems.push_back(t);
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
      }

      inline const T& operator[](int n) const
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=elems.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-elems.size(),elems.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,elems.size());
  			#endif
        if(n<elems.size()) return elems[n];
        else throw commons::IndexOutOfBoundsException(n,elems.size());
      }
      inline T& operator[](int n)
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=elems.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-elems.size(),elems.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,elems.size());
  			#endif
        if(n<elems.size()) return elems[n];
        else throw commons::IndexOutOfBoundsException(n,elems.size());
      }
      template<class U> inline bool operator==(const MatrixRow<U>& r) const
      {
        if(r.size()!=size()) return false;
        for(int i=0;i<elems.size();i++) {if(elems[i]!=r.elems[i]) return false;}
        return true;
      }
      template <class U> inline bool operator!=(const MatrixRow<U>& r) const {return !(operator==(r));}
      template<class U> MatrixRow<T>& operator=(const MatrixRow<U>& t)
      {
        elems=std::vector<T>(); for(int i=0;i<t.elems.size();i++) elems.push_back(t.elems[i]);
        flex::Iterable<T>::start=flex::IteratorPos<T>(elems[0],0);
        flex::Iterable<T>::stop=flex::IteratorPos<T>(elems[elems.size()],elems.size());
        return this;
      }
    };
    template<class T> class Matrix : public flex::Iterable<MatrixRow<T>>
    {
    protected:
      std::vector<MatrixRow<T>> rows;
      int w=0,h=0;
    public:
      Matrix() : Matrix(0,0) {cout << "Blank constructor in Matrix<T>\n";}
      Matrix(int r) : Matrix(0,r) {cout << "Blank (row) constructor in Matrix<T>\n";}
      template<class U=T>
      Matrix(int wid,int hei,const U& zero=U()) : rows(hei),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[hei],hei))
      {
        cout << "Matrix<T> row-column constructor called\n";
        w=wid; h=hei;
        for(int i=0;i<hei;i++)
        {
          std::vector<T> r;
          for(int j=0;j<wid;j++) r.push_back(T(zero));
          rows[i]=r;
        }
        if(rows.size())
        {
          flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
          flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        }
      }
      Matrix(const std::vector<T>& row) : Matrix(std::vector<std::vector<T>>(1,row)) {}
      /*Matrix(const std::vector<std::vector<T>>&  mat) : rows(mat.size()),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[mat.size()-1],mat.size()))
      {
        h=mat.size();
        for(int i=0;i<h;i++)
        {
          if(!mat[i].size()) continue;
          if(!w) w=mat[i].size();
          else {if(mat[i].size()!=w) throw NonUniformMatrixException();}
          rows[i]=mat[i];
        }
        flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
      }*/
      template<class U=T>
      Matrix(const std::vector<std::vector<T>>& mat,const U& zero=U()) : rows(mat.size()),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[mat.size()],mat.size()))
      {
        cout << "Full matrix constructor called in Matrix<T>\n";
        h=mat.size();
        for(int i=0;i<h;i++)
        {
          if(!mat[i].size()) continue;
          rows[i]=mat[i];
          if(!w) w=mat[i].size();
          else
          {
            if(w<mat[i].size()) {for(int j=0;j<i;j++) {for(int k=0;k<mat[i].size()-w;k++) rows[j].addElement(T(zero));}}
            else if(w>mat[i].size()) {for(int k=0;k<w-mat[i].size();k++) rows[i].addElement(T(zero));}
          }
        }
        if(rows.size())
        {
          flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
          flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        }
      }
      template<class U>
      Matrix(const U** matr,int r,int col) : rows(r),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[r],r))
      {
        cout << "Point constructor for matrix\n";
        h=r; w=col;
        //matr[1][0]; 2nd row 1st column
        for(int i=0;i<h;i++)
        {
          std::vector<T> row;
          for(int j=0;j<w;j++) row.push_back(T(matr[i][j]));
          rows.push_back(row);
        }
        if(rows.size())
        {
          flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
          flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        }
      }
      template<class U,class V=T> Matrix(const std::vector<MatrixRow<U>>& mrs,const V& zero=V(),bool fill=false) : rows(mrs.size()),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[mrs.size()],mrs.size()))
      {
        cout << "Rowwise constructor called in Matrix<T>\n";
        h=mrs.size();
        for(int i=0;i<h;i++)
        {
          if(!mrs[i].size()) continue;
          rows[i]=mrs[i];
          if(!w) w=mrs[i].size();
          else
          {
            if(w!=mrs[i].size())
            {
              if(fill)
              {
                if(w<mrs[i].size()) {for(int j=0;j<i;j++) {for(int k=0;k<mrs[i].size()-w;k++) rows[j].addElement(T(zero));}}
                else if(w>mrs[i].size()) {for(int k=0;k<w-mrs[i].size();k++) rows[i].addElement(T(zero));}
              }
              else throw NonUniformMatrixException();
            }
          }
        }
        if(rows.size())
        {
          flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
          flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        }
      }
      template<class U> Matrix(const Matrix<U>& m) : rows(m.rows.size()),flex::Iterable<MatrixRow<T>>(flex::IteratorPos<MatrixRow<T>>(rows[0],0),flex::IteratorPos<MatrixRow<T>>(rows[m.rows.size()],m.rows.size()))
      {
        cout << "Copy constructor for Matrix<T>\n";
        w=m.w; h=m.h;
        for(int i=0;i<m.rows.size();i++) rows[i]=m.rows[i];
        if(rows.size())
        {
          flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
          flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        }
      }

      void next(flex::IteratorPos<MatrixRow<T>>& p) override {if(p.ind<=rows.size()) p.el=&rows[++p.ind];}

      inline int getSize() const {return w*h;}
      inline int size() const {return w*h;}
      inline int getRowCount() const {return h;}
      inline int getColumnCount() const {return w;}
      inline const std::vector<MatrixRow<T>>& getRows() const {return rows;}

      MatrixRow<T>& operator[](int n)
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=rows.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-rows.size(),rows.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,rows.size());
  			#endif
        if(n<rows.size()) return rows[n];
        else throw commons::IndexOutOfBoundsException(n,rows.size());
      }
      const MatrixRow<T>& operator[](int n) const
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=rows.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-rows.size(),rows.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,rows.size());
  			#endif
        if(n<rows.size()) return rows[n];
        else throw commons::IndexOutOfBoundsException(n,rows.size());
      }
      template<class U=T> Matrix<U> transpose() const
      {
        Matrix<T> ret(h,w);
        for(int i=0;i<w;i++) for(int j=0;j<h;j++) ret[i][j]=rows[j][i];
        return ret;
      }

      Matrix<T>& operator=(const Matrix<T>& m2)
      {
        cout << "Assignment operator called in Matrix<T>\n";
        h=m2.h; w=m2.h;
        rows=std::vector<MatrixRow<T>>(m2.rows.size());
        for(int i=0;i<m2.rows.size();i++) rows[i]=m2.rows[i];
        cout << "Assigned: "<<rows.size()<< " rows\n";
        flex::Iterable<MatrixRow<T>>::start=flex::IteratorPos<MatrixRow<T>>(rows[0],0);
        flex::Iterable<MatrixRow<T>>::stop=flex::IteratorPos<MatrixRow<T>>(rows[rows.size()],rows.size());
        return *this;
      }
    };
    template<class T> class NumericMatrix : public Matrix<T>
    {
    public:
      NumericMatrix() : NumericMatrix(0,0) {}
      NumericMatrix(int r) : NumericMatrix(0,r) {}
      template<class U=T> NumericMatrix(int wid,int hei,const U& zero=U()) : Matrix<T>(wid,hei,zero) {}
      template<class U> NumericMatrix(const MatrixRow<U>& r) : Matrix<T>(r) {}
      template<class U=T> NumericMatrix(const std::vector<std::vector<T>>& mat,const U& zero=U()) : Matrix<T>(mat,zero) {}
      template<class U> NumericMatrix(const U** matr,int r,int col) : Matrix<T>(matr,r,col) {}
      template<class U,class V=T> NumericMatrix(const std::vector<MatrixRow<U>>& mrs,const V& zero=V(),bool fill=false) : Matrix<T>(mrs,zero,fill) {}
      template<class U> NumericMatrix(const Matrix<U>& m) : Matrix<T>(m) {cout << "Copy constructor (sort-of) for NumericMatrix<T>\n";}
    //private:
      template<class U> NumericMatrix(const NumericMatrix<U>& mr) :Matrix<T>(mr) {}
    public:

      //template<class U=T> inline NumericMatrix<U> transpose() const {return NumericMatrix<U>(Matrix<T>::transpose());}
      inline NumericMatrix<T>& operator=(const NumericMatrix<T>& m2) {Matrix<T>::operator=(m2); return *this;}
    };
  }
}
