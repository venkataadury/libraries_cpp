#ifndef INCLUDED_GENERAL
#define INCLUDED_GENERAL 1
#include "maths/maths.h"
#include "commons/flex.h"
#include "commons/commons.h"
#include <cmath>

using namespace std;
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
  class MatrixDimensionMismatchException : public exception
  {
  public:
    MatrixDimensionMismatchException() {cout << "Error matrix dimensions are mismatched for operation.\n";}
    MatrixDimensionMismatchException(std::pair<int,int> d1,std::pair<int,int> d2) {cout << "Matrix dimensions: ("<<get<0>(d1)<<","<<get<1>(d1)<<") do not match dimensions: ("<<get<0>(d2)<<","<<get<1>(d2)<<") for required operation\n";}
    MatrixDimensionMismatchException(std::pair<int,int> d1,std::pair<int,int> d2,std::string op) {cout << "Matrix dimensions: ("<<get<0>(d1)<<","<<get<1>(d1)<<") do not match dimensions: ("<<get<0>(d2)<<","<<get<1>(d2)<<") for required operation ("<<op<<")\n";}
  };
  class BadMatrixDimensionsForOperationException : public exception
  {
  public:
    BadMatrixDimensionsForOperationException() {cout << "Error matrix dimensions are invalid for requested operation.\n";}
    BadMatrixDimensionsForOperationException(std::pair<int,int> d1) {cout << "Matrix dimensions: ("<<get<0>(d1)<<","<<get<1>(d1)<<") are invalid for required operation\n";}
    BadMatrixDimensionsForOperationException(std::pair<int,int> d1,std::string op) {cout << "Matrix dimensions: ("<<get<0>(d1)<<","<<get<1>(d1)<<") are invalid for required operation ("<<op<<")\n";}
  };
  namespace matrices
  {
    template<class T> class Vector;
    template<class T> class MatrixRow
    {
      std::vector<T> elems;
    public:
      MatrixRow() {}
      template<class U> MatrixRow(int s,U creator)
      {
        elems=std::vector<T>(s);
        for(int i=0;i<s;i++) elems[i]=T(creator);
      }
      MatrixRow(int s) : MatrixRow(s,T()) {}
      MatrixRow(T* const arr,int s)
      {
        elems=std::vector<T>(s);
        for(int i=0;i<s;i++) elems[i]=arr[i];
      }
      template<class U> MatrixRow(const std::vector<U>& arr) : elems(arr.size())
      {
        for(int i=0;i<arr.size();i++) elems[i]=T(arr[i]);
      }
      template<class U> MatrixRow(const MatrixRow<U>& r) : MatrixRow(r.getElements()) {}
      template<class U> MatrixRow(const Vector<U>& v);

      inline int size() const {return elems.size();}
      inline int getSize() const {return elems.size();}
      const std::vector<T>& getElements() const {return elems;}
      std::vector<T>& getElements() {return elems;}

      inline void addElement(const T& t) {elems.push_back(t);}

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
        elems=std::vector<T>(); for(int i=0;i<t.size();i++) elems.push_back(T(t[i]));
        return *this;
      }

      operator std::vector<T>() const
      {
        std::vector<T> ret; for(int i=0;i<elems.size();i++) ret.push_back(elems[i]);
        return ret;
      }

      //Iterators
      inline flex::IndexIterator<T,MatrixRow> begin() {return flex::IndexIterator<T,MatrixRow<T>>(0,this);}
      inline flex::IndexIterator<T,MatrixRow> end() {return flex::IndexIterator<T,MatrixRow<T>>(elems.size(),this);}
    };
    template<class T> class Matrix
    {
    protected:
      std::vector<MatrixRow<T>> rows;
      int w=0,h=0;
    public:
      Matrix() : Matrix(0,0) {}
      Matrix(int r) : Matrix(0,r) {}
      template<class U=T>
      Matrix(int wid,int hei,const U& zero=U())
      {
        w=wid; h=hei;
        for(int i=0;i<hei;i++)
        {
          std::vector<T> r;
          for(int j=0;j<wid;j++) r.push_back(T(zero));
          rows.push_back(MatrixRow<T>(r));
        }
      }
      Matrix(const std::vector<T>& row) : Matrix(std::vector<std::vector<T>>(1,row)) {}
      /*Matrix(const std::vector<std::vector<T>>&  mat) : rows(mat.size())
      {
        h=mat.size();
        for(int i=0;i<h;i++)
        {
          if(!mat[i].size()) continue;
          if(!w) w=mat[i].size();
          else {if(mat[i].size()!=w) throw NonUniformMatrixException();}
          rows[i]=mat[i];
        }
      }*/
      template<class U=T>
      Matrix(const std::vector<std::vector<T>>& mat,const U& zero=U())
      {
        h=mat.size();
        for(int i=0;i<h;i++)
        {
          if(!mat[i].size()) continue;
          rows.push_back(mat[i]);
          if(!w) w=mat[i].size();
          else
          {
            if(w<mat[i].size()) {for(int j=0;j<i;j++) {for(int k=0;k<mat[i].size()-w;k++) rows[j].addElement(T(zero));}}
            else if(w>mat[i].size()) {for(int k=0;k<w-mat[i].size();k++) rows[i].addElement(T(zero));}
          }
        }
      }
      template<class U>
      Matrix(const U** matr,int r,int col)
      {
        h=r; w=col;
        //matr[1][0]; 2nd row 1st column
        for(int i=0;i<h;i++)
        {
          std::vector<T> row;
          for(int j=0;j<w;j++) row.push_back(T(matr[i][j]));
          rows.push_back(row);
        }
      }
      template<class U,class V=T> Matrix(const std::vector<MatrixRow<U>>& mrs,const V& zero=V(),bool fill=false)
      {
        h=mrs.size();
        for(int i=0;i<h;i++)
        {
          if(!mrs[i].size()) continue;
          rows.push_back(mrs[i]);
          if(!w) w=mrs[i].size();
          else
          {
            if(w!=mrs[i].size())
            {
              if(fill)
              {
                if(w<mrs[i].size()) {for(int j=0;j<i;j++) {for(int k=0;k<mrs[i].size()-w;k++) rows[j].addElement(T(zero));} w=mrs[i].size();}
                else if(w>mrs[i].size()) {for(int k=0;k<w-mrs[i].size();k++) rows[i].addElement(T(zero));}
              }
              else throw NonUniformMatrixException();
            }
          }
        }
      }
      template<class U> Matrix(const Matrix<U>& m) : rows(m.getRows().size())
      {
        this->w=m.getWidth(); this->h=m.getHeight();
        for(int i=0;i<m.getRows().size();i++) rows[i]=m.getRows()[i];
      }

      inline int getSize() const {return w*h;}
      inline int size() const {return w*h;}
      inline int getRowCount() const {return h;}
      inline int getColumnCount() const {return w;}
      inline int getWidth() const {return w;}
      inline int getHeight() const {return h;}
      inline const std::vector<MatrixRow<T>>& getRows() const {return rows;}
      inline std::vector<MatrixRow<T>> getRows(int s,int e=-1) const
      {
        std::vector<MatrixRow<T>> ret;
        if(e==-1) e=rows.size();
        for(int i=s;i<e;i++) ret.push_back(rows[i]);
        return ret;
      }

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
      template<class U> bool operator==(const Matrix<U>& m2) const
      {
        if(getRowCount()!=m2.getRowCount() || getColumnCount()!=m2.getColumnCount()) return false;
        for(int i=0;i<rows.size();i++)
        {
          for(int j=0;j<rows[i].size();j++) if(rows[i][j]!=m2[i][j]) return false;
        }
        return true;
      }

      Matrix<T>& operator=(const Matrix<T>& m2)
      {
        h=m2.h; w=m2.w;
        rows=std::vector<MatrixRow<T>>(m2.rows.size());
        for(int i=0;i<m2.rows.size();i++) rows[i]=m2.rows[i];
        return *this;
      }

      //Iterators
      inline flex::IndexIterator<MatrixRow<T>,Matrix<T>> begin() {return flex::IndexIterator<MatrixRow<T>,Matrix<T>>(0,this);}
      inline flex::IndexIterator<MatrixRow<T>,Matrix<T>> end() {return flex::IndexIterator<MatrixRow<T>,Matrix<T>>(rows.size(),this);}
    };
    template<class T=double> class NumericMatrix : public Matrix<T>
    {
    public:
      NumericMatrix() : NumericMatrix(0,0) {}
      NumericMatrix(int r) : NumericMatrix(0,r) {}
      template<class U=T> NumericMatrix(int wid,int hei,const U& zero=U()) : Matrix<T>(wid,hei,zero) {}
      template<class U> NumericMatrix(const MatrixRow<U>& r) : Matrix<T>(r) {}
      template<class U=T> NumericMatrix(const std::vector<std::vector<T>>& mat,const U& zero=U()) : Matrix<T>(mat,zero) {}
      template<class U> NumericMatrix(const U** matr,int r,int col) : Matrix<T>(matr,r,col) {}
      template<class U,class V=T> NumericMatrix(const std::vector<MatrixRow<U>>& mrs,const V& zero=V(),bool fill=false) : Matrix<T>(mrs,zero,fill) {}
      template<class U> NumericMatrix(const Matrix<U>& m) : Matrix<T>(m) {}
    public:
      template<class U=T> inline NumericMatrix<U> transpose() const {return NumericMatrix<U>(Matrix<T>::transpose());}
      inline NumericMatrix<T>& operator=(const NumericMatrix<T>& m2) {Matrix<T>::operator=(m2); return *this;}

      template<class U,class A=T> NumericMatrix<A> operator+(const NumericMatrix<U>& ext) const
      {
        if(ext.getWidth()!=this->getWidth() || ext.getHeight()!=this->getHeight()) throw MatrixDimensionMismatchException(make_pair(ext.getHeight(),ext.getWidth()),make_pair(this->h,this->w),"+");
        NumericMatrix<A> ret(this->w,this->h);
        for(int i=0;i<this->h;i++) {for(int j=0;j<this->w;j++) ret[i][j]=this->rows[i][j]+ext.getRows()[i][j];}
        return ret;
      }
      template<class U,class A=T> NumericMatrix<A> operator-(const NumericMatrix<U>& ext) const
      {
        if(ext.getWidth()!=this->w || ext.getHeight()!=this->h) throw MatrixDimensionMismatchException(make_pair(ext.getHeight(),ext.getWidth()),make_pair(this->h,this->w),"-");
        NumericMatrix<A> ret(this->w,this->h);
        for(int i=0;i<this->h;i++) {for(int j=0;j<this->w;j++) ret[i][j]=this->rows[i][j]-ext.getRows()[i][j];}
        return ret;
      }
      NumericMatrix<T> operator-() const
      {
        NumericMatrix<T> ret(this->w,this->h);
        for(int i=0;i<this->h;i++) for(int j=0;j<this->w;j++) ret[i][j]=-this->rows[i][j];
        return ret;
      }
      NumericMatrix<T> operator*(const T& el) const
      {
        NumericMatrix<T> ret(this->w,this->h);
        for(int i=0;i<this->h;i++) for(int j=0;j<this->w;j++) ret[i][j]=this->rows[i][j]*el;
        return ret;
      }
      NumericMatrix<T> operator/(const T& el) const
      {
        NumericMatrix<T> ret(this->w,this->h);
        for(int i=0;i<this->h;i++) for(int j=0;j<this->w;j++) ret[i][j]=this->rows[i][j]/el;
        return ret;
      }
      void operator*=(const T& el) {for(int i=0;i<this->h;i++) for(int j=0;j<this->w;j++) this->rows[i][j]=this->rows[i][j]*el;}
      void operator/=(const T& el) {for(int i=0;i<this->h;i++) for(int j=0;j<this->w;j++) this->rows[i][j]=this->rows[i][j]/el;}
      template<class U,class A=T> NumericMatrix<A> operator*(const NumericMatrix<U>& ext) const
      {
        if(this->getWidth()!=ext.getHeight()) throw MatrixDimensionMismatchException(make_pair(ext.getWidth(),ext.getHeight()),make_pair(this->h,this->w),"*");
        NumericMatrix<A> ret(ext.getWidth(),this->h);
        for(int i=0;i<this->h;i++)
        {
          for(int j=0;j<ext.getWidth();j++)
          {
            ret[i][j]=A();
            for(int k=0;k<this->w;k++) ret[i][j]=ret[i][j]+this->rows[i][k]*ext.getRows()[k][j];
          }
        }
        return ret;
      }
      template<class U=T> T average(const U& zero=U()) const
      {
        T mean=T(zero);
        for(int i=0;i<this->h;i++)
        {
          for(int j=0;j<this->w;j++)
            mean=mean+this->rows[i][j];
        }
        return mean/this->getSize();
      }
    private:
      T determinant(std::vector<int> omit) const //Needs optimization
      {
        if(this->getWidth()!=this->getHeight()) throw BadMatrixDimensionsForOperationException(make_pair(this->h,this->w),"det");
        T det=T();
        std::vector<int> allowed; for(int i=0;i<this->w;i++) {if(contains(omit,i)) continue; allowed.push_back(i);}
        if(omit.size()+1==this->w) return this->rows[this->h-1][allowed[0]];
        for(int j=0;j<allowed.size();j++)
        {
          std::vector<int> temp=omit; temp.push_back(allowed[j]);
          if(j%2!=0) det+=this->rows[omit.size()][allowed[j]]*determinant(temp);
          else det-=this->rows[omit.size()][allowed[j]]*determinant(temp);
        }
        return det;
      }
    public:
      inline T determinant() const {return determinant(std::vector<int>());}
    };
    template<class T=double> class Vector : public NumericMatrix<T>
    {
    public:
      Vector() : Vector(0) {}
      template<class U=T> Vector(int s,const U& zero=U()) : NumericMatrix<T>(1,s,zero) {}
      template<class U> Vector(const std::vector<U>& v,const U& zero=U()) : NumericMatrix<T>(1,v.size(),zero) {for(int i=0;i<v.size();i++) this->rows[i][0]=T(v[i]);}
      template<class U> Vector(const MatrixRow<U>& v,const U& zero=U()) : NumericMatrix<T>(1,v.getSize(),zero) {for(int i=0;i<v.getSize();i++) this->rows[i][0]=T(v[i]);}
      template<class U> Vector(const U* arr,int s,const U& zero=U()) : Vector(s,zero) {for(int i=0;i<arr.size();i++) this->rows[i][0]=T(arr[i]);}
      template<class U> Vector(const Matrix<U>& mat) : NumericMatrix<T>(mat)
      {
        if(mat.getWidth()!=1)
        {
          if(mat.getHeight()!=1) throw BadMatrixDimensionsForOperationException(make_pair(mat.getHeight(),mat.getWidth()),"to-vector");
          else {cout << "WARN: Transposing matrix to suit Vector\n"; operator=(this->transpose());}
        }
      }
      inline T getMagnitude() const {return norm();}
      inline T norm2() const
      {
        T ret=T();
        for(int i=0;i<this->rows.size();i++) ret=ret+this->rows[i][0]*this->rows[i][0];
        return ret;
      }
      inline T norm() const {return sqrt(norm2());}
      inline void normalize()
      {
        //cout << "Norm: " << norm() << "\n";1
        this->operator/=(norm());
      }
      template<class U> U dot(const Vector<U>& v2) const
      {
        if(v2.getSize()!=this->getSize()) throw MatrixDimensionMismatchException(make_pair(this->h,this->w),make_pair(v2.getHeight(),v2.getWidth()),"dot-product");
        U ret=this->rows[0][0]*v2(0);
        for(int i=1;i<this->rows.size();i++) ret=ret+this->rows[i][0]*v2(i);
        return ret;
      }
      template<class U> U dotIgnore(const Vector<U>& v2,const U& zero=U()) const
      {
        if(!this->size() || !v2.size()) return zero;
        U ret=this->rows[0][0]*v2(0);
        for(int i=1;i<min(this->rows.size(),v2.size());i++) ret=ret+this->rows[i][0]*v2(i);
        return ret;
      }

      T& operator()(int n)
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=this->rows.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-this->rows.size(),this->rows.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,this->rows.size());
  			#endif
        if(n<this->rows.size()) return this->rows[n][0];
        else throw commons::IndexOutOfBoundsException(n,this->rows.size());
      }
      const T& operator()(int n) const
      {
        #ifdef REVERSEPARSE
  				if(n<0) n=this->rows.size()+n;
  				if(n<0) throw commons::InvalidArrayIndexException(n-this->rows.size(),this->rows.size());
  			#else
  				if(n<0) throw commons::InvalidArrayIndexException(n,this->rows.size());
  			#endif
        if(n<this->rows.size()) return this->rows[n][0];
        else throw commons::IndexOutOfBoundsException(n,this->rows.size());
      }
      operator std::vector<T>() const
      {
        std::vector<T> ret; for(int i=0;i<this->rows.size();i++) ret.push_back(this->rows[i][0]);
        return ret;
      }
    };
  }
  namespace stats
  {
    template<class T> class Sample
    {
      std::vector<T> data;
    public:
      Sample() {}
      Sample(int s) {data=std::vector<T>(s);}
      Sample(const std::vector<T>& d) {data=d;}
      template<class U> Sample(int s,const U& zero) {for(int i=0;i<s;i++) data.push_back(T(zero));}

      const std::vector<T>& getData() const {return data;}
      inline void push_back(const T& obj) {data.push_back(obj);}
      inline int getSize() const {return data.size();}
      inline int size() const {return data.size();}

      inline auto begin() {return data.begin();}
      inline auto end() {return data.end();}

      template<class U=T> T mean(const U& zero=U()) const
      {
        T ret=T(zero);
        for(int i=0;i<data.size();i++) ret=ret+data[i];
        //cout << "Mean: " << ret<< " with size: "<<data.size()<<"\n";
        return ret/data.size();
      }
      template<class U=T> T variance(const T& avg=mean(),const U& zero=U()) const
      {
        T ret=T(zero);
        //T avg=mean();
        for(int i=0;i<data.size();i++) ret=ret+(data[i]-avg)*(data[i]-avg);
        return ret/data.size();
      }

      std::vector<T> bootstrapSample(int s=50) const
      {
        std::vector<T> ret;
        for(int i=0;i<s;i++) ret.push_back(data[randint(data.size())]);
        return ret;
      }

      T& operator[](int x) {return data[x];}
      const T& operator[](int x) const {return data[x];}
    };
  }
}
template<class T> std::ostream& operator<<(std::ostream& os,const general::matrices::MatrixRow<T>& mr)
{
  if(!mr.size()) {os << "||"; return os;}
  os << '|' << mr[0];
  for(int i=1;i<mr.size();i++) os <<"\t"<< mr[i];
  os << "|";
  return os;
}
template<class T> std::ostream& operator<<(std::ostream& os,const general::matrices::Matrix<T>& mr)
{
  if(!mr.getRowCount()) return os;
  else os << mr[0];
  for(int i=1;i<mr.getRowCount();i++) os <<"\n" << mr[i];
  return os;
}
//Adding operations to general::matrices

//Completing general::matrices
template<class T> template<class U> general::matrices::MatrixRow<T>::MatrixRow(const Vector<U>& v)
{
  for(int i=0;i<v.size();i++) elems.push_back(v(i));
}
template<class T> static std::ostream& operator<<(std::ostream& os,const general::stats::Sample<T>& s)
{
  os << "{";
  for(int i=0;i<s.size();i++) os << s[i] <<", ";
  return os <<"}";
}

template<class T> using GeneralMatrix = typename general::matrices::Matrix<T>;
template<class T> using GeneralMatrixRow = typename general::matrices::MatrixRow<T>;
template<class T=double> using NumericMatrix = typename general::matrices::NumericMatrix<T>;
template<class T=double> using GeneralVector = typename general::matrices::Vector<T>;
template<class T=double> using GeneralSample = typename general::stats::Sample<T>;
#endif
