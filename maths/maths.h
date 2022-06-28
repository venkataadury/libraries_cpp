#ifndef INCLUDED_MATHS
#define INCLUDED_MATHS 1
#define LIMIT 200
#include "commons/commons.h"
#include <assert.h>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "Eigen337/Eigen/Core"
#include "Eigen337/Eigen/Dense"
#include "Eigen337/Eigen/Geometry"

using namespace std;
namespace maths
{
	static int PRECISION=3;
	static bool radians=true;
	//Simple functions
	double pow(double n,int p)
	{
		if(p<0)
			return 1/pow(n,-p);
		if(p==0)
			return 1;
		return n*pow(n,p-1);
	}
	inline int sumFirst(int n=1,double p=1) {if(n==1) return 1; else return sumFirst(n-1,p)+pow(n,p);}
	inline double sqr(double n) {return pow(n,2);}
	inline short int sgn(double x) {return (x>0)?1:((x<0)?-1:0);}
	inline int resetRandSeed(unsigned long t=time(NULL)) {srand(t); return t;}
	inline double abs(double n) {return (n>=0)?n:-n;}
	inline int max(int n1,int n2) {return ((n1>n2)?n1:n2);} //returns maximum of two numbers
	inline int min(int n1,int n2) {return ((n1<n2)?n1:n2);} //returns minimum of two numbers
	inline double inf() {return std::numeric_limits<double>::infinity();}
	inline double nan() {return std::numeric_limits<double>::quiet_NaN();}
	inline double random() {return rand()/(double)RAND_MAX;}
	inline double randomAround(double X=0) {return random()-(X+0.5);}
	unsigned long int randint(unsigned long int l)
	{
		unsigned long int i=rand();
		while(i<l)
			i=i*RAND_MAX+rand();
		return i%l;
	}
	inline long int randint(int ll,int ul) {return randint(ul-ll)+ll;}
	inline bool toss(double p) {return (random() <= p);}
	inline bool toss() {return randint(2);}
	inline double round(double x,int err) {return (::round(x*pow(10,err)))/(double)pow(10,err);}



	inline bool close(double x1,double x2,int err=3) {return abs(x1-x2)<pow(10,-err);}
	inline void setPrecision(int n) {PRECISION=n; cout.precision(n);} //PRECISION in decimal (upto 'n' decimal places)
	//Angles -Start
	inline bool isRadians() {return radians;}
	inline bool isDegrees() {return !radians;}
	inline bool setRadians(bool v) {return (radians=v);}
	inline bool setDegrees(bool v) {return (radians=(!v));}
	inline double toRadians(double deg) {return (deg/180.0)*M_PI;}
	inline double toDegrees(double rad) {return (rad/M_PI)*180.0;}

	double sininv(double val)
	{
		double theta=asin(val);
		if(!radians)
			theta=toDegrees(theta);
		return theta;
	}
	double cosinv(double val)
	{
		double theta=acos(val);
		if(!radians)
			theta=toDegrees(theta);
		return theta;
	}
	double taninv(double val)
	{
		double theta=atan(val);
		if(!radians)
			theta=toDegrees(theta);
		return theta;
	}
	//Angles -End

	bool EQUALITY(int n1,int n2) {return n1==n2;}
	int gcd(int n1,int n2) //Finds greatest common divisor (GCD) of two numbers (positive integers only)
	{
		assert((n1>=0 && n2>=0));
		if(n1==0)
		{
			if(n2==0)
				return 1;
			else
				return n2;
		}
		if(n2==0)
			return n1;
		for(int i=min(n1,n2);i>1;i--)
		{
			if(n1%i==0 && n2%i==0)
				return i;
		}
		return 1;
	}
	int* range(int st,int end,int step=1)
	{
		if(end<st && step>0)
			step*=-1;
		int sz=(end-st)/step+1;
		int* ret=new int[sz];
		int K=0;
		if(end>=st)
		{
			for(int i=st;i<end;i+=step)
				ret[K++]=i;
		}
		else
		{
			for(int i=st;i>end;i+=step)
				ret[K++]=i;
		}
		return ret;
	}
	std::string toBinary(unsigned long int n)
	{
		std::string ret("");
		unsigned long int el,oper;
		for(int i=sizeof(n)*8-1;i>=0;i--)
		{
			oper=(1ul<<i);
			el=((n&oper)!=0);
			ret=(ret<<(int)el);
		}
		return ret;
	}

	int fact(int n)
	{
		if(n<=1)
			return 1;
		return n*fact(n-1);
	}

	inline double differentiate(double (*fx)(double),double x,double h=1e-6)
	{
		return (fx(x+h)-fx(x))/h;
	}
	template<class T,class U=T>
	U integrate(U (*fx)(T),T ll,T ul,T acc=0.001)
	{
		if(abs((double)(fx(ul)-fx(ll)))<acc)
			return fx((ul+ll)/2)*(ul-ll);
		return integrate(fx,(T)ll,(T)(ll+ul)/2,(T)acc)+integrate(fx,(T)(ll+ul)/2,(T)ul,(T)acc);
	}
	inline bool closeValues(double v1,double v2,double tol=0.5) {return abs(v2-v1)<=tol;}

	//pre-defined classes
	class MatrixRow;
	class Matrix;
	class Polynomial;
	class PolyMatrix;
	class Vector;
	//pre-defined functions
	static std::ostream& operator<<(std::ostream& os,Matrix mtr);
	static std::ostream& operator<<(std::ostream& os,Polynomial p);
	static std::ostream& operator<<(std::ostream& os,PolyMatrix pm);

	//Exceptions
	template<class T>class GeneralIndexOutOfBoundsException : public exception
	{
	public:
		GeneralIndexOutOfBoundsException()
		{
			cout <<"ERROR: Incompatible Index out of bounds operation executed"<<"\n";
		}
		GeneralIndexOutOfBoundsException(int s)
		{
			cout << "ERROR: Index "<< s <<" out of bounds\n";
		}
		GeneralIndexOutOfBoundsException(int i,int s)
		{
			cout << "ERROR: Index "<<i<<" out of bounds for size: "<<s<<"\n";
		}
	};
	class IncompatibleMatrixOperationException : public exception
	{
	public:
		IncompatibleMatrixOperationException()
		{
			cerr <<"ERROR: Incompatible Matrix operation executed"<<"\n";
		}
		IncompatibleMatrixOperationException(std::string msg) : IncompatibleMatrixOperationException() {cerr << msg<<"\n";}
	};
	class InvalidVectorOperationException : public exception
	{
	public:
		virtual const char* what() throw() {return "Invalid vector operation executed.";}
	};
	class InvalidVectorTransformException : public InvalidVectorOperationException
	{
	public:
		InvalidVectorTransformException() {cout << "Invalid transformation on vector\n";}
		InvalidVectorTransformException(commons::String str) : InvalidVectorTransformException() {cout << str.toString() << "\n";}
	};

	/*** Pure Maths ***/
	class Interval
	{
	public:
		double begin=0,end=0;
	public:
		Interval() {}
		Interval(double v) {begin=v;end=inf();}
		Interval(double s,double e)
		{
			if(s>e)
			{
				end=s;
				begin=e;
			}
			else
			{
				end=e;
				begin=s;
			}
		}
		//~Interval(): NOT Required
		inline double getLowerBound() {return begin;}
		inline double getUpperBound() {return end;}
		inline double span() {return end-begin;}
		inline bool contains(double v) const {return (v>=begin && v<=end);}
		inline bool exclusiveContains(double v) {return (v>begin && v<end);}
		inline bool inclusiveContains(double v) {return contains(v);}

		//Operator overloading
		inline bool operator==(Interval i2) {return (begin==i2.begin && end==i2.end);}
		bool operator>(Interval i2) {return (end>i2.end);}
		bool operator<(Interval i2) {return (begin<i2.begin);}
		bool operator>=(Interval i2) {return operator>(i2) || operator==(i2);}
		bool operator<=(Interval i2) {return operator<(i2) || operator==(i2);}
		inline bool operator!=(Interval i2) {return !operator==(i2);}
	};
	class BinaryNumber
	{
		long int value=0;
	public:
		BinaryNumber() {}
		BinaryNumber(long int n) {value=n;}
		BinaryNumber(std::string number) {value=std::stol(number);}

		inline long int getValue() {return value;}
		static BinaryNumber binaryString(std::string str)
		{
			long int n=0;
			for(int i=0;i<str.length();i++)
			{
				n<<=1;
				if(str.at(i)=='1')
					n|=1;
			}
			return BinaryNumber(n);
		}
		static BinaryNumber binaryString32(std::string str)
		{
			int m=1;
			char ch=str.at(0);
			str=str.substr(1);
			if(ch=='1')
				m=-1;
			int n=0;
			for(int i=0;i<str.length();i++)
			{
				n<<=1;
				if(str.at(i)=='1')
					n|=1;
			}
			return BinaryNumber(n*m);
		}

		inline std::string getBinary() {return toBinary(value);}
	};
	/*** Matrices ***/
	//Classes
	class MatrixRow
	{
		double* elems;
		int size;

	public:
		MatrixRow() {elems=new double[0]; size=0;}
		MatrixRow(int l)
		{
			size=l;
			elems=new double[l];
			std::fill(elems,elems+l,0);
		}
		MatrixRow(double*& els,int n)
		{
			size=n;
			elems=new double[n];
			for(int i=0;i<n;i++)
				elems[i]=els[i];
		}
		MatrixRow(std::string str) : MatrixRow(commons::String(str),' ') {}
		MatrixRow(commons::String str) : MatrixRow(str,' ') {}
		MatrixRow(std::string str,char ch) : MatrixRow(commons::String(str),ch) {}
		MatrixRow(const MatrixRow& ext) //Copy Constructor
		{
			size=ext.size;
			elems=new double[size];
			for(int i=0;i<size;i++)
				elems[i]=ext.elems[i];
		}
		MatrixRow(commons::String str,char ch)
		{
			int nC=str.countChar(ch);
			size=nC+1;
			elems=new double[size];
			commons::String* rws;
			rws=commons::String::cut(str,ch);
			std::string prc;
			for(int i=0;i<size;i++)
			{
				prc=rws[i].toString();
				elems[i]=std::stod(prc);
			}
		}
		~MatrixRow()
		{
			delete[] elems;
		}
		bool isZero() const
		{
			for(int i=0;i<size;i++)
			{
				if(elems[i]!=0)
					return false;
			}
			return true;
		}
		inline int getSize() const {return size;}
		double set(int n,double v)
		{
			#ifdef REVERSEPARSE
				if(n<0) n=size+n;
				if(n<0) throw commons::InvalidArrayIndexException(n-size,size);
			#else
				if(n<0) throw commons::InvalidArrayIndexException(n,size);
			#endif
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n]=v;
		}
		double get(int n) const
		{
			#ifdef REVERSEPARSE
				if(n<0) n=size+n;
				if(n<0) throw commons::InvalidArrayIndexException(n-size,size);
			#else
				if(n<0) throw commons::InvalidArrayIndexException(n,size);
			#endif
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n];
		}
		//Simple functions
		MatrixRow excludeElement(int e)
		{
			MatrixRow ret(size-1);
			for(int i=0;i<e;i++)
				ret[i]=elems[i];
			for(int i=e+1;i<size;i++)
				ret[i-1]=elems[i];
			return ret;
		}
		//operator overloading
		double& operator[](int n)
		{
			#ifdef REVERSEPARSE
				if(n<0) n=size+n;
				if(n<0) throw commons::InvalidArrayIndexException(n-size,size);
			#else
				if(n<0) throw commons::InvalidArrayIndexException(n,size);
			#endif
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n];
		}
		const double& operator[](int n) const
		{
			#ifdef REVERSEPARSE
				if(n<0) n=size+n;
				if(n<0) throw commons::InvalidArrayIndexException(n-size,size);
			#else
				if(n<0) throw commons::InvalidArrayIndexException(n,size);
			#endif
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n];
		}
		bool operator==(const MatrixRow& mr) const
		{
			if(size!=mr.size)
				return false;
			for(int i=0;i<size;i++)
			{
				if(elems[i]!=mr.elems[i])
					return false;
			}
			return true;
		}
		bool operator!=(const MatrixRow& mr) const {return !(this->operator==(mr));}
		MatrixRow operator+(const MatrixRow& mr) const
		{
			if(size!=mr.size)
				throw IncompatibleMatrixOperationException("No of colums in one row not same ");
			MatrixRow ret(size);
			for(int i=0;i<size;i++)
				ret[i]=mr[i]+elems[i];
			return ret;
		}
		void operator+=(const MatrixRow& mr)
		{
			if(size!=mr.size)
				throw IncompatibleMatrixOperationException("No of colums in one of the rows not same ");
			for(int i=0;i<size;i++)
				elems[i]+=mr.elems[i];
		}
		MatrixRow operator*(double val) const
		{
			MatrixRow ret(*this);
			for(int i=0;i<ret.getSize();i++)
				ret[i]*=val;
			return ret;
		}
		MatrixRow& operator=(const MatrixRow& row)
		{
			size=row.size;
			delete[] elems;
			elems=new double[size];
			for(int i=0;i<size;i++)
				elems[i]=row.elems[i];
			return *this;
		}
	};
	class Matrix
	{
	protected:
		MatrixRow* mrows;
		int hei=0;
	public:
		Matrix() : Matrix(0,0) {}
		Matrix(int w,int h)
		{
			mrows=new MatrixRow[h];
			for(int i=0;i<h;i++)
				mrows[i]=MatrixRow(w);
			//std::fill(mrows,mrows+h,*(new MatrixRow(w)));
			hei=h;
		}
		Matrix(double** matr,int rows,int cols)
		{
			mrows=new MatrixRow[rows];
			for(int i=0;i<rows;i++)
			{
				for(int j=0;j<cols;j++)
					mrows[i].set(j,matr[i][j]);
			}
			hei=rows;
		}
		Matrix(const Matrix& ext) //Copy Constructor
		{
			hei=ext.hei;
			mrows=new MatrixRow[hei];
			for(int i=0;i<hei;i++)
				mrows[i]=MatrixRow(ext.mrows[i]);
		}
		/*template<int r,int c> Matrix(double[r][c] mat)
		{
			hei=r;
			mrows=new MatrixRow[r];
			for(int i=0;i<hei;i++)
			{
				mrows[i]=MatrixRow(c);
				for(int j=0;j<c;j++) mrows[i][j]=mat[i][j];
			}
		}*/
		~Matrix()  {delete[] mrows;}
	protected:
		Matrix(int rws)
		{
			hei=rws;
			mrows=new MatrixRow[rws];
			for(int i=0;i<hei;i++)
				mrows[i]=MatrixRow();
		}
	public:
		//Getting data from matrix
		 inline int getRows() const {return hei;}
		 int getColumns() const
		 {
			 if(hei==0)
				 return 0;
			 return mrows[0].getSize();
		}
		inline int getSize() const {return getRows()*getColumns();}
		double set(int r,int c,double v)
		{
			if(r>=hei)
				throw commons::IndexOutOfBoundsException(r,hei);
			return mrows[r].set(c,v);
		}
		double get(int r,int c) const
		{
			if(r>=hei)
				throw commons::IndexOutOfBoundsException(r,hei);
			return mrows[r].get(c);
		}
		const MatrixRow& getRow(int i) const
		{
			if(i>=hei)
				throw commons::IndexOutOfBoundsException(i,hei);
			return mrows[i];
		}
		Vector getColumn(int n) const;
		//Checking for uniformity in shape of matrix
		bool isUniform() const//All rows have same no of columns?
		{
			int r=getColumns();
			for(int i=0;i<hei;i++)
			{
				if(mrows[i].getSize()!=r)
					return false;
			}
			return true;
		}
		inline bool isSquare() const {return (isUniform() && getRows()==getColumns());}
		//Building up
		void addRow(MatrixRow ro)
		{
			if(hei==0)
			{
				mrows=new MatrixRow[1] {ro};
				hei=1;
				return;
			}
			MatrixRow* matr=new MatrixRow[hei+1];
			std::copy(mrows,mrows+hei,matr);
			delete[] mrows;
			matr[hei]=ro;
			mrows=matr;
			hei++;
		}
		Matrix excludeRow(int r)
		{
			Matrix ret;
			for(int i=0;i<r;i++)
				ret.addRow(mrows[i]);
			for(int i=r+1;i<hei;i++)
				ret.addRow(mrows[i]);
			return ret;
		}
		Matrix excludeColumn(int c)
		{
			Matrix ret;
			for(int i=0;i<hei;i++)
				ret.addRow(mrows[i].excludeElement(c));
			return ret;
		}
		//Matrix operations
		double determinant()
		{
			if(getRows()!=getColumns())
				throw IncompatibleMatrixOperationException("Determinant of non-square matrix requested.");
			if(getRows()==0)
				return 1;

			int cC=getColumns();
			if(cC<=1)
				return mrows[0][0];
			double det=0;
			for(int i=0;i<cC;i++)
			{
				Matrix temp=excludeRow(0).excludeColumn(i);
				det+=mrows[0][i]*temp.determinant()*((i%2==0)?1:-1);
			}
			return det;
		}
		Matrix getCofactor()
		{
			Matrix ret(getRows(),getColumns());
			for(int i=0;i<hei;i++)
			{
				for(int j=0;j<mrows[i].getSize();j++)
					ret[i][j]=this->excludeRow(i).excludeColumn(j).determinant()*(((i+j)%2==0)?1:-1);
			}
			return ret;
		}
		Matrix transpose()
		{
			Matrix ret(getRows(),getColumns()); //Width-Height Initialization
			for(int i=0;i<ret.hei;i++)
			{
				for(int j=0;j<ret.mrows[i].getSize();j++)
					ret[i][j]=this->get(j,i);
			}
			return ret;
		}
		inline Matrix getAdjoint() {return getCofactor().transpose();}
		Matrix getInverse()
		{
			double det=determinant();
			if(det==0)
				throw IncompatibleMatrixOperationException("Requested inverse of a singular matrix.");
			Matrix ret=getAdjoint()/det;
			return ret;
		}
		Matrix squeeze()
		{
			Matrix ret(1,getRows()+getColumns()-1);
			for(int k=0;k<ret.hei;k++)
			{
				for(int i=0;i<hei;i++)
				{
					for(int j=0;j<mrows[i].getSize();j++)
					{
						if((i+j)==k)
							ret[k][0]+=get(i,j);
					}
				}
			}
			return ret;
		}

		Polynomial getCharDeterminant();
		std::vector<double> getEigenvalues();
		Matrix* diagonalize();
		//static methods
		static Matrix getIdentity(int s)
		{
			Matrix ret(s,s);
			for(int i=0;i<s;i++)
			{
				for(int j=0;j<s;j++)
					ret[i][j]=(i==j)?1:0;
			}
			return ret;
		}
		static Matrix getUnitMatrix(int r,int c,int w) {return getUnitMatrix(r,c,w,w);}
		static Matrix getUnitMatrix(int r,int c,int w,int h)
		{
			Matrix mat(w,h);
			mat[r][c]=1;
			return mat;
		}
		static Matrix getMatrix(double (*rule)(int,int),int w,int h) //Width-Height Instantiation
		{
			Matrix ret(w,h);
			for(int i=0;i<h;i++)
			{
				for(int j=0;j<w;j++)
					ret[i][j]=rule(i,j);
			}
			return ret;
		}
		//operator overloading
		MatrixRow& operator[](int n)
		{
			if(n>=hei)
				throw commons::IndexOutOfBoundsException(n,hei);
			return mrows[n];
		}
		const MatrixRow& operator[](int n) const
		{
			if(n>=hei)
				throw commons::IndexOutOfBoundsException(n,hei);
			return mrows[n];
		}
		bool operator==(const Matrix& m2) const
		{
			if(hei!=m2.hei)
				return false;
			for(int i=0;i<hei;i++)
			{
				if(mrows[i]!=m2.mrows[i])
					return false;
			}
			return true;
		}
		bool operator!=(const Matrix& m2) const {return !(this->operator==(m2));}
		Matrix operator+(const Matrix& m2) const
		{
			if(getRows()!=m2.getRows())
				throw IncompatibleMatrixOperationException("Number of rows don't match");
			Matrix ret(hei);
			for(int i=0;i<hei;i++)
				ret[i]=mrows[i]+m2.mrows[i];
			return ret;
		}
		void operator+=(const Matrix& m2)
		{
			for(int i=0;i<hei;i++)
				mrows[i]+=m2.mrows[i];
		}
		Matrix operator*(double val) const
		{
			Matrix m(*this);
			for(int i=0;i<m.getRows();i++)
			{
				for(int j=0;j<m[i].getSize();j++)
					m[i][j]*=val;
			}
			return m;
		}
		Matrix operator/(double val) const
		{
			Matrix m(*this);
			for(int i=0;i<m.getRows();i++)
			{
				for(int j=0;j<m[i].getSize();j++)
					m[i][j]/=val;
			}
			return m;
		}
		inline Matrix operator-(const Matrix& m2) const  {return this->operator+(m2*-1);}
		Matrix operator*(const Matrix& m2) const
		{
			if(this->getColumns()!=m2.getRows())
				throw IncompatibleMatrixOperationException("Incompatible multiplication");
			Matrix ret(m2.getColumns(),this->getRows());
			int rows=ret.getRows(),columns=ret.getColumns();
			//DEBUG: cout <<"\t"<<columns<<" "<<m2.getColumns()<<"\n";
			for(int i=0;i<rows;i++)
			{
				for(int j=0;j<columns;j++)
				{
					for(int k=0;k<this->getColumns();k++)
						ret[i][j]+=this->get(i,k)*m2.get(k,j);
				}
			}
			return ret;
		}
		Matrix operator-() {return operator*(-1);}
		Matrix& operator=(const Matrix& obj) //Assignment Operator
		{
			hei=0;
			hei=obj.hei;
			delete[] mrows;
			mrows=new MatrixRow[hei];
			for(int i=0;i<hei;i++)
				mrows[i]=obj.getRow(i);
			return *this;
		}
		std::vector<Vector> getEigenvectors();
		//static methods
		static Matrix getRotation(double theta,double phi)
		{
			/*
			 * z=rcosA
			 * x=rsinAcosB
			 * y=rsinAsinB
			 */
			if(!radians)
			{
				theta=toRadians(theta);
				phi=toRadians(phi);
			}
			Matrix rot=Matrix(3,3);
			rot[0][0]=cos(theta)*cos(phi); rot[0][1]=cos(theta)*sin(phi); rot[0][2]=-sin(theta);
			rot[1][0]=-sin(phi); rot[1][1]=cos(phi); rot[1][2]=0;
			rot[2][0]=sin(theta)*cos(phi); rot[2][1]=sin(theta)*sin(phi); rot[2][2]=cos(theta);
			return rot;
		}
	};
	class Vector : public Matrix
	{
	public:
		Vector() : Matrix() {}
		Vector(int sz) : Matrix(1,sz) {}
		explicit Vector(double* els,int sz) : Matrix(1,sz)
		{
			for(int i=0;i<hei;i++)
				mrows[i][0]=els[i];
		}
		Vector(double x,double y) : Vector(2)
		{
			mrows[0][0]=x;
			mrows[1][0]=y;
		}
		Vector(double x,double y,double z) : Vector(3)
		{
			mrows[0][0]=x;
			mrows[1][0]=y;
			mrows[2][0]=z;
		}
		Vector(std::vector<double> vct) : Matrix(1,vct.size())
		{
			for(int i=0;i<hei;i++)
				mrows[i][0]=vct[i];
		}
		Vector(Matrix m) : Matrix(m) //Refer Matrix::Matrix(const Matrix&) for working of constructor
		{
			if(m.getColumns()!=1)
				throw IncompatibleMatrixOperationException("Cannot convert multi-column matrix to Vector. Are you trying with single-row multi column? (Did you forget to transpose?)");
		}

		double getX() const
		{
			if(hei<=0)
				throw commons::IndexOutOfBoundsException(0,hei);
			return get(0,0);
		}
		double getY() const
		{
			if(hei<=1)
				return 0;
			return get(1,0);
		}
		double getZ() const
		{
			if(hei<=2)
				return 0;
			return get(2,0);
		}
		double getMagnitude() const
		{
			double mag=0;
			for(int i=0;i<hei;i++)
				mag+=pow(get(i,0),2);
			return sqrt(mag);
		}
		bool isNull() const {return getSize()==0;}
		inline double& at(int i) {return mrows[i][0];}
		inline double comp(int n) const {return mrows[n][0];}

		//Vector operations
		double dot(const Vector& v2) const
		{
			if(v2.getSize()!=getSize())
				throw InvalidVectorOperationException();
			double ret=0;
			for(int i=0;i<getSize();i++)
				ret+=get(i,0)*v2[i][0];
			return ret;
		}
		Vector cross(const Vector& v2) const // Only computes 3D cross
		{
			Vector ret(3);
			ret.at(0)=comp(1)*v2.comp(2)-comp(2)*v2.comp(1);
			ret.at(1)=-(comp(0)*v2.comp(2)-comp(2)*v2.comp(0));
			ret.at(2)=comp(0)*v2.comp(1)-comp(1)*v2.comp(0);
			return ret;
		}
		inline double angleWith(const Vector& v2) const {return cosinv(dot(v2)/(v2.getMagnitude()*getMagnitude()));}
		Vector unitVector() const {return (*this)/getMagnitude();}
		Vector forceMagnitude(double mag) const {return unitVector()*mag;}
		Vector subVector(int C) const
		{
			if(C>getSize())
				return *this;
			Vector ret(C);
			for(int i=0;i<ret.getRows();i++)
				ret[i][0]=(*this).comp(i);
			return ret;
		}
		bool isZero() const
		{
			for(int i=0;i<getSize();i++)
			{
				if(mrows[i][0]!=0)
					return false;
			}
			return true;
		}
		Vector quickTransform(Matrix mat) const {return mat*(*this);}
		Vector transform(Matrix mat) const
		{
			if(mat.getRows()>getSize())
				throw InvalidVectorTransformException("Incompatible matrix to multiply. Too large for vector.");
			return mat*subVector(mat.getRows());
		}
		Vector componentAlong(const Vector& v2)   {return v2.forceMagnitude(dot(v2)/v2.getMagnitude());}

		//Static methods
		static Matrix buildMatrix(std::vector<Vector> va)
		{
			if(va.size()==0)
				return Matrix(0,0);
			Matrix ret(va.size(),va[0].getSize());
			for(int i=0;i<va.size();i++)
			{
				for(int j=0;j<va[i].getSize();j++)
					ret[j][i]=va[i][j][0];
			}
			return ret;
		}
		//Operator overloading
		Vector operator+(Vector v2) const
		{
			int sz=max(v2.getSize(),getSize());
			Vector ret(sz);
			for(int i=0;i<sz;i++)
			{
				if(i<v2.getSize())
					ret[i][0]+=v2[i][0];
				if(i<getSize())
					ret[i][0]+=(*this)[i][0];
			}
			return ret;
		}
		bool operator==(const Vector& v) const
		{
			if(v.getSize()!=getSize())
				return false;
			for(int i=0;i<getSize();i++)
			{
				if(v[i][0]!=get(i,0))
					return false;
			}
			return true;
		}

		//Static Methods
		static Vector basis(int x,int i)
		{
			Vector ret(x);
			ret[i][0]=1;
			return ret;
		}
		/*Vector& operator=(const Vector& ext)
		{
			cout << "Vector =\n";
			dynamic_cast<Matrix*>(this)->operator=(ext);
			return *this;
		}*/
	};

	/*** Function Algebra***/
	class Polynomial
	{
		double* coeffs=new double[1] {0};
		int degree=0;
		char var='x';
		Polynomial* der;
	public:
		Polynomial() {}
		Polynomial(int deg)
		{
			degree=deg;
			coeffs=new double[deg+1];
			std::fill(coeffs,coeffs+deg+1,0);
		}
		Polynomial(double* cfs,int deg)
		{
			degree=deg;
			coeffs=cfs;
		}
		Polynomial(Vector v)
		{
			degree=v.getSize()-1;
			coeffs=new double[degree+1];
			for(int i=0;i<=degree;i++)
				coeffs[i]=v[i][0];
		}
		Polynomial(const Polynomial& poly)
		{
			degree=poly.degree;
			var=poly.var;
			der=poly.der;
			coeffs=new double[degree+1];
			for(int i=0;i<=degree;i++)
				coeffs[i]=poly.coeffs[i];
		}
		~Polynomial() {delete[] coeffs;}

		inline void setVar(char ch) {var=ch;}
		inline char getVar() const {return var;}
		double getVal(double x) const
		{
			double s=0;
			for(int i=0;i<=degree;i++)
				s+=pow(x,i)*coeffs[i];
			return s;
		}
		inline double* getCoeffs() const {return coeffs;}
		inline int getDegree() const {return degree;}
		bool operator==(const Polynomial& p2) const
		{
			if(degree!=p2.degree)
				return false;
			for(int i=0;i<=degree;i++)
			{
				if(coeffs[i]!=p2.coeffs[i])
					return false;
			}
			return true;
		}
		Polynomial derivative()
		{
			Polynomial ret(degree-1);
			ret.setVar(var);
			for(int i=0;i<=ret.degree;i++)
				ret.coeffs[i]=coeffs[i+1]*(i+1);
			return ret;
		}
		Polynomial integral(double c=0) const
		{
			Polynomial ret(degree+1);
			ret.setVar(var);
			for(int i=0;i<=degree;i++)
				ret.coeffs[i+1]=coeffs[i]/(i+1);
			ret.coeffs[0]=c;
			return ret;
		}
		inline double integrate(double ll,double ul) const
		{
			Polynomial p=integral();
			return p.getVal(ul)-p.getVal(ll);
		}

		//Solving
		double bisectionRoot(double ll,double ul,double err)
		{
			if(getVal(ul)*getVal(ll)>0)
				return ll-1;
			if(abs(getVal(ll))<err)
				return ll;
			if(abs(getVal(ul))<err)
				return ul;
			double m=(ul+ll)/2.0;
			if(abs(getVal(m))<err)
				return m;
			if(getVal(ul)*getVal(m)<0)
				return bisectionRoot(m,ul,err);
			else
				return bisectionRoot(ll,m,err);
		}
		double newtonraphsonRoot(double x,double err,double errVal=0,int iter=0)
		{
			if(iter==0)
			{
				Polynomial p=derivative();
				der=&p;
			}
			if(abs(getVal(x))<err)
				return x;
			double m=der->getVal(x);
			double c=getVal(x)-m*x;
			if(iter>LIMIT)
				return errVal;
			if(m==0)
				m=err;
			return newtonraphsonRoot(-c/m,err,errVal,++iter);
		}
		std::vector<double> getRoots(int eval=3)
		{
			double err=pow(10,-eval);
			if(degree==0)
				return std::vector<double>();
			if(degree==1)
				return std::vector<double>(1,-coeffs[0]/coeffs[1]);
			std::vector<double> dR=derivative().getRoots();
			std::vector<double> fR=sort(set(dR)); //Removes duplicates "set(vector)" and sorts
			int uR=fR.size();
			if(uR==0)
				return std::vector<double> {newtonraphsonRoot(0,err)};
			std::vector<double> ret; //=*(new std::vector<double>());
			double r;
			for(int i=0;i<fR.size()-1;i++)
			{
				r=bisectionRoot(fR[i],fR[i+1],err);
				if(r<fR[i])
					continue;
				ret.push_back(r);
			}
			double FR=fR[0],LR=fR[fR.size()-1];
			double r1=newtonraphsonRoot(LR+1,err,LR-1);
			double r2=newtonraphsonRoot(FR-1,err,FR+1);
			if(r1!=LR-1)
				ret.push_back(r1);
			if(r2!=FR+1)
				ret.push_back(r2);
			return ret;
		}
		inline Vector getVector() const {return Vector(coeffs,degree+1);}
		//Operator overloading
		inline bool operator!=(const Polynomial& p2) const {return !operator==(p2);}
		Polynomial operator+(const Polynomial& p2) const
		{
			int D=max(p2.degree,degree);
			Polynomial ret(D);
			ret.setVar(var);
			for(int i=0;i<=D;i++)
			{
				if(degree>=i)
					ret.coeffs[i]+=coeffs[i];
				if(p2.degree>=i)
					ret.coeffs[i]+=p2.coeffs[i];
			}
			return ret;
		}
		Polynomial operator*(double n) const
		{
			Polynomial P(degree);
			P.setVar(var);
			for(int i=0;i<=degree;i++)
				P.coeffs[i]=coeffs[i]*n;
			return P;
		}
		inline Polynomial operator-(const Polynomial& pn) const {return operator+(pn*-1);}
		inline Polynomial operator-() const {return operator*(-1);}
		Polynomial operator*(const Polynomial& p2) const
		{
			if(p2.degree==0)
				return operator*(p2.getVal(0));
			if(degree==0)
				return p2*((double)getVal(0));
			Vector v1=getVector();
			Matrix tp=p2.getVector().transpose();
			Matrix m=v1*tp;
			return *(new Polynomial((Vector)m.squeeze()));
		}
		void operator*=(double val)
		{
			for(int i=0;i<=degree;i++)
				coeffs[i]*=val;
		}
		Polynomial& operator=(const Polynomial& poly)
		{
			degree=poly.degree;
			var=poly.var;
			der=poly.der;
			coeffs=new double[degree+1];
			for(int i=0;i<=degree;i++)
				coeffs[i]=poly.coeffs[i];
			return *this;
		}
		//Static methods
		inline static Polynomial CharPolynomial(double val) {return Polynomial(new double[2] {val,-1},1);}
	};
	//Mixed Polynomials and Matrices
	class PolyMatrixRow
	{
		Polynomial* elems=new Polynomial[0];
		int size=0;

	public:
		PolyMatrixRow() {}
		PolyMatrixRow(int l)
		{
			size=l;
			elems=new Polynomial[l];
		}
		PolyMatrixRow(Polynomial*& els,int n)
		{
			size=n;
			elems=els;
		}

		inline int getSize() {return size;}
		Polynomial set(int n,Polynomial v)
		{
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n]=v;
		}
		Polynomial& get(int n)
		{
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n];
		}
		//Simple functions
		PolyMatrixRow excludeElement(int e)
		{
			PolyMatrixRow ret(size-1);
			for(int i=0;i<e;i++)
				ret[i]=elems[i];
			for(int i=e+1;i<size;i++)
				ret[i-1]=elems[i];
			return ret;
		}
		//operator overloading
		Polynomial& operator[](int n)
		{
			if(n>=size)
				throw commons::IndexOutOfBoundsException(n,size);
			return elems[n];
		}
		bool operator==(PolyMatrixRow mr)
		{
			if(size!=mr.size)
				return false;
			for(int i=0;i<size;i++)
			{
				if(elems[i]!=mr.elems[i])
					return false;
			}
			return true;
		}
		bool operator!=(PolyMatrixRow mr) {return !(this->operator==(mr));}
		PolyMatrixRow operator+(PolyMatrixRow mr)
		{
			if(size!=mr.size)
				throw IncompatibleMatrixOperationException("No of colums in one row not same ");
			PolyMatrixRow ret(size);
			for(int i=0;i<size;i++)
				ret[i]=mr[i]+elems[i];
			return ret;
		}
		PolyMatrixRow operator*(double val)
		{
			PolyMatrixRow ret(*this);
			for(int i=0;i<ret.getSize();i++)
				ret[i]*=val;
			return ret;
		}
	};
	class PolyMatrix
	{
		PolyMatrixRow* mrows=new PolyMatrixRow[0];
		int hei=0;
	public:
		PolyMatrix() {}
		PolyMatrix(int h)
		{
			hei=h;
			mrows=new PolyMatrixRow[h];
		}
		PolyMatrix(int w,int h) //Width-Height Instantiation
		{
			hei=h;
			mrows=new PolyMatrixRow[h];
			for(int i=0;i<h;i++)
				mrows[i]=PolyMatrixRow(w);
		}
		PolyMatrix(Matrix orig,bool (*cond)(int x,int y),Polynomial (*fx)(double val))
		{
			int wid=orig.getColumns();
			hei=orig.getRows();
			mrows=new PolyMatrixRow[hei];
			for(int i=0;i<hei;i++)
				mrows[i]=PolyMatrixRow(orig[i].getSize());
			for(int i=0;i<hei;i++)
			{
				for(int j=0;j<mrows[i].getSize();j++)
				{
					if(cond(i,j))
						mrows[i][j]=fx(orig[i][j]);
					else
						mrows[i][j]=Polynomial(new double[1] {orig[i][j]},0);
				}
			}
		}

		inline int getRows() {return hei;}
		inline int getColumns() {return ((hei==0)?0:mrows[0].getSize());}
		//Building up
		void addRow(PolyMatrixRow ro)
		{
			if(hei==0)
			{
				mrows=new PolyMatrixRow[1] {ro};
				hei=1;
				return;
			}
			PolyMatrixRow* matr=new PolyMatrixRow[hei+1];
			std::copy(mrows,mrows+hei,matr);
			matr[hei]=ro;
			mrows=matr;
			hei++;
		}
		PolyMatrix excludeRow(int r)
		{
			PolyMatrix ret;
			for(int i=0;i<r;i++)
				ret.addRow(mrows[i]);
			for(int i=r+1;i<hei;i++)
				ret.addRow(mrows[i]);
			return ret;
		}
		PolyMatrix excludeColumn(int c)
		{
			PolyMatrix ret;
			for(int i=0;i<hei;i++)
				ret.addRow(mrows[i].excludeElement(c));
			return ret;
		}
		//Matrix operations
		Polynomial determinant()
		{
			if(getRows()!=getColumns())
				throw IncompatibleMatrixOperationException("Determinant of non-square matrix requested.");
			int cC=getColumns();
			if(cC<=1)
				return mrows[0][0];
			Polynomial det(new double[1] {0},0);
			Polynomial temp;
			for(int i=0;i<cC;i++)
				det=det+((mrows[0][i]*((excludeRow(0).excludeColumn(i)).determinant()))*((i%2==0)?1:-1));
			return det;
		}
		//operator overloading
		PolyMatrixRow& operator[](int x)
		{
			if(x>=hei)
				throw commons::IndexOutOfBoundsException(x,hei);
			return mrows[x];
		}
		const PolyMatrixRow& operator[](int x) const
		{
			if(x>=hei)
				throw commons::IndexOutOfBoundsException(x,hei);
			return mrows[x];
		}
	};
	class EquationSystem
	{
		int nvar=0;
		Matrix eqns;
		Vector result;
		int K=0;
	public:
		EquationSystem() {}
		EquationSystem(int nv)
		{
			nvar=nv;
			eqns=Matrix(nv,nv);
			result=Vector(nv);
		}
		EquationSystem(Matrix es,Vector r)
		{
			nvar=r.getSize();
			result=r;
			eqns=es;
		}
		EquationSystem(Matrix es) : EquationSystem(es,Vector(new double[es.getRows()],es.getRows())) {} //TODO: Check

		inline Matrix getMatrixForm() {return eqns;}
		inline Vector getResults() {return result;}
		void addEquation(Vector coeffs,double res)
		{
			for(int i=0;i<coeffs.getSize();i++)
				eqns[K][i]=coeffs[i][0];
			result[K++][0]=res;
		}
		EquationSystem dropLastVariables()
		{
			EquationSystem ret(*this);
			while(K<ret.nvar)
			{
				ret.eqns=ret.eqns.excludeColumn(ret.nvar-1).excludeRow(ret.nvar-1);
				ret.result=ret.result.excludeRow(ret.nvar-1);
				ret.nvar--;
			}
			return ret;
		}
		EquationSystem reduce(int var,double val,int ie=-1)
		{
			if(ie==-1)
				ie=var;
			Matrix nes(nvar-1,nvar-1);
			Vector ner(nvar-1);
			int K=0;
			for(int i=0;i<eqns.getRows();i++)
			{
				if(i==ie)
					continue;
				for(int j=0;j<var;j++)
					nes[K][j]=eqns[i][j];
				for(int j=var+1;j<eqns[i].getSize();j++)
					nes[K][j-1]=eqns[i][j];
				ner[K][0]=result[i][0]-val*eqns[i][var];
				K++;
			}
			return EquationSystem(nes,ner);
		}
		inline bool solvable(int err=3)
		{
			Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic> rm(eqns.getRows(),eqns.getColumns());
			for(int i=0;i<eqns.getRows();i++) for(int j=0;j<eqns.getColumns();j++) rm(i,j)=eqns[i][j];
			return !close(rm.determinant(),0,err);
		}
		inline int getVarCount() {return nvar;}
		Vector solve()
		{
			Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic> rm(eqns.getRows(),eqns.getColumns());
			for(int i=0;i<eqns.getRows();i++) for(int j=0;j<eqns.getColumns();j++) rm(i,j)=eqns[i][j];
			//cout << rm << "\n";
			if(rm.determinant()==0)
			{
				//cout << rm << "\n";
				throw IncompatibleMatrixOperationException("Cannot solve equation system. Will have to invert singular matrix! Try reducing the system");
			}
			Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic> inv=rm.inverse();
			//cout << inv << "\n";
			Matrix im(eqns.getColumns(),eqns.getRows()); for(int i=0;i<eqns.getColumns();i++) for(int j=0;j<eqns.getRows();j++) im[i][j]=inv(i,j);
			//cout << im << "\n";
			Vector r=(Vector)(im*result);
			//cout << r << "\n";
			return r;
		}
	};

	//Statistics and Data
	class Histogram
	{
		Interval* intervals=new Interval[0];
		int* count=new int[0];
		int nI=0;
		double start=0,interval_size=1;
	public:
		Histogram() {}
		Histogram(double st) {start=st;}
		Histogram(double st,double intv) {start=st; interval_size=intv;}

		void add(double val)
		{
			bool flag=false;
			for(int i=0;i<nI;i++)
			{
				if(intervals[i].contains(val))
				{
					flag=true;
					count[i]++;
					return;
				}
			}
			if(!flag)
			{
				/*
				 * [0,0.5] 2.6
				 * 0 + (2.6-0)
				 * [1.8,0.5] 2.2
				 * 1.8 + (2.2-1.8)
				 * interval_size*? = just<v (Question)
				 * [v/interval_size] <Greatest Integer Function>
				 */
				double iS=(int)((val-start)/interval_size);
				iS*=(double)interval_size;
				iS+=start;
				Interval newI(iS,iS+interval_size);
				intervals=append(intervals,newI,nI);
				count=append(count,1,nI++);
			}
		}
		inline double getStart() {return start;}
		inline double getStep() {return interval_size;}
		inline int getIntervalCount() {return nI;}
		Interval& getIntervalFor(double v)
		{
			for(int i=0;i<nI;i++)
			{
				if(intervals[i].contains(v))
					return intervals[i];
			}
		}
		Interval getIntervalAt(int i)
		{
			if(i>=nI)
				throw commons::IndexOutOfBoundsException(i,nI);
			return intervals[i];
		}
		int getCountAt(int i)
		{
			if(i>=nI)
				throw commons::IndexOutOfBoundsException(i,nI);
			return count[i];
		}
		int getCount(Interval I)
		{
			for(int i=0;i<nI;i++)
			{
				if(intervals[i]==I)
					return count[i];
			}
			return -1;
		}
		inline Interval* getIntervals() {return intervals;}

	};
	struct Point2D
	{
		std::pair<double,double> pt;
	public:
		Point2D() : Point2D(0,0) {}
		Point2D(double x) : Point2D(x,0) {}
		Point2D(double x,double y) {pt=make_pair(x,y);}
		Point2D(const std::pair<double,double>& p) {pt=p;}
		inline double getX() const {return get<0>(pt);}
		inline double getY() const {return get<1>(pt);}
		inline void setX(double x) {pt=make_pair(x,getY());}
		inline void setY(double y) {pt=make_pair(getX(),y);}
	};
	class PolynomialGenerator
	{
		std::vector<Point2D> points;
		EquationSystem eqsys;
	public:
		PolynomialGenerator() {}
		PolynomialGenerator(const std::vector<Point2D>& p) {points=p;}
		PolynomialGenerator(const std::vector<std::pair<double,double>>& p) {points=std::vector<Point2D>(); for(int i=0;i<p.size();i++) points.push_back(p[i]);}

		inline const Point2D& getPoint(int x) const {return points[x];}
		inline Point2D& getPoint(int x) {return points[x];}

		EquationSystem& buildEquations()
		{
			int ord=points.size();
			eqsys=EquationSystem(ord);
			for(int i=0;i<ord;i++)
			{
				double x=points[i].getX(),y=points[i].getY();
				std::vector<double> v; for(int j=0;j<ord;j++) v.push_back(pow(x,j));
				//cout << "{"; for(double& d : v) cout << d<<","; cout << "}\t"<<y<<"\n";
				eqsys.addEquation(v,y);
			}
			return eqsys;
		}
		Polynomial getPolynomial()
		{
			if(!eqsys.getVarCount())  buildEquations();
			return Polynomial(eqsys.solve());
		}
	};

	//Basis for complex development
	/*template<class T> class GeneralMatrixRow
	{
		T* elems;
		int len=0;
	public:
		GeneralMatrixRow() : GeneralMatrixRow(0) {}
		GeneralMatrixRow(int s) {len=s; elems=new T[len];}
		GeneralMatrixRow(T*& els,int s)
		{
			elems=new T[s];
			len=s;
			for(int i=0;i<s;i++)
				elems[i]=els[i];
		}

		inline int size() {return len;}
		inline int getSize() {return size();}
		inline int length() {return size();}

		//operator overloading
		T& operator[](int n)
		{
			if(n>=len)
				throw commons::IndexOutOfBoundsException(n,len);
			return elems[n];
		}
		const T& operator[](int n) const
		{
			if(n>=len)
				throw commons::IndexOutOfBoundsException(n,len);
			return elems[n];
		}
	};*/
	namespace fourier
	{
		inline static void init() {radians=true;}
		template<int T> class FourierSeries
		{
			std::vector<double> sins,coss;
			double k;
		public:
			FourierSeries()
			{
				sins=std::vector<double>(T); coss=std::vector<double>(T); k=0;
				for(int i=0;i<T;i++) sins[i]=coss[i]=0;
			}
			FourierSeries(const std::vector<double>& s,const std::vector<double>& c,double a=0)
			{
				k=a;
				sins=std::vector<double>(T); coss=std::vector<double>(T);
				for(int i=0;i<min(T,s.size());i++) sins[i]=s[i];
				for(int i=0;i<min(T,c.size());i++) coss[i]=c[i];
				for(int i=s.size();i<T;i++) sins[i]=0;
				for(int i=c.size();i<T;i++) coss[i]=0;
			}
			template<int U> FourierSeries(const FourierSeries<U>& fs)
			{
				sins=std::vector<double>(T); coss=std::vector<double>(T); k=fs.getConstant();
				for(int i=0;i<min(U,T);i++)
				{
					sins[i]=fs.getSin(i);
					coss[i]=fs.getCos(i);
				}
				for(int i=U;i<T;i++) sins[i]=coss[i]=0;
			}

			inline double getConstant() const {return k;}
			inline double getSin(int x) const {return sins[x];}
			inline double getCos(int x) const {return coss[x];}
		};
		class FourierGenerator
		{
			std::vector<Point2D> points;
		public:
			FourierGenerator() {}
			FourierGenerator(const std::vector<Point2D>& p) {points=p; cout << "Constructed generator for fouriers\n";}
			FourierGenerator(const std::vector<std::pair<double,double>>& p) {points=std::vector<Point2D>(); for(int i=0;i<p.size();i++) points.push_back(p[i]);}

			inline const Point2D& getPoint(int x) const {return points[x];}
			inline Point2D& getPoint(int x) {return points[x];}

			template<int T> FourierSeries<T> getFourierSeries() const
			{
				assert(T*2<points.size());
				cout << "Get fourier series: <"<<T<<">\n";
				std::vector<double> ns,nc;
				double temp;
				Polynomial expe=PolynomialGenerator(points).getPolynomial();
				cout << expe << "\n";
				for(int i=0;i<points.size();i++) cout << points[i].getX()<<","<<expe .getVal((double)(points[i].getX())) << "\n";
				double fc=expe.integrate(-M_PI,M_PI)/(2*M_PI);
				cin >> temp;
				for(int i=1;i<=T;i++)
				{
					cout << "Begin: "<< i << "\n";
					std::vector<Point2D> sps=points,cps=points;
					for(int j=0;j<points.size();j++)
					{
						sps[j].setY(sps[j].getY()*::sin(i*sps[j].getX()));
						cps[j].setY(cps[j].getY()*::cos(i*cps[j].getX()));
					}
					cout << "\tSins\n";
					for(Point2D& p : sps) cout <<"\t\t"<< p.getX()<<","<<p.getY() << "\n";
					cout << "\tCoSins\n";
					for(Point2D& p : sps) cout <<"\t\t"<< p.getX()<<","<<p.getY() << "\n";
					cout << "\tGenerated points: \n";
					Polynomial CSP=PolynomialGenerator(sps).getPolynomial();
					cout <<"\t"<< CSP << "\n";
					ns.push_back(CSP.integrate(-M_PI,M_PI)/M_PI);
					cout << "\tsins ready\n";
					CSP=PolynomialGenerator(cps).getPolynomial();
					cout << "\t" << CSP << "\n";
					nc.push_back(CSP.integrate(-M_PI,M_PI)/M_PI);
					cout << "\tcoses ready\n";
				}
				return FourierSeries<T>(ns,nc,fc);
			}
		};
	}
}

/*namespace geom2D
{
	class Point2D;
	static std::ostream& operator<<(std::ostream& os,const Point2D& pt);
	class Point2D
	{
	public:
		double x=0,y=0;
		Point2D() {}
		Point2D(double xc,double yc=0) {x=xc;y=yc;}
		Point2D(const maths::Vector& v) {x=v.getX(); y=v.getY();}


	};
}*/
namespace geometry
{
	template<int d> class Point
	{
		std::vector<double> p;
		double dumpdim=0;
	public:
		Point()
		{
			p=std::vector<double>();
			for(int i=0;i<d;i++) p.push_back(0);
		}
		Point(double* ps,int n=d)
		{
			for(int i=0;i<n;i++) p.push_back(ps[i]);
		}
		Point(const maths::Vector& v)
		{
			for(int i=0;i<min(v.getSize(),d);i++) p.push_back(v[i][0]);
		}

		double& operator[](const int& x)
		{
			if(x<d) {return p[x];}
			throw commons::IndexOutOfBoundsException();
		}
		const double& operator[](const int& x) const
		{
			if(x<d) {return p[x];}
			else {return dumpdim;}
		}
		double distanceFromOrigin() const
		{
			double s=0;
			for(int i=0;i<d;i++)
				s+=p[i]*p[i];
			return sqrt(s);
		}

		template<int d2> maths::Vector toVector(const Point<d2>& pt=Point<d>()) const
		{
			maths::Vector ret(max(d,d2));
			for(int i=0;i<max(d,d2);i++)
				ret[i][0]=p[i]-pt[i];
			return ret;
		}
		//Template functions
		template<int d2> double dist(const Point<d2>& p2) const { return distanceFromOrigin(p2-(*this));}
		template<int d2> Point<max(d,d2)> operator-(const Point<d2>& p2) const
		{
			Point<max(d,d2)> ret;
			for(int i=0;i<max(d,d2);i++)
				ret[i]=p2[i]-p[i];
			return ret;
		}
		Point<d> operator-() const
		{
			Point<d> ret;
			for(int i=0;i<d;i++)
				ret[i]=-p[i];
			return ret;
		}
		template<int d2> Point<max(d,d2)> operator+(const Point<d2>& p2) const
		{
			Point<max(d,d2)> ret;
			for(int i=0;i<max(d,d2);i++)
				ret[i]=p2[i]+p[i];
			return ret;
		}
	};
	template<int d> static std::ostream& operator<<(std::ostream& os,const Point<d>& pt);
}
namespace geom3D
{
	typedef geometry::Point<3> Point;
	//Class declaration
	class Point3D;
	//I/O Methods
	static std::ostream& operator<<(std::ostream& os,const Point3D& pt);

	//Basic Classes
	class Point3D //: public geometry::Point<3>
	{
	public:
		double x=0,y=0,z=0;
		/*Point3D() : geometry::Point() {}
		Point3D(double xc,double yc=0,double zc=0) : geometry::Point(new double[] {xc,yc,zc}) {}
		Point3D(const maths::Vector& v1) : geometry::Point(v1) {}
		Point3D(const geometry::Point<3> pt) : geometry::Point(pt[0],pt[1],pt[2]) {}*/
		Point3D() {}
		Point3D(double xc,double yc=0,double zc=0) {x=xc;y=yc;z=zc;}
		Point3D(const maths::Vector& v1)
		{
			int s=v1.getSize();
			if(s>0) x=v1[0][0];
			if(s>1) y=v1[1][0];
			if(s>2) z=v1[2][0];
		}
		Point3D(const geometry::Point<3> pt) : Point3D(pt[0],pt[1],pt[2]) {}
		//Point3D(const Point3D& pt) DEFAULT

		double dist(const Point3D& p2) const;
		maths::Vector toVector(const Point3D& p2) const;

		//Operator overloading
		inline Point3D operator+(const Point3D& pt) const {return Point3D(x+pt.x,y+pt.y,z+pt.z);}
		inline Point3D operator-() const {return Point3D(-x,-y,-z);}
		inline Point3D operator-(const Point3D& p2) const {return Point3D(x-p2.x,y-p2.y,z-p2.z);}

	};
	class Dimension3D
	{
		double dimX,dimY,dimZ;
	public:
		Dimension3D(double x=0,double y=0,double z=0) {dimX=x; dimY=y; dimZ=z;}

		inline bool isPoint() const {return (dimX==0 && dimY==0 && dimZ==0);}
		inline bool isLine()  const {return (((dimX==0 || dimY==0) && dimZ==0) || (dimX==0 && dimY==0));}
		inline bool isPlane() const {return (dimX==0 || dimY==0 || dimZ==0);}
		inline bool isSpace() const {return (dimX!=0 && dimY!=0 && dimZ!=0);}
		inline double getX() const {return dimX;}
		inline double getY() const {return dimY;}
		inline double getZ() const {return dimZ;}
	};
	class EulerAngle
	{
		double xy,yz,zx;
	public:
		EulerAngle() {xy=yz=zx=0;}
		EulerAngle(double a1,double a2=0,double a3=0,bool conv=true)
		{
			if(maths::isDegrees() && conv)
			{
				a1=maths::toRadians(a1);
				a2=maths::toRadians(a2);
				a3=maths::toRadians(a3);
			}
			xy=a1; yz=a2; zx=a3;
		}
		/*EulerAngle(Vector v) //Gives the Euler angle required for turning x axis to vector v
		{
			v=v/v.getMagnitude();
			xy=cosinv(v.getX());
		}*/

		maths::Matrix* getRotationMatrix() const;
		//Operator overloading
		EulerAngle operator-() const {return EulerAngle(-xy,-yz,-zx,false);}
	};
}

//NAMESPACE maths
using namespace maths;
//Other functions from above classes
Vector Matrix::getColumn(int n) const
{
	Vector v(getRows());
	for(int i=0;i<hei;i++)
		v[i][0]=mrows[i][n];
	return v;
}
static Matrix operator*(const double& d,const Matrix& m) {return m*d;}
inline Polynomial Matrix::getCharDeterminant()  {return (new PolyMatrix(*this,EQUALITY,Polynomial::CharPolynomial))->determinant();}
std::vector<double> Matrix::getEigenvalues()
{
	if(getRows()!=getColumns())
		throw IncompatibleMatrixOperationException("Requested eigenvalue calculation for rectangular matrix.");
	return getCharDeterminant().getRoots();
}
std::vector<Vector> Matrix::getEigenvectors()
{
	std::vector<double> EVs=getEigenvalues();
	std::vector<Vector> ret;
	Matrix m;
	Vector r(getRows()),J;
	Vector res;
	EquationSystem ES,oES;
	bool solv=false;
	int K=0;
	for(double ev : EVs)
	{
		J=Vector(getRows());
		m=(*this)-(getIdentity(getColumns())*ev);
		ES=EquationSystem(m,r);
		if(ES.solvable())
		{
			cout << m << "\n";
			cout << "Eigenvalue: "<<ev<<": Non-zero determinant: "<<m.determinant()<<" to zero value.\n";
			continue;
		}
		for(int k=0;k<getRows();k++)
		{
			oES=ES.reduce(k,1,k);
			if(!oES.solvable())
				continue;
			res=oES.solve();
			for(int i=0;i<k;i++)
				J[i][0]=res[i][0];
			J[k][0]=1;
			for(int i=k+1;i<getRows();i++)
				J[i][0]=res[i-1][0];
			ret.push_back(J);
			solv=true;
			break;
		}

		/*while(ES.getVarCount()!=1)
		{
			int k=0;
			while(k<getRows())
			{
				if(mrows[k].isZero())
					break;
				k++;
			}
			if(k>=getRows())
				k=0;

			J[k][0]=1;
		}*/
		if(solv)
			continue;
		int i;
		for(i=0;i<getColumns();i++)
		{
			if(m.getColumn(i).isZero())
				break;
		}
		if(i==getColumns())
			i=0;
		res=Vector(getRows());
		res[i][0]=1;
		J=res;
		//cout << J << "\n";
		ret.push_back(J);
		/*
		 * (A-l)v=0
		 * 00 10
		 * 01 11
		 */
	}
	return ret;
}
Matrix* Matrix::diagonalize()
{
	std::vector<Vector> EVs=getEigenvectors();
	std::vector<double> evs=getEigenvalues();
	cout << evs.size() << "\n";
	Matrix P=Vector::buildMatrix(EVs);
	//cout << P << "\n";
	Matrix PI=P.getInverse();
	Matrix D=(getRows(),getColumns());
	int K=0;
	for(int i=0;i<D.getRows();i++)
	{
		for(int j=0;j<D.getColumns();j++)
		{
			if(i!=j)
				continue;
			D[i][j]=evs[K++];
		}
	}
	return new Matrix[3] {P,D,PI};
}
//Input methods
static std::istream& operator>>(std::istream& is,Matrix& mtr)
{
	std::string in="";
	commons::String str;
	int nC=0;
	Matrix ret;
	while(in!="done")
	{
		if(!std::getline(is,in))
			break;
		if(in=="done")
			break;
		ret.addRow(MatrixRow(in));
	}
	mtr=ret;
	return is;
}
static std::ostream& maths::operator<<(std::ostream& os,Matrix mtr)
{
	if(mtr.getRows()==0)
	{
		os << "<Empty Matrix>";
		return os;
	}
	double acc=pow(10,-PRECISION);
	for(int i=0;i<mtr.getRows();i++)
	{
		for(int j=0;j<mtr[i].getSize();j++)
			os << ((abs(mtr[i][j])>=acc)?mtr[i][j]:0)<<"\t";
		os <<"\n";
	}
	return os;
}
static std::ostream& maths::operator<<(std::ostream& os,PolyMatrix mtr)
{
	for(int i=0;i<mtr.getRows();i++)
	{
		for(int j=0;j<mtr[i].getSize();j++)
			cout << mtr[i][j]<<"\t";
		cout <<"\n";
	}
	return os;
}
static std::ostream& operator<<(std::ostream& os,Interval in)
{
	os << "("<<in.getLowerBound()<<","<<in.getUpperBound()<<")";
	return os;
}
static std::ostream& operator<<(std::ostream& os,Histogram hs)
{
	os << "Histogram\n";
	Interval* it=hs.getIntervals();
	int lC=hs.getIntervalCount();
	it=sort(it,lC);
	for(int i=0;i<lC;i++)
	{
		os << it[i]<<":\t";
		os << hs.getCount(it[i]) <<"\n";
	}
	return os;
}
inline static std::ostream& operator<<(std::ostream& os,const Point2D& pt) {return os << "("<<pt.getX()<<","<<pt.getY()<<")";}
template<int T> static std::ostream& operator<<(std::ostream& os,const maths::fourier::FourierSeries<T>& fs)
{
	os << fs.getConstant();
	for(int i=0;i<T;i++)
	{
		if(fs.getSin(i)) os <<"+"<< fs.getSin(i) << ".sin("<<(i?to_string(i+1):"")<<"x) ";
		if(fs.getCos(i)) os <<"+"<< fs.getSin(i) << ".cos("<<(i?to_string(i+1):"")<<"x) ";
	}
	return os;
}
static std::ostream& maths::operator<<(std::ostream& os,Polynomial p)
{
	char ch=p.getVar();
	double* cfs=p.getCoeffs();
	const int deg=p.getDegree();
	if(deg==0)
	{
		os << p.getVal(0);
		return os;
	}
	for(int i=deg;i>=0;i--)
	{
		if(cfs[i]==0)
			continue;
		switch(i)
		{
			case 0:
				os << (cfs[0]>0?" +":" ")<<cfs[0];
				break;
			case 1:
				os << (cfs[1]>0?" +":" -");
				if(cfs[1]!=1 && cfs[1]!=-1)
					os << maths::abs(cfs[1]);
				os << ch;
				break;
			default:
				if(i!=deg)
					os << (cfs[i]>0?" +":" -");
				else
					os << ((cfs[i]>0)?"":"-");
				if(cfs[i]!=1 && cfs[i]!=-1)
					os << maths::abs(cfs[i]);
				os << ch<<"^"<<i;
				break;
		}
	}
	return os;
}
static std::ostream& geom3D::operator<<(std::ostream& os, const geom3D::Point3D& p)
{
	double acc=pow(10,-maths::PRECISION);;
	os << "("<<((maths::abs(p.x)>=acc)?p.x:0)<<","<<((maths::abs(p.y)>=acc)?p.y:0)<<","<<((maths::abs(p.z)>=acc)?p.z:0)<<")";
	return os;
}
//Functions from geom3D
double geom3D::Point3D::dist(const Point3D& p2=Point3D()) const {return sqrt((p2.x-x)*(p2.x-x)+(p2.y-y)*(p2.y-y)+(p2.z-z)*(p2.z-z));}
maths::Vector geom3D::Point3D::toVector(const geom3D::Point3D& p2=geom3D::Point3D())const {return Vector(x-p2.x,y-p2.y,z-p2.z);}

//Functions for vector rotation
static Matrix rotateZAxisTo(const Vector& V)
{
	Vector v=V/V.getMagnitude();
	double tanphi=v.getY()/v.getX();
	if(tanphi!=tanphi) //tanphi is NaN
	{return Matrix::getRotation(0,0);}
	double phi=taninv(tanphi);
	double tphi=phi;
	if(maths::isDegrees())
		tphi=toRadians(phi);
	//cout << sin(phi) << " "  << cos(phi) << " " <<phi <<"\n";
	double tantheta=(v.getX()*cos(tphi)+v.getY()*sin(tphi))/v.getZ();
	if(tantheta!=tantheta) // tantheta is NaN
		tantheta=(v.getX()<=0)?-1:1;
	double theta=taninv(tantheta);
	return Matrix::getRotation(theta,phi);
}
inline static Matrix rotateXAxisTo(const Vector& v) {return rotateZAxisTo(Vector(1,0,0)).getInverse()*rotateZAxisTo(v);}
inline static Matrix rotateYAxisTo(const Vector& v) {return rotateZAxisTo(Vector(0,1,0)).getInverse()*rotateZAxisTo(v);}

inline static Matrix getTransformation(const geom3D::Point3D& p1,const geom3D::Point3D& p2=geom3D::Point3D()){return rotateZAxisTo(p2.toVector(p1));}
static Vector transformVector(const Vector& v,const geom3D::Point3D& p1,geom3D::Point3D p2=geom3D::Point3D()) //Returns the vector as it would be in frame when p1 is origin and p1-p2 is Z-axis
{
	//geom3D::Point3D ORIGIN=p1;
	Vector V=v-p1.toVector();
	//cout  << p1 << "\t" << p2 << "\n";
	Matrix M=getTransformation(p1,p2);
	//cout << M << "\n";
	V=V.transform(M);
	return V;
}
inline static Matrix getReverseTransformation(const geom3D::Point3D& p1,const geom3D::Point3D& p2=geom3D::Point3D()) {return getTransformation(p1,p2).getInverse();}
static Vector reverseTransformVector(const Vector& v,const geom3D::Point3D& p1,const geom3D::Point3D& p2=geom3D::Point3D()) //Returns the vector in normal co-ordinates if it was v when p1-p2 is Z-axis
{
	//Vector V=v-p1.toVector();
	Vector V;
	Matrix rot=getReverseTransformation(p1,p2);
	try {V=v.transform(rot);}
	catch(IncompatibleMatrixOperationException ex) {cout << "No such possible config.\n";}
	return (V+p1.toVector());
}
static Matrix rotate2D(double angle) //Gives 2D Rotation about z-axis by angle "theta"
{
	if(maths::isDegrees())
		angle=toRadians(angle);
	Matrix ret(3,3);
	ret[2][2]=1;
	ret[0][0]=cos(angle);
	ret[1][0]=sin(angle);
	ret[0][1]=-sin(angle);
	ret[1][1]=cos(angle);
	return ret;
}
static Matrix getRotationTransformAbout(const geom3D::Point3D& p1,const geom3D::Point3D& p2,double theta) {return getReverseTransformation(p1,p2)*rotate2D(theta)*getTransformation(p1,p2);}

static Vector rotateVectorAbout(const Vector& V,const geom3D::Point3D& p1,const geom3D::Point3D& p2,double theta)
{
	Vector v=V-p1.toVector();
	/*Vector v1=transformVector(v,p1,p2);
	//cout << v1 << "\n";
	v1=v1.transform(rotate2D(theta));
	//cout << v1 << "\n";
	return reverseTransformVector((v1),p1,p2);*/
	Matrix rot=getRotationTransformAbout(p1,p2,theta);
	Vector w=v.transform(rot);
	return (w+p1.toVector());
}
static Vector rotateVectorAbout(const Vector& V,const Vector& dir,double theta) {return rotateVectorAbout(V,geom3D::Point3D(0,0,0),geom3D::Point3D(dir),theta);}
static Matrix rotateFrame(const Vector& X,const Vector& Y)
{
	Vector Z=X.cross(Y);
	//cout << Z << "\n"; //Correct
	Matrix R=rotateZAxisTo(Y);
	Vector t=R*X;
	t=t/t.getMagnitude();
	double a=cosinv(t.dot(Vector(1,0,0)));
	return rotate2D(a).getInverse()*R;
}

//Missing methods in geom3D
Matrix* geom3D::EulerAngle::getRotationMatrix() const
{
	Matrix* ret=new Matrix(3,3);
	double XY=xy,YZ=yz,ZX=zx;
	if(maths::isDegrees())
	{
		XY=toDegrees(XY);
		YZ=toDegrees(YZ);
		ZX=toDegrees(ZX);
	}
	//cout << "XY: " << XY<<"\t"<<"YZ: "<<YZ<<"\t"<<"ZX: "<<ZX<<"\n";
	geom3D::Point3D ORIGIN=*(new Point3D(0,0,0));
	*ret=getRotationTransformAbout(ORIGIN,geom3D::Point3D(0,0,1),XY);
	*ret=getRotationTransformAbout(ORIGIN,geom3D::Point3D(1,0,0),YZ)*(*ret);
	*ret=getRotationTransformAbout(ORIGIN,geom3D::Point3D(0,1,0),ZX)*(*ret);
	return ret;
}

//Other rotations
static Vector rotateVector(const Vector& v,const geom3D::EulerAngle& EA,const geom3D::Point3D& about=geom3D::Point3D())
{
	Vector V=v-about.toVector();
	Matrix* rot=EA.getRotationMatrix();
	Vector w=V.transform(*rot);
	delete rot;
	return (w+about.toVector());
}

//Angle relations
double angleDegrees(double deg)
{
	if(maths::radians)
		return toRadians(deg);
	else
		return deg;
}

//Other operations
static double sum(const Matrix& m)
{
	double s=0;
	for(int i=0;i<m.getRows();i++)
	{
		for(int j=0;j<m.getColumns();j++)
			s+=m[i][j];
	}
	return s;
}
#endif
