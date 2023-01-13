#ifndef INCLUDED_IMAGES
#define INCLUDED_IMAGES 1
#include "general/general.h"
//#include "draw/draw.h"

namespace draw
{
  class BigImageException : public std::exception {};
  namespace images
  {
    inline static double matrixAverage(const NumericMatrix<double>& m) {return m.average();}
    class AssignmentInSubdividedImageException : public exception {};
    class AlreadySubdividedImageException : public exception {};
    class ImageDimensionsInsufficientException : public exception {};
    class BlackAndWhiteImage
    {
      NumericMatrix<double> imgdata;
      bool sub=false;
      int subc=1,subr=1;
      double drawco=0.35;
    public:
      BlackAndWhiteImage() {}
      BlackAndWhiteImage(int w,int h) {imgdata=NumericMatrix<double>(w,h);}
      BlackAndWhiteImage(const GeneralVector<double>& data,int w,int h) : BlackAndWhiteImage(w,h)
      {
        int K=0;
        for(int i=0;i<h;i++)
        {
          for(int j=0;j<w;j++) imgdata[i][j]=data(K++);
        }
      }
      explicit BlackAndWhiteImage(const NumericMatrix<double>& dat) {imgdata=dat;}

      inline const NumericMatrix<double>& getImageData() const {return imgdata;}
      inline double at(int x, int y) const
      {
        if(!sub) {return imgdata[y][x];}
        else {return imgdata[y/subr][x/subc];}
      }
      inline void set(int x,int y,double v)
      {
        if(sub) throw AssignmentInSubdividedImageException();
        else imgdata[y][x]=v;
      }
      inline int getWidth() const
      {
        if(!sub) return imgdata.getColumnCount();
        else return imgdata.getColumnCount()*subc;
      }
      inline int getHeight() const
      {
        if(!sub) return imgdata.getRowCount();
        else return imgdata.getRowCount()*subr;
      }
      GeneralVector<double> unwindData() const
      {
        GeneralVector<double> ret(this->getWidth()*this->getHeight());
        int K=0;
        for(int i=0;i<this->getHeight();i++)
        {
          for(int j=0;j<this->getWidth();j++) ret(K++)=at(j,i);
        }
        return ret;
      }


      inline void useCutoff(double v) {drawco=v;}
      bool drawValAt(int x,int y) const {return at(x,y)>=drawco;}
      inline void subdivideRows(int n) {subdivide(n,1);}
      inline void subdivideColumns(int n) {subdivide(1,n);}
      inline void subdivide(int r,int c) {subc=c; subr=r; sub=true;}
      inline void snap() {sub=false; subc=subr=1;}
      BlackAndWhiteImage getSubImage(int x,int y,int w,int h) const
      {
        BlackAndWhiteImage ret(w,h);
        for(int i=x;i<x+w;i++)
        {
          for(int j=y;j<y+h;j++) {ret.set(i-x,j-y,this->at(i,j));}
        }
        return ret;
      }

      BlackAndWhiteImage rescale(int x,int y,double (*pf)(const NumericMatrix<double>&)=&matrixAverage)
      {
        BlackAndWhiteImage ret(x,y);
        if(sub) throw AlreadySubdividedImageException();
        int ow=this->getWidth(),oh=this->getHeight();
        this->subdivide(x,y);
        for(int i=0;i<x;i++)
        {
          //for(int j=0;j<y;j++) ret.set(i,j,this->getSubImage(i*ow,j*oh,ow,oh).getImageData().average());
          for(int j=0;j<y;j++) ret.set(i,j,pf(this->getSubImage(i*ow,j*oh,ow,oh).getImageData()));
        }
        this->snap();
        return ret;
      }
      BlackAndWhiteImage& pasteImage(const BlackAndWhiteImage& img,int x,int y)
      {
        if(x+img.getWidth()>=this->getWidth() || y+img.getHeight()>=this->getHeight()) throw ImageDimensionsInsufficientException();
        for(int i=x;i<x+img.getWidth();i++)
        {
          for(int j=y;j<y+img.getHeight();j++) this->set(i,j,img.at(i-x,j-y));
        }
        return *this;
      }
      BlackAndWhiteImage erasePadding(double cut=0.1,bool rev=false)
      {
        int huc=this->getHeight(),hlc=0;
        int vlc=this->getWidth(),vuc=0;
        bool vin=false,hin=false,hf=false;
        for(int j=0;j<this->getHeight();j++)
        {
          hf=false;
          for(int i=0;i<this->getWidth();i++)
          {
            double v=this->at(i,j);
            if((!rev && v>cut) || (rev && v<cut))
            {
              hf=true;
              //cout <<"\t"<< i <<","<< j << "\n";
              if(vlc>i) {vlc=i;}
              else if(vuc<i) {vuc=i;}
            }
          }
          if(!hf)  //Empty row
          {if(!hin) hlc=j;}
          else
          {
            //cout << "Filled row: "<<j << "\n";
            huc=j;
            hin=true;
          }
        }
        //cout << "Image taken from ("<<vlc<<"x"<<(hlc+1)<<") to ("<<vuc<<"x"<<huc<<")\n";
        BlackAndWhiteImage ret;
        if(vuc<vlc) ret=BlackAndWhiteImage(0,0);
        else ret=getSubImage(vlc,hlc+1,vuc-vlc,huc-hlc-1);
        ret.useCutoff(drawco);
        return ret;
      }
      BlackAndWhiteImage addPadding(int fbuf,bool bidir=true,int bbuf=0)
      {
        if(bbuf) bidir=false;
        BlackAndWhiteImage ret; ret.useCutoff(drawco);
        if(bidir) ret=BlackAndWhiteImage(this->getWidth()+2*fbuf,this->getHeight()+2*fbuf);
        else ret=BlackAndWhiteImage(this->getWidth()+fbuf+bbuf,this->getHeight()+fbuf+bbuf);
        ret.pasteImage(*this,fbuf,fbuf);
        return ret;
      }
    };
    class ImageProcessor
    {
    public:
      ImageProcessor() {}

      virtual BlackAndWhiteImage processImage(const BlackAndWhiteImage& img) const=0;
    };
  }
}

static std::ostream& operator<<(std::ostream& os,draw::images::BlackAndWhiteImage img)
{
  for(int i=0;i<img.getHeight();i++)
  {
    for(int j=0;j<img.getWidth();j++) os << (int)(img.drawValAt(j,i)) << " ";
    os << "\n";
  }
  return os;
}
#endif
