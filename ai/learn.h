#ifndef INCLUDED_AILEARN
#define INCLUDED_AILEARN 1
#include "maths/stats.h"
#include "draw/images.h"
#include "general/stats.h"
#ifndef EPSILON
#define EPSILON 0.02
#endif

namespace learning
{
  static double unittol=0.1;
  typedef std::vector<std::vector<double>> DataSet;
  inline static double euclideanDistanceSquare(const GeneralVector<double>& v1,const GeneralVector<double>& v2) {return GeneralVector<double>(v2-v1).norm2();}
  inline static double classificationErrorMeasure(const GeneralVector<double>& res,const GeneralVector<double>& out)
  {
    double m=-1; int ind=-1,expind=-1;
    for(int i=0;i<res.size();i++)
    {
      if(m<res(i)) {m=res(i); ind=i;}
    }
    for(int j=0;j<out.size();j++)
    {
      if(out(j)) {expind=j; break;}
    }
    return !(expind==ind);
  }
  inline static double binaryErrorMeasure(const GeneralVector<double>& res,const GeneralVector<double>& out)
  {
    if(res.getSize()!=out.getSize()) return max(res.getSize(),out.getSize());
    for(int i=0;i<res.size();i++)
    {
      if(::abs(res(i)-out(i))>unittol) return 1;
    }
    return 0;
  }
  class ParameterHolder
  {
  public:
    inline virtual void refresh()=0;
    inline virtual void scaleLearnRates(double s=1)=0;
  };
  class ValuePredictor : public ParameterHolder
  {
  public:
    virtual double evaluate(const GeneralVector<double>&) const=0;
  };
  class Predictor : public flex::Cloneable,public ParameterHolder
  {
    public:
      virtual GeneralVector<double> evaluate(const GeneralVector<double>&)const =0;
      virtual void train(const GeneralVector<double>& in,const GeneralVector<double>& exp)=0;
      virtual Predictor* clone() const override=0;
  };
  class PreservablePredictor : public Predictor
  {
  public:
    PreservablePredictor() : Predictor() {}

    virtual void save(std::ostream&) const =0;
    virtual void load(std::istream&,std::ostream&)=0;
  };
  class NoModelsForBaggingException : public exception {};
  class NoModelsForPredictorException : public exception {};
  class UnPackedNeuralNetworkException : public exception {};
  class NotAllModelsCanBeSavedException : public exception {};
  class UndeterminedLoadTypesException : public exception {};
  class InsufficientDataToTrainException : public exception {};
  //void fitstep(std::vector<double>& weights,double (T::*f)(const std::vector<double>&) const,double exp,const T& obj,bool twosided=true,double D=0.01)
  template<class T> class NumericModel;
  template<class T> void fitstep(NumericModel<T>& model,double exp,bool twosided=true,double D=1e-3,double L=0.05,int acc=5);

  template<class T> class NumericModel
  {
  protected:
    T data;
    std::vector<double> params;
  public:
    double k=0;
    NumericModel() {}

    virtual T& loadData(const T& d) {return data=d;}
    inline const T& getData() const {return data;}
    double predict(const T& d) {loadData(d); return predict();}

    virtual double predict() const=0;

    void learnstep(const T& d,double ac)
    {
      loadData(d);
      fitstep(*this,ac);
    }
    inline void learn(const std::vector<T>& datas,const std::vector<double>& vals) {for(int i=0;i<min(datas.size(),vals.size());i++) learnstep(datas[i],vals[i]);}
    inline void train(const std::vector<T>& datas,const std::vector<double>& vals) {learn(datas,vals);}
    inline std::vector<double>& getParams() {return params;}
    inline const std::vector<double>& getParams() const {return params;}
    virtual std::vector<double> gradient(double exp,double D=1e-6,int acc=3)
    {
      std::vector<double> grad;
      double opred=predict();
      //double opred=exp;
      for(int i=0;i<params.size();i++)
      {
        params[i]+=D;
        //cout << predict() << "\t" << opred << "\n";
        grad.push_back(maths::round((predict()-opred),acc)/D);
        params[i]-=D;
      }
      return grad;
    }
    std::vector<double> gradientAt(const T& d,double D=1e-6,int acc=3)
    {
      loadData(d);
      return gradient(predict(),D,acc);
    }
    void reduce(int acc=3)
    {
      for(int i=0;i<params.size();i++) params[i]=maths::round(params[i],acc);
    }
  };
  template<class T> void fitstep(NumericModel<T>& model,double exp,bool twosided,double D,double L,int acc)
  {
    //cout << "Fitstep\n";
    std::vector<double> nw;
    std::vector<double> dev,grad;
    double opred=model.predict();
    //double orerr=(opred-exp)*L;
    cout <<"Prediction: " << opred << "\t for expected: " << exp << "\n";// " using L*orerr = "<<orerr<<"\n";
    if(::abs(opred-exp)<pow(10,-acc)) return;
    std::vector<double> gradvec=model.gradient(exp,D,acc);
    cout << "Gradient Vector: "<< gradvec << "\n";
    for(int i=0;i<model.getParams().size();i++) {grad.push_back(L*(opred-exp)*min(max(gradvec[i],-1),1));}
    /*model.k+=D;
    double g=L*(model.predict()-exp)/D;
    model.k-=D; model.k-=g;*/
    cout << "Gradient: "<<grad; // << "\n";
    GeneralVector<double> gv(grad); //if(!gv.norm()) return;
    //gv.normalize(); gv=gv*L;
    cout << gv << "\n";
    //cout << "Deviation: "<<dev; // << "\n";
    for(int i=0;i<model.getParams().size();i++) model.getParams()[i]-=gv(i);
    cout << "Fixed: "<<model.getParams();
  }


  template<int T=1> class VectorDotModel : public NumericModel<GeneralVector<double>>
  {
  public:
    VectorDotModel(int psize=T) : NumericModel<GeneralVector<double>>()
    {
      this->params=std::vector<double>(psize); for(int i=0;i<psize;i++) params[i]=1;
    }
    inline virtual double predict(const GeneralVector<double>& d) {return NumericModel<GeneralVector<double>>::predict(d);}
    inline virtual double predict() const override {return GeneralVector<double>(params).dotIgnore(data)+this->k;}

    void fitData(const std::vector<GeneralVector<double>>& datas,const std::vector<double>& vals)
    {
      int tempv;
      Eigen::MatrixXd mymat(params.size(),params.size());
      for(int i=0;i<params.size();i++)
      {
        for(int j=i;j<params.size();j++)
        {
          mymat(i,j)=0;
          for(int k=0;k<datas.size();k++)
          {
            const GeneralVector<double>& x=datas[k];
            mymat(i,j)+=x(i)*x(j);
          }
        }
      }
      for(int i=0;i<params.size();i++)
      {
        for(int j=0;j<i;j++) mymat(i,j)=mymat(j,i);
      }
      cout << mymat << "\n";
      cout << mymat.determinant() << "\n";
      cin >> tempv;
      Eigen::VectorXd ansv(params.size());
      for(int i=0;i<params.size();i++)
      {
        ansv(i)=0;
        for(int j=0;j<datas.size();j++) ansv(i)+=vals[j]*datas[j](i);
      }
      cout << ansv << "\n\n\n";
      cout << mymat.inverse() << "\n";
      Eigen::VectorXd ret=mymat.inverse()*ansv;
      for(int i=0;i<params.size();i++) params[i]=ret(i);
    }
    std::vector<double> gradient(double exp,double D=1e-13,int acc=4) override
    {
      return data;
      /*double p=predict();
      std::vector<double> npar;
      for(int i=0;i<params.size();i++) npar.push_back(data[i]*/
      //return GeneralVector(npar).dotIgnore(data);
    }
  };
  template<int T> class MatrixDotModel : public NumericModel<NumericMatrix<double>>
  {
  public:
    MatrixDotModel() : NumericModel<NumericMatrix<double>>() {this->params=std::vector<double>(2*T); for(int i=0;i<2*T;i++) params[i]=1;}

    inline virtual double predict(const NumericMatrix<double>& d) {return NumericModel<NumericMatrix<double>>::predict(d);}
    inline virtual double predict() const override
    {
      auto p=generateVectors(); GeneralVector<double> v1=get<0>(p),v2=get<1>(p);
      //cout << v1.size() <<"\t"<<v2.size()<<"\n";
      return (v1.transpose()*data*v2)[0][0];
      //return GeneralVector<double>(params).dotIgnore(data)+this->k;
    }
    inline virtual std::pair<GeneralVector<double>,GeneralVector<double>> generateVectors() const
    {
      GeneralVector<double> v1(params.size()/2),v2(params.size()/2);
      for(int i=0;i<params.size()/2;i++) v1(i)=params[i];
      for(int i=params.size()/2;i<params.size();i++) v2(i-params.size()/2)=params[i];
      return make_pair(v1,v2);
    }
  };
  namespace imgproc
  {
    class SimpleSlider : public draw::images::ImageProcessor
    {
      int dimX=1,dimY=1;
      int stepX=1,stepY=1;
    public:
      SimpleSlider(int dx,int dy) : SimpleSlider(dx,dy,1,1) {}
      SimpleSlider(int dx,int dy,int sx,int sy)
      {
        dimX=dx; dimY=dy;
        stepX=sx; stepY=sy;
      }

      draw::images::BlackAndWhiteImage processImage(const draw::images::BlackAndWhiteImage& img) const override
      {
        double iw=(img.getWidth()-img.getWidth()%dimX)/stepX,ih=(img.getHeight()-img.getHeight()%dimY)/stepY;
        draw::images::BlackAndWhiteImage ret(iw,ih);
        int i=0,j=0;
        int R=0,C=0;
        //cout << "Choosing final size: "<<ret.getWidth()<<","<<ret.getHeight()<<"\n";
        while(i+dimY<=img.getHeight())
        {
          j=0;
          C=0;
          while(j+dimX<=img.getWidth())
          {
            ret.set(C++,R,evaluate(img.getSubImage(j,i,dimX,dimY)));
            j+=stepX;
            //cout << j << "\t" << img.getWidth() << "\n";
          }
          R++;
          i+=stepY;
          //cout <<"\t"<< i << "\t" << img.getHeight() << "\n\n";
        }
        return ret;
      }
      inline virtual double evaluate(const draw::images::BlackAndWhiteImage& img) const {return img.getImageData().average();}
    };
    namespace pooling
    {
      class MaxPooling : public SimpleSlider
      {
      public:
        MaxPooling(int dx,int dy) : MaxPooling(dx,dy,1,1) {}
        MaxPooling(int dx,int dy,int sx,int sy) : SimpleSlider(dx,dy,sx,sy) {}

        inline virtual double evaluate(const draw::images::BlackAndWhiteImage& img) const {return general::stats::matrixMax(img.getImageData());}
      };
      class MinPooling : public SimpleSlider
      {
      public:
        MinPooling(int dx,int dy) : MinPooling(dx,dy,1,1) {}
        MinPooling(int dx,int dy,int sx,int sy) : SimpleSlider(dx,dy,sx,sy) {}

        inline virtual double evaluate(const draw::images::BlackAndWhiteImage& img) const {return general::stats::matrixMin(img.getImageData());}
      };
    }
  }
  //template<int t> class ReducedParameterMatrixDotModel : public
  namespace neuralnetworks
  {
    class Neuron : public flex::Cloneable,public ValuePredictor
    {
    protected:
      int sz;
      GeneralVector<double> weights;
      double cutoff,bias=0,biasweight=0;
      double learnrate=0.01;
    public:
      int type=0;
    public:
      Neuron(int s,double v=0,double lr=0.1) : Neuron(s,GeneralVector<double>(s),v,lr) {}
      Neuron(int s,GeneralVector<double> extw,double v=0,double lr=0.1) {sz=s; cutoff=v; weights=extw; learnrate=lr;}

      inline int getInputCount() const {return sz;}
      inline void setLearnRate(double lr) {learnrate=lr;}
      inline double getLearnRate() const {return learnrate;}
      inline void addBias(double b=1) {bias=b;}
      inline double getBias() const {return bias;}
      inline void setBiasWeight(double bw) {biasweight=bw;}
      inline double getBiasWeight() const {return biasweight;}
      void negate()
      {
        for(int i=0;i<weights.getSize();i++) weights(i)*=-1;
        biasweight*=-1;
      }
      inline GeneralVector<double> getWeights() const
      {
        if(bias)
        {
          std::vector<double> ret=weights; ret.push_back(biasweight);
          return ret;
        }
        else return weights;
      }
      inline void randomizeWeights(double ll=0,double ul=1) {for(int i=0;i<weights.getSize();i++) weights(i)=randgen::throwarandompoint(ll,ul); biasweight=randgen::throwarandompoint(ll,ul);}

      inline virtual double evaluate(const GeneralVector<double>& in) const override {return weights.dot(in)+biasweight*bias;}
      inline virtual double operator()(const GeneralVector<double>& in) const {return evaluate(in)>=cutoff;}
      inline virtual double activatorDerivative(double x) const {return 1;} //Here x is the OUTPUT of the activaor (not the input)
      inline virtual Neuron* clone() const override
      {
        Neuron* ret= new Neuron(sz,weights,cutoff,learnrate);
        ret->addBias(bias);
        ret->setBiasWeight(this->getBiasWeight());
        return ret;
      }
      inline void refresh() override {this->randomizeWeights(-1,1);}
      inline void scaleLearnRates(double v=1) override {learnrate*=v;}
      inline virtual Neuron* reinit(int s,double v=0,double lr=0) const
      {
        if(!lr) lr=learnrate;
        return new Neuron(s,v,lr);
      }
      inline virtual Neuron* reinit(int s,GeneralVector<double> extw,double v=0,double lr=0) const
      {
        if(!lr) lr=learnrate;
        return new Neuron(s,extw,v,lr);
      }

      double adjust(const GeneralVector<double>& data,double exp)
      {
        double myval=this->operator()(data),der=activatorDerivative(myval)+EPSILON;
        double error=((exp-myval)*der)/weights.size();
        //cout << "Evaluated: "<<myval<<" with derivative: "<<der << "\n";
        if(error)
        {
          for(int i=0;i<weights.getSize();i++) weights(i)+=error*data(i)*(learnrate);
          if(bias) biasweight+=error*bias*learnrate;
        }
        return error;
      }
      double onlyAdjust(const GeneralVector<double>& data,double e)
      {
        //cout << "Adjust: "<<(data*e).transpose()<<"\n";
        weights=weights+(data*(e*learnrate));
        //for(int i=0;i<weights.getSize();i++) weights(i)+=e*data(i)*learnrate;
        if(bias) biasweight+=e*bias*learnrate;
        return e;
      }
    };
    class PlainNeuron : public Neuron
    {
      double scale=1;
    public:
      PlainNeuron(int s,double v=0,double lr=0.1) : Neuron(s,v,lr) {type=3;}
      PlainNeuron(int s,std::vector<double> w,double v=0,double lr=0.1) : Neuron(s,w,v,lr) {type=3;}

      inline double operator()(const GeneralVector<double>& in) const override {return scale*evaluate(in);}
      inline void setScale(double v) {scale=v;}
      inline double getScale() const {return scale;}
      inline virtual double activatorDerivative(double x) const override {return scale;}
      inline virtual PlainNeuron* clone() const override //Not a perfect clone (Returns a new neuron with same weights and no bias)
      {
        PlainNeuron* ret= new PlainNeuron(sz,weights,cutoff,learnrate);
        ret->addBias(bias);
        ret->setBiasWeight(this->getBiasWeight());
        return ret;
      }

      inline virtual PlainNeuron* reinit(int s,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        return new PlainNeuron(s,v,lr);
      }
      inline virtual PlainNeuron* reinit(int s,GeneralVector<double> extw,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        return new PlainNeuron(s,extw,v,lr);
      }
    };
    class SigmoidNeuron : public Neuron
    {
      double scale=1;
    public:
      SigmoidNeuron(int s,double v=0,double lr=0.1) : Neuron(s,v,lr) {type=1;}
      SigmoidNeuron(int s,std::vector<double> w,double v=0,double lr=0.1) : Neuron(s,w,v,lr) {type=1;}

      inline double operator()(const GeneralVector<double>& in) const override
      {
        double v=evaluate(in);
        //cout << v << " was evaluated\n";
        return scale/(1+exp(-v));
      }
      inline void setScale(double v) {scale=v;}
      inline double getScale() const {return scale;}
      inline virtual double activatorDerivative(double x) const override {return x*(1-x/scale);}
      inline virtual SigmoidNeuron* clone() const override //Not a perfect clone (Returns a new neuron with same weights and no bias)
      {
        SigmoidNeuron* ret= new SigmoidNeuron(sz,weights,cutoff,learnrate);
        ret->addBias(bias);
        ret->setBiasWeight(this->getBiasWeight());
        return ret;
      }

      inline virtual SigmoidNeuron* reinit(int s,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        return new SigmoidNeuron(s,v,lr);
      }
      inline virtual SigmoidNeuron* reinit(int s,GeneralVector<double> extw,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        return new SigmoidNeuron(s,extw,v,lr);
      }
    };
    class LinearNeuron : public Neuron
    {
    protected:
      double top=1,bot=0;
    public:
      LinearNeuron(int s,double v=0,double lr=0.1,double t=1,double b=0) : Neuron(s,v,lr) {type=2; top=t; bot=b;}
      LinearNeuron(int s,std::vector<double> w,double v=0,double lr=0.1,double t=1,double b=0) : Neuron(s,w,v,lr) {type=2; top=t; bot=b;}

      inline double operator()(const GeneralVector<double>& in) const override
      {
        double v=evaluate(in);
        if(v<bot) return bot;
        if(v>top) return top;
        return v;
      }
      inline void setMaxValue(double d) {top=d;}
      inline void setMinValue(double d) {bot=d;}
      inline double getMaxValue() const {return top;}
      inline double getMinValue() const {return bot;}
      inline virtual double activatorDerivative(double x) const override {return (x==bot || x==top)?0:(1-EPSILON);}
      inline virtual LinearNeuron* clone() const override
      {
        LinearNeuron* ret= new LinearNeuron(sz,weights,cutoff,learnrate,top,bot);
        ret->addBias(bias);
        ret->setBiasWeight(this->getBiasWeight());
        return ret;
      }
      inline virtual LinearNeuron* reinit(int s,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        LinearNeuron* ret=new LinearNeuron(s,v,lr);
        ret->setMaxValue(top);
        ret->setMinValue(bot);
        return ret;
      }
      inline virtual LinearNeuron* reinit(int s,GeneralVector<double> extw,double v=0,double lr=0) const override
      {
        if(!lr) lr=this->learnrate;
        LinearNeuron* ret=new LinearNeuron(s,extw,v,lr);
        ret->setMaxValue(top);
        ret->setMinValue(bot);
        return ret;
      }
    };
    static std::vector<Neuron*> NEURONTYPES;
    class NeuronLayer : public flex::Cloneable
    {
    protected:
      std::vector<Neuron*> neurons;
      std::vector<std::vector<int>> insplit;
      int seq=0,inp=0;
    public:
      NeuronLayer() {}
      NeuronLayer(int in) {inp=in;}
      NeuronLayer(int nc,Neuron* sample,double lrw=-1,double urw=1,bool bias=true,double b=1)
      {
        for(int i=0;i<nc;i++)
        {
          neurons.push_back(sample->clone());
          neurons[i]->randomizeWeights(lrw,urw);
          if(bias) neurons[i]->addBias(b);
        }
      }


      inline int getNeuronCount() const {return neurons.size();}
      inline void addNeuron(Neuron* n) {neurons.push_back(n);}
      inline const std::vector<Neuron*>& getNeurons() const {return neurons;}
      inline void refresh() {for(Neuron* n : neurons) n->refresh();}

      GeneralVector<double> evaluate(const GeneralVector<double>& in) const
      {
        std::vector<double> result;
        std::vector<GeneralVector<double>> inputs=getInputs(in);
        for(int i=0;i<inputs.size();i++) result.push_back(neurons[i]->operator()(inputs[i]));
        return result;
      }
      virtual std::vector<GeneralVector<double>> getInputs(const GeneralVector<double>& in) const //Use this to split up selected parts of the input to pass to different neurons (values can repeat between neurons)
      {
        //cout << inp << "\t"<<neurons.size()<<"\n";
        std::vector<GeneralVector<double>> ret;
        for(int i=0;i<neurons.size();i++) ret.push_back(in);
        return ret;
      }
      virtual std::vector<std::vector<int>> calculateInputSplit(int insize=0) //This gives (with weights) how the inputs are split
      {
        if(!insize) insize=(inp)?inp:neurons.size();
        insplit=std::vector<std::vector<int>>();
        std::vector<int> v; for(int j=0;j<insize;j++) v.push_back(j);
        for(int i=0;i<neurons.size();i++) insplit.push_back(v);
        return insplit;
      }
      inline const std::vector<std::vector<int>>& getInputSplit() const {return insplit;}
      inline void setInputSplit(const std::vector<std::vector<int>>& v) {insplit=v;}

      void sreyasBackPropagateAdjust(std::vector<NeuronLayer*> nl,const std::vector<GeneralVector<double>>& dataline,int ind, const GeneralVector<double>& exp)
      {
        //GeneralVector<double> error=dataline[ind+1]-exp;
        //cout << "Neurons\n";
        //for(int i=0;i<neurons.size();i++) cout << i << "\t"<<(long int)(neurons[i])<<"\n";
        std::vector<GeneralVector<double>> nodeins=getInputs(dataline[ind]);
        std::vector<double> errs;
        //cout <<"Sizes: "<< exp.size() << "\t"<<neurons.size() << "\n";
        for(int i=0;i<exp.size();i++) {errs.push_back(neurons[i]->adjust(nodeins[i],exp(i)));}
        //double err=neurons[seq]->adjust(nodeins[seq],exp(seq));
        if(ind>0)
        {
          //std::vector<std::vector<int>> insplit=getInputSplit();
          if(!insplit.size()) throw UnPackedNeuralNetworkException();
          GeneralVector<double> nexp=dataline[ind];
          for(int i=0;i<neurons.size();i++)
          {
            int K=0;
            for(int j=0;j<insplit[i].size();j++)
              nexp(j)+=errs[i]/neurons[i]->getWeights()(K++);
          }
          nl[ind-1]->sreyasBackPropagateAdjust(nl,dataline,ind-1,nexp);
        }
        seq++; seq%=neurons.size();
      }
      void derivativeBackPropagateAdjust(std::vector<NeuronLayer*> nl,const std::vector<GeneralVector<double>>& dataline,int ind, GeneralVector<double> errv)
      {
        for(int i=0;i<neurons.size();i++)
        {
          double der=neurons[i]->activatorDerivative(dataline[ind+1](i))+EPSILON;
          errv(i)*=(der/neurons[i]->getWeights().size());
          //neurons[i]->onlyAdjust(dataline[ind],errv(i));
        }
        /*cout <<errv.transpose()<<"\n";
        cout <<dataline[ind].transpose()<<"\n\n";*/
        if(ind>0)
        {
          if(!insplit.size()) throw UnPackedNeuralNetworkException();
          GeneralVector<double> nerrv(nl[ind-1]->getNeuronCount());
          for(int i=0;i<neurons.size();i++)
          {
            int K=0;
            for(int j=0;j<insplit[i].size();j++) nerrv(j)+=errv(i)*neurons[i]->getWeights()(K++);
            //for(int j=0;j<insplit[i].size();j++) nerrv(j)+=errv(i)*neurons[i]->getLearnRate()*neurons[i]->getWeights()(K++);
          }
          //cout << "\n\n";
          nl[ind-1]->derivativeBackPropagateAdjust(nl,dataline,ind-1,nerrv);
        }
        for(int i=0;i<neurons.size();i++) neurons[i]->onlyAdjust(dataline[ind],errv(i));
      }

      inline NeuronLayer* clone() const override
      {
        NeuronLayer* ret=new NeuronLayer();
        for(int i=0;i<neurons.size();i++) ret->addNeuron(neurons[i]->clone());
        ret->insplit=insplit;
        ret->inp=inp;
        return ret;
      }
      inline void scaleLearnRates(double v=1) {for(Neuron* n : neurons) n->scaleLearnRates(v);}
      //Iterators
      inline auto begin() {return neurons.begin();}
      inline auto end() {return neurons.end();}
    };
    class SplitDataNeuronLayer : public NeuronLayer
    {
    public:
      SplitDataNeuronLayer() : NeuronLayer() {}
      SplitDataNeuronLayer(int n) : NeuronLayer(n) {}

      std::vector<std::vector<int>> calculateInputSplit(int insize=0) override
      {
        int K=0;
        if(!insize) insize=(inp)?inp:neurons.size();
        insplit=std::vector<std::vector<int>>();
        for(int i=0;i<neurons.size();i++)
        {
          std::vector<int> tv;
          for(int j=0;j<neurons[i]->getInputCount();j++) tv.push_back(K++);
          insplit.push_back(tv);
        }
        return insplit;
      }
      std::vector<GeneralVector<double>> getInputs(const GeneralVector<double>& in) const override //Use this to split up selected parts of the input to pass to different neurons (values can repeat between neurons)
      {
        //cout << inp << "\t"<<neurons.size()<<"\n";
        std::vector<GeneralVector<double>> ret;
        int K=0;
        for(int i=0;i<neurons.size();i++)
        {
          GeneralVector<double> t(neurons[i]->getInputCount());
          for(int j=0;j<neurons[i]->getInputCount();j++) t(j)=in(K++);
          ret.push_back(t);
        }
        return ret;
      }
    };
    class NeuralNetwork : public PreservablePredictor
    {
      std::vector<NeuronLayer*> layers;
    public:
      NeuralNetwork() {}
      NeuralNetwork(const std::vector<NeuronLayer*>& l)  {layers=l;}
      NeuralNetwork(std::istream& fn,std::ostream& logf) {load(fn,logf);}
      NeuralNetwork(const std::vector<int> layerdat,int ins,Neuron* ntype=nullptr,bool autopack=true,double lwl=-1,double uwl=1,double E=1.5,double S=1)
      {
        double LR=S;
        if(!ntype) {cout << "Using default neurons!\n"; ntype=new Neuron(10);}
        int i=0;
        cout << "Building NN with "<<layerdat.size()<<" layers with neuron counts:\n";
        while(i<(int)layerdat.size())
        {
          int c=(i==0)?ins:layerdat[i-1];
          cout << layerdat[i]<<" taking input of size: "<<c << "\n";
          layers.push_back(new NeuronLayer(layerdat[i],ntype->reinit(c,0,LR),lwl,uwl));
          LR/=E;
          i++;
        }
        if(autopack) {this->pack(ins); cout << "Packed!\n";}
        cout << "NN ready\n";
      }


      inline const std::vector<NeuronLayer*>& getLayers() const {return layers;}
      //inline std::vector<NeuronLayer*>& getLayers() {return layers;}
      inline int getLayerCount() const {return layers.size();}
      inline void addLayer(NeuronLayer* l) {layers.push_back(l);}
      inline void train(const GeneralVector<double>& data,const GeneralVector<double>& exp) {derivativeBackPropagateAdjust(data,exp);}
      inline void pack(int s)
      {
        for(int i=1;i<layers.size();i++) {layers[i]->calculateInputSplit(layers[i-1]->getNeuronCount());}
        layers[0]->calculateInputSplit(s);
      }

      inline GeneralVector<double> evaluate(const GeneralVector<double>& input) const {return evaluate(input,layers.size()-1);}
      GeneralVector<double> evaluate(const GeneralVector<double>& input,int stop) const
      {
        GeneralVector<double> in=input;
        for(int i=0;i<=stop;i++) in=layers[i]->evaluate(in);
        return in;
      }

      void sreyasBackPropagateAdjust(const GeneralVector<double>& data,const GeneralVector<double>& exp)
      {
        std::vector<GeneralVector<double>> dataline(1,data);
        for(int i=0;i<layers.size();i++) dataline.push_back(layers[i]->evaluate(dataline[i]));
        //cout << "Received: "<<exp.transpose()<<" as correction for predicted: "<<dataline[layers.size()].transpose()<<"\n";
        if(exp==dataline[layers.size()]) return;
        layers[layers.size()-1]->sreyasBackPropagateAdjust(layers,dataline,layers.size()-1,exp);
      }
      void derivativeBackPropagateAdjust(const GeneralVector<double>& data,const GeneralVector<double>& exp)
      {
        std::vector<GeneralVector<double>> dataline(1,data);
        for(int i=0;i<layers.size();i++) dataline.push_back(layers[i]->evaluate(dataline[i]));
        GeneralVector<double> errv=exp-dataline[layers.size()];
        double err=errv.norm2()/2.0;
        if(!err) return;
        /*cout << "Expected: "<<exp.transpose()<<"\n";
        cout << "Got: "<<dataline[layers.size()].transpose()<<"\n";*/
        layers[layers.size()-1]->derivativeBackPropagateAdjust(layers,dataline,layers.size()-1,errv);
        /*cout << "Corrected!!!\n";
        std::vector<GeneralVector<double>> dataline2(1,data);
        for(int i=0;i<layers.size();i++) dataline2.push_back(layers[i]->evaluate(dataline2[i]));
        cout << "Expected: "<<exp.transpose()<<"\n";
        cout << "Got: "<<dataline2[layers.size()].transpose()<<"\n\n\n---------\n\n";*/
      }
      inline GeneralVector<double> operator()(const GeneralVector<double>& in)const {return evaluate(in);}
      inline void save(std::ostream& os) const override {writeNNData(os);}
      inline void load(std::istream& fn,std::ostream& logf=std::cout) override
      {
        layers=std::vector<NeuronLayer*>();
        int lc; fn >> lc;
        logf << "Structure has "<<lc<<" layers\n";
        NeuronLayer* lr=nullptr;
        stringstream ss;
        for(int ln=0;ln<lc;ln++)
        {
          lr=new NeuronLayer();
          logf << "Loading layer... "<<(ln+1)<<"\n";
          int nc; fn >> nc;
          std::string inl;
          logf << "Layer has "<<nc<<" neurons\n"; getline(fn,inl);
          std::vector<std::vector<int>> splv;
          for(int cr1=0;cr1<nc;cr1++)
          {
            getline(fn,inl);
            ss=stringstream(inl);
            std::vector<int> sv; int sp;
            while(true)
            {
              ss >> sp;
              if(ss.eof()) break;
              sv.push_back(sp);
            }
            splv.push_back(sv);
          }
          logf << "Loaded insplits\n";
          for(int cr1=0;cr1<nc;cr1++)
          {
            int wc,tp;  double lrt;
            fn >> tp; fn >> wc; fn>>lrt;
            logf << "Neuron of type "<<tp<<" has "<<wc<<" weights\n";
            std::vector<double> wts;
            double n;
            for(int cr2=0;cr2<wc;cr2++)
            {
              fn >> n;
              wts.push_back(n);
            }
            Neuron* nn;
            switch(tp)
            {
              case 1:
                nn=new SigmoidNeuron(wts.size(),wts,0,lrt);
                break;
              case 2:
                nn=new LinearNeuron(wts.size(),wts,0,lrt);
                break;
              default:
                nn=new Neuron(wts.size(),wts,0,lrt);
            }
            double bs,bw; fn >> bs; fn >> bw;
            nn->addBias(bs); nn->setBiasWeight(bw);
            lr->addNeuron(nn);
          }
          lr->setInputSplit(splv);
          layers.push_back(lr);
          logf << "Done loading layer\n";
        }
      }
      void writeNNData(std::ostream& os) const
      {
        os << layers.size() << "\n";
        for(int i=0;i<layers.size();i++)
        {
          os << " " << layers[i]->getNeuronCount() << "\n";
          auto is=layers[i]->getInputSplit();
          for(int j=0;j<layers[i]->getNeuronCount();j++)
          {
            for(auto& v : is[j]) os << v << " ";
            os << "\n";
          }
          for(Neuron* n : layers[i]->getNeurons())
          {
            os <<"  "<<n->type<<" "<<(n->getWeights().getSize()-1) <<" "<<n->getLearnRate()<< "\n";
            for(int i=0;i<(n->getWeights().getSize()-1);i++) os << "  " << n->getWeights()(i) << "\n";
            os << "  "<<n->getBias()<<" "<<n->getBiasWeight()<<"\n";
          }
        }
      }
      virtual NeuralNetwork* clone() const override
      {
        std::vector<NeuronLayer*> newlays;
        for(int i=0;i<layers.size();i++) newlays.push_back(layers[i]->clone());
        return new NeuralNetwork(newlays);
      }
      inline void refresh() override {for(NeuronLayer* nl : layers) nl->refresh();}
      inline void scaleLearnRates(double v=1) override {for(NeuronLayer* l : layers) l->scaleLearnRates(v);}
    };
    class SingleNeuronModel : public PreservablePredictor
    {
      Neuron* heldNeuron;
    public:
      SingleNeuronModel(Neuron* src) {heldNeuron=src;}

      inline GeneralVector<double> evaluate(const GeneralVector<double>& v) const override {return std::vector<double>(1,heldNeuron->evaluate(v));}
      inline Neuron* getNeuron() {return heldNeuron;}
      inline const Neuron* getNeuron() const {return heldNeuron;}

      inline void train(const GeneralVector<double>& in,const GeneralVector<double>& out) override {heldNeuron->adjust(in,out(0));}
      inline void refresh() override {heldNeuron->refresh();}
      inline SingleNeuronModel* clone() const override {return new SingleNeuronModel(heldNeuron->clone());}
      void save(std::ostream& os) const override {}
      void load(std::istream& is,std::ostream& log) override {}
      inline void scaleLearnRates(double v=1) override {heldNeuron->scaleLearnRates(v);}
      //inline double evaluate(const GeneralVector<double>& v) const override {return heldNeuron->evaluate(v);}
    };
  }
  static void trainModelByData(Predictor& model,const DataSet& data,const DataSet& labels,const int repeats=1,const int rf=500)
  {
    cout << "Started training...\n";
    for(int r=0;r<repeats;r++)
    {
      for(int i=0;i<data.size();i++)
      {
        model.train(data[i],labels[i]);
        if((i+1)%rf==0) cout << (i+1) << " trained\n";
        //if((i+1)%5==0) exit(1);
      }
      cout << "Completed "<<(r+1)<<" of "<<repeats<<" iterations\n";
    }
    cout << "Done\n";
  }
  static double testModelForError(const Predictor& model,const DataSet& data,const DataSet& labels,double (*errf)(const GeneralVector<double>& v1,const GeneralVector<double>& v2)=euclideanDistanceSquare,const int rf=1000)
  {
    cout << "Started testing...\n";
    double err=0;
    double tev=0;
    for(int i=0;i<data.size();i++)
    {
      GeneralVector<double> res=model.evaluate(data[i]);
      tev=errf(res,labels[i]);
      err+=tev;
      if(tev>1) cout << "tev>1: "<<tev<<"\n";
      if((i+1)%rf==0) cout << (i+1) << " tested\n";
    }
    cout << err << " of "<<data.size()<<" samples\n";
    cout << "Testing completed\n";
    return err/data.size();
  }

  namespace enhance
  {
    class Bagging : public PreservablePredictor
    {
    protected:
      std::vector<Predictor*> models;
      std::vector<double> weights;
      std::vector<PreservablePredictor*> loadorder;
    public:
      Bagging() {}
      Bagging(std::vector<Predictor*> preds) {models=preds;}

      GeneralVector<double> evaluate(const GeneralVector<double>& in) const override
      {
        if(!models.size()) throw NoModelsForBaggingException();
        GeneralVector<double> ret=(models[0]->evaluate(in))*weights[0];
        for(int i=1;i<models.size();i++) ret=ret+(models[i]->evaluate(in)*weights[i]);
        return ret;
      }

      inline void normalizeWeights()
      {
        double s=0; for(double& d : weights) s+=d;
        for(double& d : weights) d/=s;
      }

      inline void setLoadOrder(std::vector<PreservablePredictor*> preds) {loadorder=preds;}
      inline void save(std::ostream& os) const override
      {
        os << models.size() << "\n";
        for(Predictor* p : models) ((PreservablePredictor*)p)->save(os);
        for(const double& d : weights) os << d << " ";
        os << "\n";
      }
      inline void load(std::istream& is,std::ostream& log=std::cout) override
      {
        int mc; is >> mc;
        cout << mc << " models to load\n";
        if(loadorder.size()<mc) throw UndeterminedLoadTypesException();
        for(int i=0;i<mc;i++)
        {
          PreservablePredictor* npp=(PreservablePredictor*)(loadorder[i]->clone());
          npp->load(is,log);
          models.push_back(npp);
        }
        cout << "Loading weights\n";
        std::vector<double> wts;
        for(int i=0;i<models.size();i++)
        {
          double w; is >> w;
          wts.push_back(w);
        }
        cout << "Done\n";
        weights=wts;
      }
      inline void scaleLearnRates(double v=1) override {for(auto* p : models) p->scaleLearnRates(v);}

      /*void buildForest(Predictor* pr,const std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>>& traindata,int bcount=4,int splitc=0)
      {
        const std::vector<std::vector<double>>& ins=get<0>(traindata);
        const std::vector<std::vector<double>>& outs=get<1>(traindata);
        if(!splitc) splitc=ins.size()/(2*bcount);
        if(!splitc) throw InsufficientDataToTrainException();
        std::vector<double> corr;
        for(int m=0;m<bcount;m++)
        {
          Predictor* p=pr->clone();
          for(int i=0;i<splitc;i++)
          {
            int ind=randgen::throwarandompoint(0,ins.size());
            p->train(ins[ind],outs[ind]);
            if((i+1)%500==0) cout << (i+1) << " trained out of "<<splitc<<"\n";
          }
          int cc=0;
          for(int i=0;i<3*splitc;i++)
          {
            int ind=randgen::throwarandompoint(0,ins.size());
            if(testParameter(p->evaluate(ins[ind]),outs[ind])) cc++;
          }
          cout << "Model "<<m<<" identified "<<cc<<" correctly out of "<<(3*splitc)<<" tried.\n";
          models.push_back(p);
          corr.push_back(cc);
        }
        weights=corr;
        normalizeWeights();
      }*/
      virtual void templateModel(Predictor* pr,const std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>>& traindata,int bcount=4,int splitc=0)
      {
        int tvar;
        const std::vector<std::vector<double>>& ins=get<0>(traindata);
        const std::vector<std::vector<double>>& outs=get<1>(traindata);
        if(!splitc) splitc=ins.size()/(2*bcount);
        if(!splitc) throw InsufficientDataToTrainException();
        std::vector<double> corr;
        cout << "Loading model for Bagging with "<<bcount<<" branches\n";
        for(int m=0;m<bcount;m++)
        {
          Predictor* p=pr->clone();
          for(int i=0;i<splitc;i++)
          {
            int ind=randgen::throwarandompoint(0,ins.size());
            p->train(ins[ind],outs[ind]);
            if((i+1)%500==0) cout << (i+1) << " trained out of "<<splitc<<"\n";
          }
          cout << "Trained on "<<splitc<<" samples\n";
          cout << "Started testing\n";
          int cc=0;
          for(int i=0;i<2*splitc;i++)
          {
            int ind=randgen::throwarandompoint(0,ins.size());
            if(testParameter(p->evaluate(ins[ind]),outs[ind])) cc++;
          }
          cout << "Model "<<(m+1)<<" identified "<<cc<<" correctly out of "<<(2*splitc)<<" tried.\n";
          models.push_back(p);
          corr.push_back(cc);
        }
        weights=corr;
        normalizeWeights();
        cout << weights << "\n";
      }
      virtual bool testParameter(std::vector<double> out,std::vector<double> exp)
      {
        double m=-1; int ind=-1,expind=-1;
        for(int i=0;i<out.size();i++)
        {
          if(m<out[i]) {m=out[i]; ind=i;}
        }
        for(int j=0;j<exp.size();j++)
        {
          if(exp[j]>0.6) {expind=j; break;}
        }
        return(expind==ind);
      }
      virtual Bagging* clone() const override
      {
        Bagging* ret=new Bagging();
        ret->weights=weights;
        for(Predictor* p : models) ret->models.push_back(p->clone());
        ret->loadorder=loadorder;
        return ret;
      }
      const Predictor* getModel(int x) const {return models[x];}
      void train(const GeneralVector<double>& in,const GeneralVector<double>& exp) {for(Predictor* p : models) p->train(in,exp);}
      void refresh() override {for(Predictor* p : models) p->refresh(); /*Refresh weights?*/}
    };
    class BoostedBagging : public Bagging
    {
    public:
      BoostedBagging() : Bagging() {}
      BoostedBagging(std::vector<Predictor*> preds) : Bagging(preds) {}

      virtual void templateModel(Predictor* pr,const std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>>& traindata,int bcount=4,double pe=0.01)
      {
        int tvar;
        const std::vector<std::vector<double>>& ins=get<0>(traindata);
        const std::vector<std::vector<double>>& outs=get<1>(traindata);
        std::vector<double> corr;
        cout << "Loading model for boosted bagging with "<<bcount<<" branches\n";
        std::vector<int> incorrs; for(int i=0;i<ins.size();i++) incorrs.push_back(i);
        int tott;
        for(int m=0;m<bcount;m++)
        {
          Predictor* p=pr->clone(); p->refresh();
          tott=0;
          while(tott<ins.size()*(1-pe))
          {
            cout << "Training (further) as "<<tott<<" < "<<(int)(ins.size()*(1-pe))<<"\n";
            for(int i=0;i<incorrs.size();i++)
            {
              p->train(ins[incorrs[i]],outs[incorrs[i]]);
              if((i+1)%500==0) cout << (i+1) << " trained out of "<<incorrs.size()<<"\n";
              tott++;
            }
          }
          cout << "Trained on "<<tott<<" samples\n";
          cout << "Full data training (repeat)\n";
          for(int i=0;i<ins.size();i++)
          {
            p->train(ins[i],outs[i]);
            if((i+1)%500==0) cout << (i+1) << " trained out of "<<ins.size()<<"\n";
          }
          cout << "Started testing\n";
          int cc=0;
          incorrs=std::vector<int>();
          for(int i=0;i<ins.size();i++)
          {
            if(testParameter(p->evaluate(ins[i]),outs[i])) cc++;
            else incorrs.push_back(i);
          }
          cout << "Model "<<(m+1)<<" identified "<<cc<<" correctly out of "<<(ins.size())<<" tried.\n";
          models.push_back(p);
          corr.push_back(cc);
        }
        weights=corr;
        normalizeWeights();
        cout << weights << "\n";
      }

      virtual BoostedBagging* clone() const override
      {
        BoostedBagging* ret=new BoostedBagging();
        ret->weights=weights;
        for(Predictor* p : models) ret->models.push_back(p->clone());
        ret->loadorder=loadorder;
        return ret;
      }
    };

    /*static std::pair<std::vector<double>,std::vector<double>> adjustParameter(double& par,Predictor& model,const DataSet& traindata,const DataSet& trainlabels,const DataSet& testdata,const DataSet& testlabels,double (*errf)(const GeneralVector<double>& v1,const GeneralVector<double>& v2)=euclideanDistanceSquare,double ecut=0.01,bool minimize=true,double D=0.01)
    {
      static std::vector<double> vallist,errlist;
      vallist.push_back(par);
      double G=ecut+1;
      while(G>ecut)
      {
        double err=0;

      }
    }*/
  }
  namespace support
  {
    static std::vector<GeneralVector<double>> extractFeatures(const std::vector<GeneralVector<double>>& rawdata,std::vector<int> features)
    {
      std::vector<GeneralVector<double>> ret;
      for(const GeneralVector<double>& rv : rawdata)
      {
        std::vector<double> v;
        for(int i : features) v.push_back(rv(i));
        ret.push_back(v);
      }
      return ret;
    }
    static std::vector<GeneralVector<double>> removeFeatures(const std::vector<GeneralVector<double>>& rawdata,std::vector<int> features)
    {
      std::vector<int> nfeat;
      for(int i=0;i<rawdata[0].getSize();i++)
      {
        if(contains(features,i)) continue;
        nfeat.push_back(i);
      }
      return extractFeatures(rawdata,nfeat);
    }
  }
  namespace enhance
  {
    class FeatureOptimization // : public PreservablePredictor
    {
      std::vector<PreservablePredictor*> preds;
      std::vector<std::vector<int>> insplits;
      std::vector<double> weights;
    public:
      FeatureOptimization(PreservablePredictor* samplepred,int parc,const std::vector<GeneralVector<double>>& data,std::vector<GeneralVector<double>> values,int rep=1,int bcount=4,double frac=0.8,const std::vector<int>& statf=std::vector<int>())
      {
        std::vector<int> nf;
        int featcount=data[0].getSize();
        assert(featcount>parc);
        int ind;
        for(int i=0;i<bcount;i++)
        {
          nf=statf;
          while(nf.size()<parc)
          {
            int rv=randgen::throwarandompoint(0,featcount);
            if(!contains(nf,rv)) nf.push_back(rv);
          }
          PreservablePredictor* np=(PreservablePredictor*)(samplepred->clone()); np->refresh();
          std::vector<GeneralVector<double>> exf=support::extractFeatures(data,nf);
          for(int r=0;r<rep;r++)
          {
            for(int t=0;t<exf.size()*frac;t++)
            {
              ind=randgen::throwarandompoint(exf.size());
              np->train(exf[ind],values[ind]);
            }
          }
          double e=0; int K=0;
          for(int t=exf.size()*frac;t<exf.size();t++)
          {
            ind=randgen::throwarandompoint(exf.size());
            auto pv=np->evaluate(exf[ind]);
            e+=GeneralVector<double>(values[ind]-pv).norm()/values[ind].norm();
            K++;
          }
          e/=K;
          cout << "Model "<<(i+1)<<" predicted with error average of "<<(e*100)<<"%\n";
          weights.push_back(1-e);
          preds.push_back(np);
          insplits.push_back(nf);
        }
        double v=0; for(double w : weights) v+=w;
        for(double& w : weights) w/=v;
        cout << "Weights: "<< weights<<"\n";
      }
      /*FeatureOptimization(PreservablePredictor* samplepred,int parc,const std::vector<GeneralVector<double>>& data,std::vector<double> values,int rep=1,int bcount=4,double frac=0.8)
      {

        FeatureOptimization(samplepred,parc,data,nvvec,rep,bcount,frac);
      }*/

      GeneralVector<double> evaluate(const GeneralVector<double>& in) const //override
      {
        if(!preds.size()) throw NoModelsForPredictorException();
        std::vector<double> fext; for(int i=0;i<insplits[0].size();i++) fext.push_back(in(i));
        GeneralVector<double> ret=(preds[0]->evaluate(fext))*weights[0];
        for(int i=1;i<preds.size();i++)
        {
          fext=std::vector<double>();
          for(int j=0;j<insplits[i].size();j++) fext.push_back(in(insplits[i][j]));
          auto pv=preds[i]->evaluate(fext);
          cout << "Model "<<(i+1)<<" predicts "<<pv.transpose()<<" with insplits:\t"<<insplits[i]<<"\n";
          ret=ret+pv*weights[i];
        }
        return ret;
      }

      inline void normalizeWeights()
      {
        double s=0; for(double& d : weights) s+=d;
        for(double& d : weights) d/=s;
      }
      inline void refresh() {for(PreservablePredictor* p : preds) p->refresh();}
    };
  }
}
#endif
