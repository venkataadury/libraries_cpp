#ifndef QUEUING_H
#define QUEUING_H 1
#include "maths/stats.h"
#include <queue>
#define TIMES(n) for(int TV=0;TV<(n);TV++)
#define E 2.718281828
namespace queuing
{
  class KendallQueue;
  class EndOfQueueException : public exception {};
  class ServerFullException : public exception {};
  class Customer
  {
  public:
    long long intime=0,outqtime=-1,outtime=-1;
    int job=0;
    Customer() {}
    Customer(long long t,int j=0) {intime=t;job=j;}

    inline int getJobID() const {return job;}
    inline long long getInTime() const {return intime;}
    inline long long getQueueTime() const {if(outqtime!=-1) return outqtime-intime; else return -1;}
    inline long long getTotalTime() const {if(outtime!=-1) return outtime-intime; else return -1;}
    inline long long getServiceTime() const {if(outtime!=-1) return outtime-outqtime; else return -1;}
  };
  class QueueLine
  {
    std::vector<Customer> q;
    double lambda=0.1;
    int droppedArrivals=0;
    double qlim=1.0/0.0;
  public:
    QueueLine(double r,double lim=1.0/0.0) {lambda=r; qlim=lim;}
    /*Using default constructors*/

    inline int size() const {return q.size();}
    inline double getArrivalRate() const {return lambda;}
    inline void setArrivalRate(double r) {lambda=r;}
    inline int getDroppedArrivals() const {return droppedArrivals;}

    void simulate(const std::vector<double>& jobs,const KendallQueue& kq);

    virtual bool isEmpty() const {return q.empty();}
    virtual Customer getNextCustomer() //Queue Discipline
    {
      if(!q.empty())
      {
        Customer ret=q.front();
        q.erase(q.begin());
        return ret;
      }
      else throw EndOfQueueException();
    }
  };
  class Server
  {
    double relR=1;
    double timeunit;
    Customer c;
    bool hascust=false;
    std::vector<int> offeredjobs;
  public:
    Server(double TU,double r=1) {timeunit=TU; if(TU!=1) relR=(1-::pow(E,-r*TU))/(1-::pow(E,-TU));}

    inline bool hasCustomer() const {return hascust;}
    inline Customer* getCustomerReference() {return &c;}
    inline virtual void simulate(const std::vector<double>& jobrates,KendallQueue& kq);
    inline void takeCustomer(const Customer& cs,const KendallQueue& kq);
  };
  class Metric
  {
    std::string name;
  public:
    Metric(const std::string& mn="metric") {name=mn;}

    virtual double measure(const KendallQueue& kq) const {return 0;}
  };
  class KendallQueue
  {
    long long time=0;
    std::vector<double> services; //Probabilities of various services
    std::vector<double> servicerates; //Rates of various services
    std::vector<QueueLine> ql;
    std::vector<Server> servers;
  public:
    std::vector<std::vector<double>> measurements;
    std::vector<const Metric*> metrics;
    std::vector<std::vector<double>> custnettimes;
    long int totcust=0,syscust=0,tsyscust;
    double timeunit=1,lim=1.0/0.0;

  public:
    KendallQueue(double l,double m,int sc=1,int qc=1,double li=(double)1.0/0,double qlim=(double)1.0/0)
    {
      ql=std::vector<QueueLine>();
      TIMES(qc) ql.push_back(QueueLine(l,qlim));
      servers=std::vector<Server>();
      TIMES(sc) servers.push_back(Server(timeunit));
      servicerates.push_back(m);
      services.push_back(1);
      lim=li;
    }
    inline long long getTime() const {return time;}

    void simulate(long long tL,bool reset=true)
    {
      if(reset) totcust=0;
      while(time<tL)
      {
        tsyscust=0;
        custnettimes=std::vector<std::vector<double>>();
        TIMES(3) custnettimes.push_back(std::vector<double>());
        for(QueueLine& q : ql) {tsyscust+=q.size(); q.simulate(services,*this);}
        for(Server& s : servers)
        {
          Customer* nc=nullptr;
          if(s.hasCustomer()) {nc=s.getCustomerReference(); tsyscust++;}
          s.simulate(servicerates,*this);
          //cout << s.hasCustomer() << "\n";
          if(!s.hasCustomer())
          {
            if(nc)
            {
              custnettimes[0].push_back(nc->getQueueTime());
              custnettimes[1].push_back(nc->getServiceTime());
              custnettimes[2].push_back(nc->getTotalTime());
              totcust++;
            }
            std::vector<QueueLine*> freeqs;
            for(QueueLine& q : ql)
            {
              if(q.isEmpty()) continue;
              freeqs.push_back(&q);
            }
            if(freeqs.size()) s.takeCustomer(randgen::choice(freeqs)->getNextCustomer(),*this);
          }
        }
        syscust=tsyscust;
        for(int i=0;i<metrics.size();i++) measurements[i].push_back(metrics[i]->measure(*this));
        time++;
        if(time%(tL/100)==0) cout << time << "\n";
        //cout << time <<"\t"<< totcust << " "<<syscust<<"\n";
      }
    }
    inline bool verify(int err=4) const
    {
      return ((totcust-syscust)>=err && ((double)(syscust-(long)ql.size()-(long)servers.size())/totcust<0.05));
    }

    inline const std::vector<Server>& getServers() const {return servers;}
    inline const std::vector<QueueLine>& getQueues() const {return ql;}

    void addMetric(const Metric& met)
    {
      metrics.push_back(&met);
      measurements.push_back(std::vector<double>());
    }
  };

  namespace metrics
  {
    class QueueTimeMetric : public Metric //This is a sum of queue times at the for customers leaving at the given frame. To get the average, sum over all measurements in KendallQueue and divide by total customer count
    {
    public:
      QueueTimeMetric() : Metric("queue-time") {}

      double measure(const KendallQueue& kq) const override
      {
        double s=0;
        for(double wt : kq.custnettimes[0]) {if(wt!=-1) s+=wt;}
        return s;
      }
    } QUEUETIME;
    class ServiceTimeMetric : public Metric //This is a sum of service times at the for customers leaving at the given frame. To get the average, sum over all measurements in KendallQueue and divide by total customer count
    {
    public:
      ServiceTimeMetric() : Metric("service-time") {}

      double measure(const KendallQueue& kq) const override
      {
        double s=0;
        for(double wt : kq.custnettimes[1]) {if(wt!=-1) s+=wt;}
        return s;
      }
    } SERVICETIME;
    class TotalTimeMetric : public Metric //This is a sum of total times at the for customers leaving at the given frame. To get the average, sum over all measurements in KendallQueue and divide by total customer count
    {
    public:
      TotalTimeMetric() : Metric("sojourn-time") {}

      double measure(const KendallQueue& kq) const override
      {
        double s=0;
        for(double wt : kq.custnettimes[2]) {if(wt!=-1) s+=wt;}
        return s;
      }
    } TOTALTIME;
    class CustomerCountMetric : public Metric
    {
    public:
      CustomerCountMetric() : Metric("customer-count") {}

      inline double measure(const KendallQueue& kq) const override {return kq.syscust;}
    } CUSTCOUNT;
    class UtilizationMetric : public Metric
    {
    public:
      UtilizationMetric() : Metric("utilization") {}

      inline double measure(const KendallQueue& kq) const override //This is a utilization at that instant of the simulation. To get the expectation, you have to take the average in KendallQueue for all time
      {
        int occ=0;
        for(const Server& s : kq.getServers()) {if(s.hasCustomer()) occ++;}
        return (double)occ/kq.getServers().size();
      }
    } UTILIZATION;
    class DropMetric : public Metric
    {
    public:
      DropMetric() : Metric("dropped-arrivals") {}

      inline double measure(const KendallQueue& kq) const override //This records the total dropouts until time 't' at each step
      {
        int occ=0;
        for(const auto& s : kq.getQueues()) {occ+=s.getDroppedArrivals();}
        return occ;
      }
    } DROPS;
  }
}
void queuing::QueueLine::simulate(const std::vector<double>& jobs,const KendallQueue& kq)
{
  double p=randgen::throwarandompoint(0.0,1.0);
  if(p>lambda*kq.timeunit) return;
  if(kq.syscust>=kq.lim || q.size()>=qlim) {droppedArrivals++; return;}
  p=randgen::throwarandompoint(0.0,1.0);
  int i=0;
  for(;i<jobs.size();i++)
  {
    p-=jobs[i];
    if(p<0) break;
  }
  //cout << "Entry into system\n";
  q.push_back(Customer(kq.getTime(),i));
}
inline void queuing::Server::simulate(const std::vector<double>& jobrates,KendallQueue& kq)
{
  if(!hascust) return;
  double p=randgen::throwarandompoint(0,1);
  if(p<jobrates[c.getJobID()]*relR)
  {
    //cout << p << " vs " << jobrates[c.getJobID()] << "\n";
    //cout << "Exit from system\n";
    c.outtime=kq.getTime();
    hascust=false;
  }
}
inline void queuing::Server::takeCustomer(const Customer& cs,const KendallQueue& kq) {c=cs; c.outqtime=kq.getTime(); hascust=true;}
#endif
