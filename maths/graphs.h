#ifndef INCLUDED_GRAPHS
#define INCLUDED_GRAPHS
#include "maths/maths.h"
namespace graphs
{
  class EdgeDoesNotConnectVertexException : public std::exception {};
  class BadEdgeDesignException : public exception {};
  class VertexNotInGraphException : public exception {};
  class AlgorithmFailStateException : public std::exception {};
  class NoPathPossibleException : public AlgorithmFailStateException {};
  class NoCyclePossibleException : public AlgorithmFailStateException {};
  class NoEulerPathPossibleException : public NoPathPossibleException {};
  class InvalidFlowException : public std::exception
  {
    std::string errmes;
  public:
    InvalidFlowException() : InvalidFlowException("Invalid flow encountered in flow algorithm") {}
    InvalidFlowException(std::string m) {errmes=m; cerr << errmes <<"\n";}

    virtual const char* what() throw() {return errmes.c_str();}
  };

  template<class T,class U=T> class Edge;
  template<class T,class U=T> class Vertex;
  template <class T,class U=T> class EdgeEnd
  {
  public:
    Vertex<U,T>* v;
    bool enter=true,leave=true;
  public:
    bool flag=false;

    EdgeEnd() {v=nullptr;}
    EdgeEnd(Vertex<U,T>& vr) {v=&vr;}
    EdgeEnd(Vertex<U,T>& vr,bool b1,bool b2) : EdgeEnd(vr) {enter=b1;leave=b2;}

    bool enters(const Vertex<U,T>& vr) const {return (v==&vr && enter);}
    bool leaves(const Vertex<U,T>& vr) const {return (v==&vr && leave);}
    inline bool enters() const {return enters(*v);}
    inline bool exits() const {return leaves(*v);}
    inline bool leaves() const {return leaves(*v);}
    bool connects(const Vertex<U,T>& vr) const {return (v==&vr);}
    Vertex<U,T>& getVertex() {return *v;}
    const Vertex<U,T>& getVertex() const {return *v;}

     inline operator Vertex<U,T>() const {return getVertex();}
     inline operator Vertex<U,T>() {return getVertex();}
  };

  template<class T,class U> class Vertex
  {
    T data;
    std::vector<Edge<U,T>*> edges;
  public:
    bool flag=false;

  public:
    Vertex() {}
    Vertex(const T& d) {data=d;}

    void setData(const T& d) {data=d;}
    const T& getData() const {return data;}
    T& getData() {return data;}

    void addEdge(Edge<U,T>& e) {edges.push_back(&e);}
    const std::vector<Edge<U,T>*>& getAllEdges() const {return edges;}
    std::vector<Edge<U,T>*>& getAllEdges() {return edges;}
    std::vector<const Edge<U,T>*> getLeavingEdges() const ;
    std::vector<Edge<U,T>*> getLeavingEdges();
    std::vector<const Edge<U,T>*> getEnteringEdges() const ;
    std::vector<Edge<U,T>*> getEnteringEdges() ;
    int getInflux() const {return getEnteringEdges().size();}
    int getOutflux() const {return getLeavingEdges().size();}

    inline Edge<U,T>* joinToVertex(Vertex<T,U>& v,const U& data=U(),bool to=true,bool from=true);
  };

  template<class T,class U> class Edge
  {
    EdgeEnd<T,U> v1,v2;
    T data;
  public:
    bool flag=false;
  protected:
    Edge(Vertex<U,T>& vr1,Vertex<U,T>& vr2,bool direct=false) {v1=EdgeEnd<T,U>(vr1,true,!direct); v2=EdgeEnd<T,U>(vr2,!direct,true); }
    Edge() {}
    Edge(const EdgeEnd<T,U>& e1,const EdgeEnd<T,U>& e2)
    {
      if(!e1.enters()) {v1=e2; v2=e1;}
      else if(e2.enters()) {v1=e1;v2=e2;}
      else throw BadEdgeDesignException();
    }
  public:

    void setData(const T& d) {data=d;}
    const T& getData() const {return data;}
    T getData() {return data;}

    const EdgeEnd<T,U>& getEnd1() const {return v1;}
    const EdgeEnd<T,U>& getEnd2() const {return v2;}
    const EdgeEnd<T,U>& getSecondEnd(const Vertex<U,T>& fv) const
    {
      if(v1.connects(fv)) return v2;
      else if(v2.connects(fv)) return v1;
      else throw EdgeDoesNotConnectVertexException();
    }
    EdgeEnd<T,U>& getEnd1()  {return v1;}
    EdgeEnd<T,U>& getEnd2()  {return v2;}
    EdgeEnd<T,U>& getSecondEnd(const Vertex<U,T>& fv)
    {
      if(v1.connects(fv)) return v2;
      else if(v2.connects(fv)) return v1;
      else throw EdgeDoesNotConnectVertexException();
    }

    inline bool enters(const Vertex<U,T>& v) const {return v1.enters(v) || v2.enters(v);}
    inline bool exits(const Vertex<U,T>& v) const {return v2.leaves(v) || v1.leaves(v);}
    inline bool leaves(const Vertex<U,T>& v) const {return v2.leaves(v) || v1.leaves(v);}

    friend class graphs::Vertex<U,T>;
  };

  template<class V,class E> struct Graph
  {
    std::vector<Vertex<V,E>*> vertices;
    std::vector<Edge<E,V>*> edges;

  public:
    Graph(std::vector<Vertex<V,E>*> v,std::vector<Edge<E,V>*> e) {vertices=v; edges=e;}
  //public:
    Graph() {}
    Graph(Vertex<V,E>& ver)
    {
      //vertices[0]->flag=false;
      //growGraph(*vertices[0]);
      growGraph(ver);
      resetVertexFlags();
      resetEdgeFlags();
    }



    inline void deleteGraphComponents()
    {
      for(auto e : edges) delete e;
      for(auto v : vertices) delete v;
    }

    void growGraph(Vertex<V,E>& ver)
    {
      vertices.push_back(&ver);
      ver.flag=true;
      for(Edge<E,V>* edg : ver.getAllEdges())
      {
        Edge<E,V>& ed=*edg;
        if(!ed.flag)
        {
          edges.push_back(&ed); ed.flag=true;
          Vertex<V,E>* tv=&(ed.getSecondEnd(ver).getVertex());
          if(!tv->flag) growGraph(*tv);
        }
      }
    }
    inline void resetVertexFlags() {for(Vertex<V,E>* ver : vertices) ver->flag=false;}
    inline void resetEdgeFlags() {for(Edge<E,V>* ed : edges) ed->flag=false;}

    inline int indexOf(const Vertex<V,E>& v) const {return ::indexOf(&v,vertices);}

    inline const std::vector<Vertex<V,E>*>& getVertices() const {return vertices;}
    inline const std::vector<Edge<E,V>*>& getEdges() const {return edges;}

    inline const Vertex<V,E>& operator[](int i) const {return *vertices[i];}
    inline Vertex<V,E>& operator[](int i) {return *vertices[i];}

    const Vertex<V,E>* getVertexByData(const V& d) const
    {
      for(int i=0;i<vertices.size();i++)
      {
        if(vertices[i]->getData()==d) return vertices[i];
      }
      return nullptr;
    }
    Vertex<V,E>* getVertexByData(const V& d)
    {
      for(int i=0;i<vertices.size();i++)
      {
        if(vertices[i]->getData()==d) return vertices[i];
      }
      return nullptr;
    }

    virtual int getEdgeCount() const {return edges.size();}
    virtual int getVertexCount() const {return vertices.size();}

    template<class Vc=V,class Ec=E> Graph<Vc,Ec> copy(const Vc& tt=Vc(),const Ec& uu=Ec()) const //Copy while casting objects of type V to Vc and E to Ec (Useful for up/down scaling. See Residual Network)
    {
      std::vector<Vertex<Vc,Ec>*> verts;
      std::vector<Edge<Ec,Vc>*> edgs;
      for(auto v : vertices) verts.push_back(new Vertex<Vc,Ec>((Vc)(v->getData())));
      for(int i=0;i<verts.size();i++)
      {
        auto el=vertices[i]->getAllEdges();
        for(auto ed : el)
        {
          Vertex<V,E>* vp = &(ed->getSecondEnd(*vertices[i]).getVertex());
          int ind=::indexOf<Vertex<V,E>*>(vp,vertices);
          if(ind==-1) throw VertexNotInGraphException();
          if(ind>=i) edgs.push_back(verts[i]->joinToVertex(*verts[ind],(Ec)(ed->getData()),ed->enters(*vertices[ind]),ed->enters(*vertices[i])));
        }
      }

      return Graph<Vc,Ec>(verts,edgs);
    }
  };
}
//Completing missing functions
template<class T,class U> std::vector<const graphs::Edge<U,T>*> graphs::Vertex<T,U>::getEnteringEdges() const //This could be wrong syntactically. But it must return a list of edge POINTERS only
{
  std::vector<const graphs::Edge<U,T>*> ret;
  for(const graphs::Edge<U,T>* ed : edges)
  {
    if((&(ed.getEnd1()->getVertex())==this && ed->getEnd1().enters()) || (&(ed->getEnd2().getVertex())==this && ed->getEnd2().enters())) ret.push_back(ed);
  }
  return ret;
}
template<class T,class U> std::vector<graphs::Edge<U,T>*> graphs::Vertex<T,U>::getEnteringEdges() //(Copy of above - the "const" part) This could be wrong syntactically. But it must return a list of edge POINTERS only
{
  std::vector<graphs::Edge<U,T>*> ret;
  for(graphs::Edge<U,T>* ed : edges)
  {
    if((&(ed->getEnd1().getVertex())==this && ed->getEnd1().enters()) || (&(ed->getEnd2().getVertex())==this && ed->getEnd2().enters())) ret.push_back(ed);
  }
  return ret;
}
template<class T,class U> std::vector<const graphs::Edge<U,T>*> graphs::Vertex<T,U>::getLeavingEdges() const
{
  std::vector<const graphs::Edge<U,T>*> ret;
  for(const graphs::Edge<U,T>* ed : edges)
  {
    if((&(ed->getEnd1().getVertex())==this && ed->getEnd1().leaves()) || (&(ed->getEnd2().getVertex())==this && ed->getEnd2().leaves())) ret.push_back(ed);
  }
  return ret;
}
template<class T,class U> std::vector<graphs::Edge<U,T>*> graphs::Vertex<T,U>::getLeavingEdges()
{
  std::vector<graphs::Edge<U,T>*> ret;
  for(graphs::Edge<U,T>* ed : edges)
  {
    if((&(ed->getEnd1().getVertex())==this && ed->getEnd1().leaves()) || (&(ed->getEnd2().getVertex())==this && ed->getEnd2().leaves())) ret.push_back(ed);
  }
  return ret;
}
template<class T,class U> inline graphs::Edge<U,T>* graphs::Vertex<T,U>::joinToVertex(Vertex<T,U>& v,const U& d,bool to,bool from)
{
  Edge<U,T>* ne;
  if(to && from) ne=new Edge<U,T>(*this,v);
  else if(to && !from) ne=new Edge<U,T>(v,*this,true);
  else if(from && !to) ne=new Edge<U,T>(*this,v,true);
  else {throw BadEdgeDesignException();}
  v.addEdge(*ne);
  if(this!=&v) this->addEdge(*ne); //Is the check needed? Duplicates are fine (right?)
  ne->setData(d);
  return ne;
}

template<class T,class U=T> std::ostream& operator<<(std::ostream& os,const graphs::Vertex<T,U> ver) {return (os << ver.getData());}
template<class T,class U=T> std::ostream& operator<<(std::ostream& os,const graphs::EdgeEnd<T,U> ed) {return (os << ed.getVertex().getData());}
template<class T,class U=T> std::ostream& operator<<(std::ostream& os,const graphs::Edge<T,U> ed) {os << ed.getEnd1().getVertex().getData()<<" - "<<ed.getEnd2().getVertex().getData(); return os;}
template<class T,class U=T> std::ostream& operator<<(std::ostream& os,const graphs::Graph<T,U> grp)
{
  os << "Graph\n";
  for(const auto& v : grp.getVertices())
  {
    os << *v << ":\t";
    for(const auto& e : v->getAllEdges())
    {
      //cout <<"  "<<e.enters(*v)<<e.leaves(*v)<<"-";
      os << e->getSecondEnd(*v) << (e->leaves(*v)?((e->enters(*v))?"":"(to)"):"(fro)")<<"("<<e->getData()<<"), ";
    }
    os << "\n";
  }
  return os;
}

namespace graphs
{
  namespace algo
  {
    struct FlowEdge
    {
      double flow,cap;
      FlowEdge() {flow=0; cap=0;}
      explicit FlowEdge(double v) {flow=v; cap=v;}
      FlowEdge(double f,double c) {flow=f; cap=c;}
      template<class T> FlowEdge(const std::pair<T,T>& p) {flow=get<0>(p); cap=get<1>(p); assert(flow<=cap);}

      inline bool check() const {return flow<=cap;}

      operator double() const {return flow;}
    };
    static inline std::ostream& operator<<(std::ostream& os,const FlowEdge& fl) {return (os << fl.flow <<"/"<<fl.cap);}

    template<class T,class U=T> std::vector<Edge<U,T>> eulerPathInternal(Graph<T,U>& G,Vertex<T,U>& start, std::vector<Edge<U,T>> path=std::vector<Edge<U,T>>())
    {
      if(path.size()>G.getEdges().size()) throw NoEulerPathPossibleException();
      for(Edge<U,T>* e : start.getAllEdges())
      {
        if(e->flag) continue;
        e->flag=true;
        path.push_back(*e);
        Vertex<T,U>& nv=e->getSecondEnd(start).getVertex();
        return eulerPathInternal(G,nv,path);
      }
      return path;
    }
    //Wrong algorithm
    template<class T,class U=T> std::vector<Edge<U,T>> eulerPath(Graph<T,U>& G,Vertex<T,U>& start, std::vector<Edge<U,T>> path=std::vector<Edge<U,T>>()) {G.resetEdgeFlags(); return eulerPathInternal(G,start,path);}

    //template<class T> Graph<T,double> Graph<T,FlowEdge>::copy<T,double>();
    template<class T> static Graph<T,double> getResidualNetwork(Graph<T,FlowEdge> gr)
    {
      Graph<T,FlowEdge> ret=gr.copy();
      std::vector<int> edgecounts;
      for(auto& v : ret.getVertices()) edgecounts.push_back(v->getLeavingEdges().size());
      //cout << "B has "<<ret.getVertexByData("B")->getLeavingEdges().size()<<" outflowing edges\n";
      int K=0;
      for(auto& v : ret.getVertices())
      {
        auto edgs=v->getLeavingEdges();
        for(int j=0;j<edgecounts[K];j++)
        {
          auto& e=edgs[j];
          //cout <<"\t"<< e->getSecondEnd(*v).getVertex().getData()<<" is sub of "<<v->getData()<<"\n";
          if(e->getData().flow<=0) continue;
          e->getSecondEnd(*v).getVertex().joinToVertex(*v,FlowEdge(e->getData().flow),true,false);
        }
        K++;
      }
      for(auto e : ret.getEdges()) e->setData(FlowEdge(e->getData().cap-e->getData().flow));
      Graph<T,double> ret2=ret.copy(T(),(double)0.0d);
      ret.deleteGraphComponents();
      return ret2;
    }

    template<class T,class U> static bool edgeIsTraversible(const Edge<T,U>& e) {return e.getData()!=0.0d;}
    template<class T,class U> static bool alwaysTrueTraversible(const Edge<T,U>& e) {return true;}
    template<class T> static bool edgeIsTraversibleFlowNetwork(const Edge<FlowEdge,T>& e) {return e.getData().flow<e.getData().cap;}


    template<class T,class U=double> std::pair<std::vector<Vertex<T,U>*>,std::vector<Edge<U,T>*>> findPathInGraphInternal(Graph<T,U>& gr,Vertex<T,U>* s,Vertex<T,U>* t,bool (*trav)(const Edge<U,T>&)=edgeIsTraversible,int v=0,std::vector<Vertex<T,U>*> path=std::vector<Vertex<T,U>*>(),std::vector<Edge<U,T>*> pathedges=std::vector<Edge<U,T>*>())
    {
      if(v==gr.getVertices().size()) throw NoPathPossibleException();
      v+=1;
      s->flag=true;
      path.push_back(s);
      auto eds=s->getLeavingEdges();
      for(auto e : eds)
      {
        if(!trav(*e)) continue;
        auto& ver=e->getSecondEnd(*s).getVertex();
        if(ver.flag) continue;
        pathedges.push_back(e);
        if(&ver==t) {path.push_back(t); return make_pair(path,pathedges);}
        try {return findPathInGraphInternal(gr,&ver,t,trav,v,path,pathedges);}
        catch(NoPathPossibleException ex) {continue;}
      }
      throw NoPathPossibleException();
    }
    template<class T,class U=double> std::pair<std::vector<Vertex<T,U>*>,std::vector<Edge<U,T>*>> findPathInGraph(Graph<T,U>& gr,Vertex<T,U>& s,Vertex<T,U>& t,bool (*trav)(const Edge<U,T>&)=edgeIsTraversible)
    {
      gr.resetVertexFlags();
      return findPathInGraphInternal(gr,&s,&t,trav,0);
    }
    template<class T> double checkAndComputeFlow(Graph<T,FlowEdge>& gr,Vertex<T,FlowEdge>& s,Vertex<T,FlowEdge>& t,double err=1e-3)
    {
      if(s.getEnteringEdges().size() || t.getLeavingEdges().size()) throw InvalidFlowException("Source/Target must have only outgoing/incoming edges respectively");
      double ret=0;
      for(auto v : gr.getVertices())
      {
        double infl=0,outfl=0;
        for(auto e : v->getLeavingEdges())
        {
          if(!e->getData().check()) throw InvalidFlowException("Flow in some edge is more than the capacity!");
          outfl+=e->getData().flow;
        }
        for(auto e : v->getEnteringEdges())
        {
          if(!e->getData().check()) throw InvalidFlowException("Flow in some edge is more than the capacity!");
          infl+=e->getData().flow;
        }
        if((v!=&s && v!=&t) && ::abs(infl-outfl)>err) throw InvalidFlowException("Net inflow and outflow at some vertex don't match");
        if(v==&s) ret=outfl;
        if(v==&t) ret=infl;
      }
      return ret;
    }

    template<class V,class E> class DFSTree : public Graph<V,E>
    {
      int trees=0;
    public:
      template<class T=V,class U=E> DFSTree(Graph<T,U>& g,Vertex<T,U>& start,bool autostop=false,bool reverse=false,bool (*trav)(const Edge<U,T>&)=alwaysTrueTraversible<U,T>) : Graph<V,E>()
      {
        g.resetVertexFlags();
        this->growGraph(*generateTreeInternal(g,start,trav,autostop,reverse));
      }
    private:
      template<class T,class U> Vertex<V,E>* generateTreeInternal(Graph<T,U>& g,Vertex<T,U>& start,bool (*trav)(const Edge<U,T>&),bool autostop=false,bool reverse=false,Vertex<V,E>* parentlocal=nullptr,const E& dt=E(),bool e1=true,bool e2=true,int vc=0)
      {
        start.flag=true;
        Vertex<V,E>* nver=new Vertex<V,E>((V)start.getData());
        vc++;
        if(parentlocal) parentlocal->joinToVertex(*nver,dt,e1,e2);
        auto edges=(reverse)?start.getEnteringEdges():start.getLeavingEdges();
        bool proc=false;
        for(auto& e : edges)
        {
          auto& cver=e->getSecondEnd(start).getVertex();
          if(cver.flag || !trav(*e)) continue;
          generateTreeInternal(g,cver,trav,autostop,reverse,nver,(E)e->getData(),e->enters(cver),e->enters(start));
        }
        if(vc>=g.getVertexCount() || autostop) return nver;
        for(auto& ver : g.getVertices())
        {
          if(!ver->flag) {generateTreeInternal(g,*ver,trav,autostop,reverse); this->growGraph(*ver);}
        }
        return nver;
      }
    public:
      //~DFSTree() {this->deleteGraphComponents();}
    };

    template<class T> std::pair<double,DFSTree<T,double>> findMaxFlowAndMinCut(Graph<T,FlowEdge>& gr,Vertex<T,FlowEdge>& s,Vertex<T,FlowEdge>& t,double err=1e-3,bool verbose=false)
    {
      double flow=checkAndComputeFlow(gr,s,t,err);
      double oldflow=flow-2*err;
      int sind,tind;
      sind=gr.indexOf(s); tind=gr.indexOf(t);
      while(::abs(flow-oldflow)>err)
      {
        oldflow=flow;
        Graph<T,double> resnet=getResidualNetwork(gr);

        if(verbose)
        {
          cout << "Network:\n\n";
          for(const auto& vert : gr.getVertices())
          {
            cout << vert->getData() << ":\t";
            for(const auto& e : vert->getAllEdges())
            {
              //cout <<"  "<<e.enters(*v)<<e.leaves(*v)<<"-";
              cout << e->getSecondEnd(*vert).getVertex().getData() << (e->leaves(*vert)?((e->enters(*vert))?"":"(to)"):"(fro)")<<"("<<e->getData()<<"), ";
            }
            cout << "\n";
          }
          cout <<"\n\nResidual Network:\n\n";
          //graphs::algo::operator<<(std::cout,resnet);
          for(const auto& vert : resnet.getVertices())
          {
            cout << vert->getData() << ":\t";
            for(const auto& e : vert->getAllEdges())
            {
              //cout <<"  "<<e.enters(*v)<<e.leaves(*v)<<"-";
              cout << e->getSecondEnd(*vert).getVertex().getData() << (e->leaves(*vert)?((e->enters(*vert))?"":"(to)"):"(fro)")<<"("<<e->getData()<<"), ";
            }
            cout << "\n";
          }
          cout <<"--------------\n";
        }

        try
        {
          auto pth=findPathInGraph<T,double>(resnet,resnet[sind],resnet[tind]);
          std::vector<Vertex<T,double>*> pverts=get<0>(pth);
          std::vector<Edge<double,T>*> peds=get<1>(pth);
          double delta=1.0/0.0;
          for(auto ed : peds) {if(ed->getData()<delta) delta=ed->getData();}
          if(verbose) cout << delta << " is delta\n";
          std::vector<int> trace;
          for(auto v : pverts)  trace.push_back(resnet.indexOf(*v));
          for(int i=0;i<trace.size()-1;i++)
          {
            Vertex<T,FlowEdge>& vtemp=gr[trace[i]];
            for(auto& e : vtemp.getAllEdges())
            {
              if(gr.indexOf(e->getSecondEnd(vtemp).getVertex())==trace[i+1])
              {
                if(e->leaves(vtemp)) e->setData(make_pair(e->getData().flow+delta,e->getData().cap));
                else e->setData(make_pair(e->getData().flow-delta,e->getData().cap));
                break;
              }
            }
          }
          flow=checkAndComputeFlow(gr,s,t,err);
          cout << "New flow: "<<flow<<"\n";
        }
        catch(NoPathPossibleException ex) {break;}
      }
      auto resnet=getResidualNetwork(gr);
      /*cout << "Final Resnet:\n\n";
      for(const auto& vert : resnet.getVertices())
      {
        cout << vert->getData() << ":\t";
        for(const auto& e : vert->getAllEdges())
        {
          //cout <<"  "<<e.enters(*v)<<e.leaves(*v)<<"-";
          cout << e->getSecondEnd(*vert).getVertex().getData() << (e->leaves(*vert)?((e->enters(*vert))?"":"(to)"):"(fro)")<<"("<<e->getData()<<"), ";
        }
        cout << "\n";
      }
      cout <<"--------------\n";*/
      return make_pair(flow,DFSTree<T,double>(resnet,resnet[sind],true,false,edgeIsTraversible));
    }

    struct WeightedEdge
    {
      double wt;
      WeightedEdge() {wt=0;}
      explicit WeightedEdge(double v) {wt=v;}
      template<class T> WeightedEdge(const T& p) {wt=p;}


      operator double() const {return wt;}
    };

    template<class V,class E> std::vector<Vertex<V>> getCycleInternal(Graph<V,E>& g,Vertex<V>& start,std::vector<Vertex<V>> visited,const Vertex<V>& endat)
    {
      visited.push_back(start);
      for(const Edge<E,V>& ed : start.getLeavingEdges())
      {
        auto tempv=ed.getSecondEnd(start);
        if(&tempv==&endat)
        {
          if(visited.size()>=3) return visited;
          else continue;
        }
        if(contains(visited,tempv)) continue;
        try {return getCycleInternal(g,tempv,visited,start);}
        catch(NoCyclePossibleException& ex) {continue;}
      }
      throw NoCyclePossibleException();
    }
    template<class V,class E> inline std::vector<Vertex<V>> getCycle(Graph<V,E>& g,Vertex<V>& start) {return getCycleInternal(g,start,std::vector<Vertex<V>>(),start);}
    template<class T> double removeCyclicFlows(Graph<T,WeightedEdge>& gr,double err=1e-3,bool verbose=false)
    {

    }
  }
}
#endif
