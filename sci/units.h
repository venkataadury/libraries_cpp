#pragma once
#include <iostream>
#include <vector>
//#include "commons/commons.h"

#ifndef HALF_PRECISION_UNITS
typedef double numeric;
static const double ZERO=0.0,ONE=1.0,EPS=1e-10;
#else
typedef float numeric;
static const float ZERO=0.0f,ONE=1.0f,EPS=1e-10f;
#endif

namespace unitsystems
{
    class UnitIsNotPrimitiveException : public std::exception {};
    class PrimitiveUnitTypeMismatchException : public std::exception {};
    class UnitTypeMismatchException : public std::exception {};
    class NoPrimitiveConversionException : public std::exception {};

    class PrimitiveUnit;
    class Unit;
    static int primitiveID=0;

    struct StrippedUnit
    {
        std::vector<int> typeIDs;
        std::vector<int> powers;
        std::vector<numeric> values;
    };

    StrippedUnit simplifyUnit(const Unit& u,bool use_conv=false);

    class Unit
    {
    public:
        numeric extra=ZERO;
    protected:
        //numeric v;
        std::vector<PrimitiveUnit*> primitives;
        std::vector<int> powers;
        bool simplified=false;
        Unit() {}

    public:
        Unit(const std::vector<const PrimitiveUnit*>& prim,const std::vector<int>& powr); //Implemented
        Unit(const Unit& u);
        ~Unit(); //Implemented

        inline virtual numeric getVal() const; //Implemented
        inline virtual numeric toSI() const; //Implemented
        virtual bool hasSameUnit(const Unit& u2)const; //Implemented
        std::string unitString()const; //Implemented

        virtual inline void setVal(numeric val)
        {
            numeric cval=getVal();
            extra=val-cval;
        }

        inline int getPrimitiveIndex(const PrimitiveUnit* p)const; //Implemented
        inline int getPrimitiveIndex(int pid)const; //Implemented
        inline bool isPrimitive()const {return (primitives.size()==1 && powers[0]==1);}
        inline PrimitiveUnit makePrimitive()const; //Implemented
        Unit& simplify(bool use_conv=false); //Implemented
        Unit& convertTo(const PrimitiveUnit& pu);

        inline const std::vector<PrimitiveUnit*>& getPrimitives()const {return primitives;}
        inline const std::vector<int>& getPowers()const {return powers;}

        inline void operator*=(numeric num);
        inline void operator/=(numeric num);
        inline Unit operator*(numeric num)const;
        inline Unit operator/(numeric num)const;

        Unit& operator=(const Unit& u);

        friend class PrimitiveUnit;
        friend StrippedUnit simplifyUnit(const Unit& u,bool use_conv);
    };

    class PrimitiveUnit
    {
    public:
        PrimitiveUnit() = delete;
    protected:
        numeric v;
        const int typeID;
        PrimitiveUnit(int tID) : typeID(tID) {}
    public:
        PrimitiveUnit(int tID,numeric val) : PrimitiveUnit(tID) {v=val;}

        inline int getTypeID() const {return typeID;}
        virtual inline bool hasSameUnit(const PrimitiveUnit& u2)const {return typeID==u2.typeID;}
        inline numeric getVal() const {return v;}
        inline void setVal(numeric val) {v=val;}

        virtual inline void setValFromSI(numeric val,int pow=1) {v=val;}
        virtual inline numeric toSI(int pow=1) const {return v;}
        inline operator Unit() {return Unit({this},{1});}

        virtual inline PrimitiveUnit* clone()const {return new PrimitiveUnit(typeID,v);}

        inline PrimitiveUnit operator+(numeric num)const {return PrimitiveUnit(typeID,v+num);} //Implicitly assume same type
        inline PrimitiveUnit operator-(numeric num)const {return PrimitiveUnit(typeID,v-num);}
        inline PrimitiveUnit operator*(numeric num)const {return PrimitiveUnit(typeID,v*num);}
        inline PrimitiveUnit operator/(numeric num)const {return PrimitiveUnit(typeID,v/num);}
        inline void operator+=(numeric num) {v+=num;}
        inline void operator-=(numeric num) {v-=num;}
        inline void operator*=(numeric num) {v*=num;}
        inline void operator/=(numeric num) {v/=num;}

        inline PrimitiveUnit operator+(const PrimitiveUnit& p2)const;
        inline PrimitiveUnit operator-(const PrimitiveUnit& p2)const;

        inline Unit operator*(const PrimitiveUnit& p2)const {return Unit({this,&p2},{1,1});}
        inline Unit operator/(const PrimitiveUnit& p2)const
        {
            PrimitiveUnit* temp=p2.clone();
            temp->v=1/temp->v;
            Unit ret({this,temp},{1,-1});
            delete temp;
            return ret;
        }
        inline Unit operator*(const Unit& p2)const
        {
            int pidx=p2.getPrimitiveIndex(this);
            Unit ret=p2;
            if(pidx==-1) {ret.primitives.push_back(this->clone()); ret.powers.push_back(1);}
            else {(*ret.primitives[pidx])*=this->getVal(); ret.powers[pidx]++;}
            return ret;
        }
        inline Unit operator/(const Unit& p2)const
        {
            int pidx=p2.getPrimitiveIndex(this);
            Unit ret=p2;
            for(int i=0;i<ret.powers.size();i++) ret.powers[i]*=-1;
            if(pidx==-1) {ret.primitives.push_back(this->clone()); ret.powers.push_back(1);}
            else
            {
                (*ret.primitives[pidx])*=this->getVal();
                ret.powers[pidx]++;
            }
            return ret;
        }

        virtual std::string unitString()const {return std::to_string(typeID);}
    };
    class PrimitiveConversion
    {
    public:
        int fro,to;
        numeric a,b,c=ONE; //conversion is y=(ax/c+b)
        PrimitiveConversion(int f,int t) {fro=f; to=t;}

        inline numeric convert(numeric v) const {return (a*v)/c+b;}
        inline numeric invert(numeric v)const {return (c*(v-b))/a;}
    };
    static std::vector<const PrimitiveConversion*> conversions=std::vector<const PrimitiveConversion*>();
    static void addPrimitiveConversion(const PrimitiveConversion& pc) {conversions.push_back(&pc);}
    static void addPrimitiveConversion(int fro,int to,numeric a,numeric b,numeric c=ONE)
    {
        PrimitiveConversion* ta=new PrimitiveConversion(fro,to);
        ta->a=a; ta->b=b; ta->c=c;
        addPrimitiveConversion(*ta);
    }

    numeric Unit::getVal()const
    {
        if(!primitives.size()) return ZERO;
        numeric ret=1;
        for(int i=0;i<primitives.size();i++) ret*=primitives[i]->getVal();
        return ret+extra;
    }
    numeric Unit::toSI()const
    {
        if(!primitives.size()) return ZERO;
        numeric ret=1;
        for(int i=0;i<primitives.size();i++) ret*=primitives[i]->toSI(powers[i]);
        return ret+extra;
    }
    inline PrimitiveUnit Unit::makePrimitive()const
    {
        if(!this->isPrimitive()) throw UnitIsNotPrimitiveException();
        return PrimitiveUnit(primitives[0]->getTypeID(),primitives[0]->getVal());
    }

    Unit& Unit::simplify(bool use_conv)
    {
        StrippedUnit info=simplifyUnit(*this,use_conv);
        std::vector<PrimitiveUnit*> oldprimitives=primitives; primitives=std::vector<PrimitiveUnit*>(); powers=std::vector<int>();
        for(int idx=0;idx<info.typeIDs.size();idx++)
        {
            PrimitiveUnit* tocopy=nullptr;
            for(PrimitiveUnit* pu : oldprimitives)
            {
                if(pu->getTypeID()==info.typeIDs[idx])
                {
                    pu->setVal(info.values[idx]);
                    primitives.push_back(pu);
                    break;
                }
            }
            powers.push_back(info.powers[idx]);
        }
        return (*this);
    }
    Unit& Unit::convertTo(const PrimitiveUnit& pu)
    {
        std::vector<PrimitiveUnit*> oldprimitives=primitives;
        std::vector<int> oldpowers=powers;
        primitives=std::vector<PrimitiveUnit*>();
        powers=std::vector<int>();
        PrimitiveUnit* temp=pu.clone(); temp->setVal(ONE);
        primitives.push_back(temp);
        powers.push_back(0);
        for(int i=0;i<oldprimitives.size();i++)
        {
            primitives.push_back(oldprimitives[i]);
            powers.push_back(oldpowers[i]);
        }
        return this->simplify(true);
    }

    std::string Unit::unitString()const
    {
        std::string ret="";
        for(int i=0;i<primitives.size()-1;i++)
        {
            if(powers[i])
            {
                ret+=primitives[i]->unitString();
                if(powers[i]!=1) ret+=std::to_string(powers[i]);
                if(powers[i+1]) ret+=".";
            }
        }
        if(powers[primitives.size()-1])
        {
            ret+=primitives[primitives.size()-1]->unitString();
            if(powers[powers.size()-1]!=1) ret+=std::to_string(powers[powers.size()-1]);
        }
        return ret;
    }

    bool Unit::hasSameUnit(const Unit& u2)const
    {
        StrippedUnit u1s=simplifyUnit(*this),u2s=simplifyUnit(u2);
        for(int i=0;i<u1s.typeIDs.size();i++)
        {
            int j=0;
            for(;j<u2s.typeIDs.size();j++)
            {
                if(u1s.typeIDs[i]==u2s.typeIDs[j])
                {
                    if(u1s.powers[i]!=u2s.powers[j]) return false;
                }
            }
            if(j>=u2s.typeIDs.size() && u1s.powers[i]!=0) return false;
        }
        return true;
    }

    inline int Unit::getPrimitiveIndex(const PrimitiveUnit* p)const {return this->getPrimitiveIndex(p->getTypeID());}
    inline int Unit::getPrimitiveIndex(int pid)const
    {
        for(int i=0;i<primitives.size();i++)
        {
            PrimitiveUnit* pu=primitives[i];
            if(pu->getTypeID()==pid) return i;
        }
        return -1;
    }

    namespace unitconversions
    {
        static numeric convertPrimitive(numeric val,int f,int t,int pow=1,const std::vector<const PrimitiveConversion*>& skipconvs=std::vector<const PrimitiveConversion*>())
        {
            std::vector<const PrimitiveConversion*> usedconvs=skipconvs;
            for(const PrimitiveConversion* conv : conversions)
            {
                bool found=false;
                for(const PrimitiveConversion* uc : skipconvs)
                {
                    if(uc==conv) {found=true; break;}
                }
                if(found) continue;

                if((conv->fro==conv->to) || (conv->fro!=f && conv->to!=f))continue;
                if(conv->fro==f)
                {
                    numeric nxt=val; if(pow>0) for(int j=0;j<pow;j++) nxt=conv->convert(nxt); else for(int j=0;j>pow;j--) nxt=conv->invert(nxt);
                    if(conv->to==t) return nxt;
                    try
                    {
                        usedconvs.push_back(conv);
                        return convertPrimitive(nxt,conv->to,t,pow,usedconvs);
                    }
                    catch(const NoPrimitiveConversionException& ex) {}
                }
                if(conv->to==f)
                {
                    numeric nxt=val; if(pow>0) for(int j=0;j<pow;j++) nxt=conv->invert(nxt); else for(int j=0;j>pow;j--) nxt=conv->convert(nxt);
                    if(conv->fro==t) return nxt;
                    try
                    {
                        usedconvs.push_back(conv);
                        return convertPrimitive(nxt,conv->fro,t,pow,usedconvs);
                    }
                    catch(const NoPrimitiveConversionException& ex) {}
                }
            }
            throw NoPrimitiveConversionException();
        }
    }
}

unitsystems::Unit::Unit(const std::vector<const PrimitiveUnit*>& prim,const std::vector<int>& powr)
{
    for(const PrimitiveUnit* pu : prim) primitives.push_back(pu->clone());
    powers=powr;
}
unitsystems::Unit::Unit(const unitsystems::Unit& u)
{
    extra=u.extra;
    powers=u.powers;
    simplified=u.simplified;
    for(unitsystems::PrimitiveUnit* pu : u.primitives) primitives.push_back(pu->clone());
}
unitsystems::Unit::~Unit() {for(PrimitiveUnit* pu : primitives) delete pu;}
unitsystems::Unit& unitsystems::Unit::operator=(const unitsystems::Unit& u)
{
    extra=u.extra;
    powers=u.powers;
    simplified=u.simplified;
    for(unitsystems::PrimitiveUnit* pu : primitives) delete pu;
    primitives=std::vector<unitsystems::PrimitiveUnit*>();
    for(unitsystems::PrimitiveUnit* pu : u.primitives) primitives.push_back(pu->clone());
    return *this;
}

inline unitsystems::PrimitiveUnit unitsystems::PrimitiveUnit::operator+(const unitsystems::PrimitiveUnit& p2)const
{
    if(typeID!=p2.typeID) throw PrimitiveUnitTypeMismatchException();
    return PrimitiveUnit(typeID,getVal()+p2.getVal());
}
inline unitsystems::PrimitiveUnit unitsystems::PrimitiveUnit::operator-(const unitsystems::PrimitiveUnit& p2)const
{
    if(typeID!=p2.typeID) throw PrimitiveUnitTypeMismatchException();
    return PrimitiveUnit(typeID,getVal()-p2.getVal());
}


unitsystems::Unit operator+(const unitsystems::Unit& u1,const unitsystems::Unit& u2)
{
    if(!u1.hasSameUnit(u2)) throw unitsystems::UnitTypeMismatchException();
    unitsystems::Unit ret=u1;
    ret.extra+=u2.getVal();
    return ret;
}
unitsystems::Unit operator-(const unitsystems::Unit& u1,const unitsystems::Unit& u2)
{
    if(!u1.hasSameUnit(u2)) throw unitsystems::UnitTypeMismatchException();
    unitsystems::Unit ret=u1;
    ret.extra-=u2.getVal();
    return ret;
}
unitsystems::Unit operator*(const unitsystems::Unit& u1,const unitsystems::Unit& u2)
{
    std::vector<const unitsystems::PrimitiveUnit*> pul; std::vector<int> pup;
    for(const auto& pu : u1.getPrimitives()) pul.push_back(pu);
    for(const auto& pu : u2.getPrimitives()) pul.push_back(pu);
    for(const auto& p : u1.getPowers()) pup.push_back(p);
    for(const auto& p : u2.getPowers()) pup.push_back(p);
    return unitsystems::Unit(pul,pup);
}
unitsystems::Unit operator/(const unitsystems::Unit& u1,const unitsystems::Unit& u2)
{
    std::vector<const unitsystems::PrimitiveUnit*> pul; std::vector<int> pup;
    for(const auto& pu : u1.getPrimitives()) pul.push_back(pu);
    std::vector<unitsystems::PrimitiveUnit*> todel;
    for(const auto& pu : u2.getPrimitives())
    {
        unitsystems::PrimitiveUnit* npu=pu->clone(); npu->setVal(ONE/pu->getVal()); //new unitsystems::PrimitiveUnit(pu->getTypeID(),ONE/pu->getVal());
        pul.push_back(npu);
        todel.push_back(npu);
    }
    for(const auto& p : u1.getPowers()) pup.push_back(p);
    for(const auto& p : u2.getPowers()) pup.push_back(-p);
    unitsystems::Unit ret(pul,pup);
    for(unitsystems::PrimitiveUnit* dpu : todel) delete dpu;
    return ret;
}

unitsystems::Unit unitsystems::Unit::operator*(numeric num)const
{
    unitsystems::Unit ret(*this);
    ret.extra*=num;
    if(primitives.size()) (*(ret.primitives[0]))*=num;
    return ret;
}
unitsystems::Unit unitsystems::Unit::operator/(numeric num)const
{
    unitsystems::Unit ret(*this);
    ret.extra/=num;
    if(primitives.size()) (*(ret.primitives[0]))/=num;
    return ret;
}
static inline unitsystems::Unit operator*(numeric num,const unitsystems::Unit& u) {return u*num;}

void unitsystems::Unit::operator*=(numeric num)
{
    extra*=num;
    if(primitives.size()) (*primitives[0])*=num;
}
void unitsystems::Unit::operator/=(numeric num)
{
    extra/=num;
    if(primitives.size()) (*primitives[0])/=num;
}

std::ostream& operator<<(std::ostream& os,const unitsystems::PrimitiveUnit& pu)
{
    os << pu.getVal()<<" "<<pu.unitString();
    return os;
}
std::ostream& operator<<(std::ostream& os,const unitsystems::Unit& u)
{
    os << u.getVal()<<" "<<u.unitString();
    return os;
}

namespace unitsystems
{
    namespace SIUnits
    {
        class Length : public PrimitiveUnit //ID=1
        {
        public:
            Length(numeric val) : PrimitiveUnit(1,val) {}

            std::string unitString()const override {return "m";}
            virtual inline Length* clone()const override{return new Length(v);}

            inline Length operator+(numeric num)const {return Length(v+num);} //Implicitly assume same type
            inline Length operator-(numeric num)const {return Length(v-num);}
            inline Length operator*(numeric num)const {return Length(v*num);}
            inline Length operator/(numeric num)const {return Length(v/num);}
        } METRE(ONE);
        class Mass : public PrimitiveUnit //ID=2
        {
        public:
            Mass(numeric val) : PrimitiveUnit(2,val) {}

            std::string unitString()const override {return "kg";}
            virtual inline Mass* clone()const override{return new Mass(v);}

            inline Mass operator+(numeric num)const {return Mass(v+num);} //Implicitly assume same type
            inline Mass operator-(numeric num)const {return Mass(v-num);}
            inline Mass operator*(numeric num)const {return Mass(v*num);}
            inline Mass operator/(numeric num)const {return Mass(v/num);}
        } KILOGRAM(ONE);
        class Time : public PrimitiveUnit //ID=3
        {
        public:
            Time(numeric val) : PrimitiveUnit(3,val) {}

            std::string unitString()const override {return "s";}
            virtual inline Time* clone()const override{return new Time(v);}

            inline Time operator+(numeric num)const {return Time(v+num);} //Implicitly assume same type
            inline Time operator-(numeric num)const {return Time(v-num);}
            inline Time operator*(numeric num)const {return Time(v*num);}
            inline Time operator/(numeric num)const {return Time(v/num);}
        } SECOND(ONE);
        class Temperature : public PrimitiveUnit //ID=4
        {
        public:
            Temperature(numeric val) : PrimitiveUnit(4,val) {}

            std::string unitString()const override {return "K";}
            virtual inline Temperature* clone()const override{return new Temperature(v);}

            inline Temperature operator+(numeric num)const {return Temperature(v+num);} //Implicitly assume same type
            inline Temperature operator-(numeric num)const {return Temperature(v-num);}
            inline Temperature operator*(numeric num)const {return Temperature(v*num);}
            inline Temperature operator/(numeric num)const {return Temperature(v/num);}
        } KELVIN(ONE);
        class Charge : public PrimitiveUnit //ID=5
        {
        public:
            Charge(numeric val) : PrimitiveUnit(5,val) {}

            std::string unitString()const override {return "C";}
            virtual inline Charge* clone()const override{return new Charge(v);}

            inline Charge operator+(numeric num)const {return Charge(v+num);} //Implicitly assume same type
            inline Charge operator-(numeric num)const {return Charge(v-num);}
            inline Charge operator*(numeric num)const {return Charge(v*num);}
            inline Charge operator/(numeric num)const {return Charge(v/num);}
        } COULOMB(ONE);
        class Moles : public PrimitiveUnit //ID=6
        {
        public:
            Moles(numeric val) : PrimitiveUnit(6,val) {}

            std::string unitString()const override {return "mol";}
            virtual inline Moles* clone()const override{return new Moles(v);}

            inline Moles operator+(numeric num)const {return Moles(v+num);} //Implicitly assume same type
            inline Moles operator-(numeric num)const {return Moles(v-num);}
            inline Moles operator*(numeric num)const {return Moles(v*num);}
            inline Moles operator/(numeric num)const {return Moles(v/num);}
        } MOLES(ONE);
    }

    namespace atomic_units
    {
        class Length : public PrimitiveUnit //ID=7
        {
        public:
            Length(numeric val) : PrimitiveUnit(7,val) {}

            std::string unitString()const override {return "nm";} //Nanometers
            virtual inline Length* clone()const override{return new Length(v);}

            virtual inline void setValFromSI(numeric val,int powr=1) {v=val*::pow(1e-9,-powr);}
            virtual inline numeric toSI(int powr=1) const {return v*::pow(1e-9,powr);}

            inline Length operator+(numeric num)const {return Length(v+num);} //Implicitly assume same type
            inline Length operator-(numeric num)const {return Length(v-num);}
            inline Length operator*(numeric num)const {return Length(v*num);}
            inline Length operator/(numeric num)const {return Length(v/num);}
        } NANOMETRE(ONE);
        class Mass : public PrimitiveUnit //ID=8
        {
        public:
            Mass(numeric val) : PrimitiveUnit(8,val) {}

            std::string unitString()const override {return "amu";}
            virtual inline Mass* clone()const override{return new Mass(v);}

            virtual inline void setValFromSI(numeric val,int powr=1) {v=val*::pow(1.6603145e-27,-powr);}
            virtual inline numeric toSI(int powr=1) const {return v*::pow(1.6603145e-27,powr);}

            inline Mass operator+(numeric num)const {return Mass(v+num);} //Implicitly assume same type
            inline Mass operator-(numeric num)const {return Mass(v-num);}
            inline Mass operator*(numeric num)const {return Mass(v*num);}
            inline Mass operator/(numeric num)const {return Mass(v/num);}
        } ATOMICMASS(ONE);
        class Time : public PrimitiveUnit //ID=9
        {
        public:
            Time(numeric val) : PrimitiveUnit(9,val) {}

            std::string unitString()const override {return "ps";} //Picoseconds
            virtual inline Time* clone()const override{return new Time(v);}

            virtual inline void setValFromSI(numeric val,int powr=1) {v=val*::pow(1e-12,-powr);}
            virtual inline numeric toSI(int powr=1) const {return v*::pow(1e-12,powr);}

            inline Time operator+(numeric num)const {return Time(v+num);} //Implicitly assume same type
            inline Time operator-(numeric num)const {return Time(v-num);}
            inline Time operator*(numeric num)const {return Time(v*num);}
            inline Time operator/(numeric num)const {return Time(v/num);}
        } PICOSECOND(ONE);
        class Temperature : public PrimitiveUnit //ID=4
        {
        public:
            Temperature(numeric val) : PrimitiveUnit(4,val) {}

            std::string unitString()const override {return "K";}
            virtual inline Temperature* clone()const override{return new Temperature(v);}

            inline Temperature operator+(numeric num)const {return Temperature(v+num);} //Implicitly assume same type
            inline Temperature operator-(numeric num)const {return Temperature(v-num);}
            inline Temperature operator*(numeric num)const {return Temperature(v*num);}
            inline Temperature operator/(numeric num)const {return Temperature(v/num);}
        } KELVIN(ONE);
        class Charge : public PrimitiveUnit //ID=11
        {
        public:
            Charge(numeric val) : PrimitiveUnit(11,val) {}

            std::string unitString()const override {return "e";} //Charge on an electron
            virtual inline Charge* clone()const override{return new Charge(v);}

            virtual inline void setValFromSI(numeric val,int powr=1) {v=val*::pow(1.602e-19,-powr);}
            virtual inline numeric toSI(int powr=1) const {return v*::pow(1.602e-19,powr);}

            inline Charge operator+(numeric num)const {return Charge(v+num);} //Implicitly assume same type
            inline Charge operator-(numeric num)const {return Charge(v-num);}
            inline Charge operator*(numeric num)const {return Charge(v*num);}
            inline Charge operator/(numeric num)const {return Charge(v/num);}
        } ELECTRONCHARGE(ONE);
        class Moles : public PrimitiveUnit //ID=12
        {
        public:
            Moles(numeric val) : PrimitiveUnit(12,val) {}

            std::string unitString()const override {return "at";}
            virtual inline Moles* clone()const override{return new Moles(v);}

            virtual inline void setValFromSI(numeric val,int powr=1) {v=val*::pow(1.6603e-24,-powr);}
            virtual inline numeric toSI(int powr=1) const {return v*::pow(1.6603e-24,powr);}

            inline Moles operator+(numeric num)const {return Moles(v+num);} //Implicitly assume same type
            inline Moles operator-(numeric num)const {return Moles(v-num);}
            inline Moles operator*(numeric num)const {return Moles(v*num);}
            inline Moles operator/(numeric num)const {return Moles(v/num);}

        } ATOM(ONE);
    }

    /*namespace SYSTEM_NAME
    {

    }*/ //Add a unit-system here if needed

    static void setupStandardConversions()
    {
        addPrimitiveConversion(7,1,1e-9,ZERO);
        addPrimitiveConversion(8,2,1.6603145e-27,ZERO);
        addPrimitiveConversion(9,3,1e-12,ZERO);
        addPrimitiveConversion(11,5,1.602e-19,ZERO);
        addPrimitiveConversion(12,6,1.6603e-24,ZERO);

        //Add more conversions
    }

    StrippedUnit simplifyUnit(const Unit& u,bool use_conv)
    {
        std::vector<int> typeIDs1,powers1;
        std::vector<numeric> vals1;
        for(int idx=0;idx<u.primitives.size();idx++)
        {
            PrimitiveUnit* pu = u.primitives[idx];
            int fid=-1;
            numeric cV=ZERO;
            for(int i=0;i<typeIDs1.size();i++)
            {
                if(typeIDs1[i]==pu->getTypeID()) {fid=i; cV=pu->getVal(); break;}
                else if(use_conv)
                {
                    try {cV = unitconversions::convertPrimitive(pu->getVal(),pu->getTypeID(),typeIDs1[i],u.powers[idx]);}
                    catch(const NoPrimitiveConversionException& ex) {continue;}
                    fid=i;
                    break;
                }
            }
            if(fid==-1)
            {
                typeIDs1.push_back(pu->getTypeID());
                powers1.push_back(u.powers[idx]);
                vals1.push_back(pu->getVal());
            }
            else
            {
                powers1[fid]+=u.powers[idx];
                vals1[fid]*=cV;
            }
        }
        return StrippedUnit({typeIDs1,powers1,vals1});
    }
}


