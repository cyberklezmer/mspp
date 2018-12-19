#ifndef PROCESS_H
#define PROCESS_H

#include "mspp/random.h"

namespace mspp
{

/// \addtogroup processes Proceses
/// @{

/// the distribution as a whole is not conditioned
template <typename D, typename Z>
class processdistribution :
    public uijdistribution<
         diracdistribution<typename D::I_t>,
         ivdistribution<D,Z>,Z>,
    virtual public vdistribution<typename D::I_t,novalue>
{
public:
    using X_t = typename D::I_t;
    using D_t = D;
    using Z_t = Z;
    processdistribution(const X_t& xi0, const D& d, unsigned int T) :
      uijdistribution<diracdistribution<X_t>,ivdistribution<D,Z>,Z>
         (diracdistribution<X_t>(xi0),ivdistribution<D,Z>(d,T)),
         vdistribution<typename D::I_t,novalue>(T)
    {
        static_assert(std::is_base_of<
             zeta<typename Z::X_t,typename Z::R_t>,Z>::value);
        assert(T>=1);
    }
    processdistribution(const X_t& xi0, const vector<D>& d) :
      uijdistribution<diracdistribution<X_t>,ivdistribution<D,Z>,Z>
         (diracdistribution<X_t>(xi0),ivdistribution<D,Z>(d)),
      vdistribution<typename D::I_t,novalue>(d.size()+1)
    {
        static_assert(std::is_base_of<
             zeta<typename Z::X_t,typename Z::R_t>,Z>::value);
        assert(d.size());
    }
    unsigned int T() const { return this->second().dim();}
    const D& d(unsigned int i) const
    { assert(i<=T()); assert(i>0); return this->second().d(i-1); }
    const X_t x0() const
    {
        diracdistribution<X_t> dd=this->first();
        X_t x = dd.x();
        return x;
    }
};

template<typename X, typename S=void>
class tdcallback : public object
{
public:
    virtual void callback(const scenario<X>& path,
                          probability up, S* state=0) const = 0;
};


template <typename X, typename E=atom<X>, typename A=novalue>
class treedistribution: virtual public vdistribution<X,novalue>
{
public:
    using X_t = X;
    using E_t = E;
    using A_t = A;

    treedistribution(unsigned int T) : vdistribution<X,novalue>(T+1) {}

    template <typename S=void>
    void foreachnode(const tdcallback<X, S> *callee, S* state=0) const
    {
        vector<E> e;
        scenario<X> s;
        A mystate;
        beforefe(mystate);
        doforeachnode(callee, e, s, state,mystate);
        afterfe(mystate);
    }

    void branches(const vector<E>& e, vector<E>& es) const
    {
        branches_are(e, es);
    }
private:

    template <typename S=void>
    void doforeachnode(const tdcallback<X,S> *callee,
                       const vector<E>& ae,
                       const scenario<X> as,
                       S* state,
                       A& mystate) const
    {
        unsigned int k=ae.size();
        vector<E> e(ae);
        vector<E> es;

        branches_are(e,es);

        scenario<X> s(as);
        s.resize(k+1);
        e.resize(k+1);
        bool rec = k < this->dim()-1;

        for(unsigned int i=0;i<es.size();i++)
        {
            e[k]=es[i];
            probability up;
            index2sinfo(e,s[k],up,mystate);
            callee->callback(s,up,state);
            if(rec)
                doforeachnode(callee, e, s, state, mystate);
        }
    }

    virtual void branches_are(const vector<E>& e, vector<E>& es) const = 0;
    virtual void beforefe(A&) const {}
    virtual void afterfe(A&) const {}
protected:
    virtual void index2sinfo(const vector<E>& e,
                                X& s,
                                probability& up,
                                A& ) const
    {
       if constexpr(std::is_same<E,atom<X_t>>::value)
       {
           assert(e.size());
           up = 1;
           for(unsigned int i=0; i<e.size(); i++)
               up *= e[i].p;
           s = e[e.size()-1].x;
       }
       else
            assert(0);
    }
};


template <typename D, typename Z>
class fdprocessdistribution :
        public processdistribution<D,Z>,
        public treedistribution<typename D::I_t>
{
public:
    using X_t=typename D::I_t;
    fdprocessdistribution(const X_t& xi0, const D& d, unsigned int T) :
       treedistribution<X_t>(T),
       processdistribution<D, Z> (xi0, d, T),
       vdistribution<X_t,novalue>(T+1)
    {
        static_assert(
           std::is_base_of<fdistribution<X_t,typename D::C_t,
                    D::flistdef>,D>::value
                    );
        assert(T);
    }
    fdprocessdistribution(const X_t& xi0, const vector<D>& d) :
       treedistribution<X_t>(d.size()+1),
       processdistribution<D, Z> (xi0, d),
       vdistribution<X_t,novalue>(d.size()+1)
    {
        static_assert(
           std::is_base_of<fdistribution<X_t,typename D::C_t>,D>::value
                    );
    }
private:
    virtual void branches_are(const vector<atom<X_t>>& e,
         vector<atom<X_t>>& es) const
    {
        unsigned int k=e.size();
        assert(this->T());
        if(k)
        {
            scenario<X_t> s;
            for(unsigned int i=0; i<k; i++)
                s.push_back(e[i].x);
            this->second().atoms(s,es);
        }
        else
            this->first().atoms(es,novalue());
    }
};


/// @} - processes

} // namespace

#endif // PROCESS_H
