#ifndef PROCESS_H
#define PROCESS_H

#include "mspp/random.h"

namespace mspp
{




/// \addtogroup processes Proceses
/// @{



template <typename D, typename Z=noxi<typename D::I_t>>
class processdistribution : public object
{
public:
    using Z_t = Z;
    using D_t = D;

    processdistribution(const typename D::I_t& xi0,
                        const D& d, unsigned int T)
      : fxi0(xi0), fd(T,d)
    {
    }
    processdistribution(const D& d, unsigned int T)
      : fxi0(0), fd(T,d)
    {
    }

    processdistribution(const typename D::I_t& xi0, const vector<D>& d)
      : fxi0(xi0), fd(d)
    {
    }

    typename D::I_t xi0() const { return fxi0; }
    const D& d(unsigned int k) const { assert(k); assert(k<=fd.size()); return fd[k-1]; }
    unsigned int T() const { return fd.size(); }
    typename D::C_t c(const scenario<typename D::I_t>& s) const
    {
        return Z(s);
    }
private:
    typename D::I_t fxi0;
    vector<D> fd;
};


template <typename I>
struct indexedatom{
//    indexedatom(const I& ax, probability ap, unsigned int ai) :
//        x(ax), p(ap), i(ai)
//    {}
    I x;
    probability p;
    unsigned int i;
};




template <typename I>
class indexedpath : public path<indexedatom<I>>
{
public:
    double uncprob() const
    {
        assert(this->size());
        assert((*this)[0].p==1);
        double r = 1.0;
        for(unsigned int i=1; i<this->size(); i++)
            r *= (*this)[i].p;
        return r;
    }

    path<I> pth() const
    {
        path<I> r(this->size());
        for(unsigned int i=0; i<this->size(); i++)
            r[i]=(*this)[i].x;
        return r;
    }
};


template<typename I, typename S=void>
class ptreecallback
{
public:
    virtual void callback(const indexedpath<I>& path, S* state=0) const = 0;
};

template <typename I>
class ptree: public object
{
public:
    using I_t = I;
public:
    unsigned int T() const { return T_is(); }

    template <typename S=void>
    void foreachnode(const ptreecallback<I,S> *callee, S* state=0) const
    {
        indexedpath<I> s;
        s.push_back({root_is(),1,0});
        callee->callback(s, state);
        if(T_is())
            doforeachnode(callee, s, state);
    }
private:
    virtual unsigned int T_is() const = 0;
    virtual I root_is() const = 0;
    virtual void branches_are(vector<atom<I>>& bchs, const indexedpath<I>& s) const = 0;

    template <typename S=void>
    void doforeachnode(const ptreecallback<I,S> *callee,
                          const indexedpath<I>& as,
                          S* state) const
    {
        unsigned int k=as.size();
        indexedpath<I> s(as);
        vector<atom<I>> itms;
        branches_are(itms,s);
        s.resize(k+1);
        bool cb = k < this->T();
        unsigned int i=0;
        for(typename vector<atom<I>>::iterator it = itms.begin();
              it != itms.end();
              it++,i++)
        {
            s[k]={it->x,it->p, i};
            callee->callback(s,state);
            if(cb)
                doforeachnode(callee, s, state);
        }
    }
};


template<typename I, typename S=void>
using sctreecallback=ptreecallback<I,S>;


template <typename I=double>
class scenariotree: public ptree<I>
{
public:
    using I_t = I;
};

template <typename D,typename Z=noxi<typename D::I_t>>
class distrscenariotree : public scenariotree<typename D::I_t>
{
    using I=typename D::I_t;
public:
    distrscenariotree(const processdistribution<D,Z>& p):fp(p)
    {}
    distrscenariotree(const D& d, unsigned int T) :
        fp(processdistribution<D,Z>(d,T))
    {}
    distrscenariotree(const I& xi0, const D& d, unsigned int T) :
        fp(processdistribution<D,Z>(xi0,d,T))
    {}

    distrscenariotree(const I& xi0, const vector<D>&d) :
        fp(processdistribution<D,Z>(xi0,d))
    {}

private:
    virtual void branches_are(vector<atom<I>>& bchs,
                              const indexedpath<I>& s) const
    {
        unsigned int k=s.size();
        assert(k);
        assert(fp.T());

        typename D::C_t c = fp.c(s.pth());
        fp.d(k).ddistribution<typename D::I_t,typename D::C_t>::atoms(bchs,c);
    }
    virtual typename D::I_t root_is() const
    {
        return fp.xi0();
    }
    virtual unsigned int  T_is() const { return fp.T(); }
private:
    processdistribution<D,Z> fp;
};

template <typename X>
using gdscenariotree=distrscenariotree<gddistribution<X>>;

template <typename X>
using gmddscenariotree=distrscenariotree<gmdddistribution<X>>;


/// @} - processes

} // namespace

#endif // PROCESS_H
