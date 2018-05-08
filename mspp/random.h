#ifndef RANDOM_H
#define RANDOM_H

#include "mspp/commons.h"

namespace mspp
{

/// \addtogroup Distributions
/// @{

using probability = double;

template <typename I>
struct atom
{
    using I_t = I;
    I x;
    probability p;
};

template <typename X=double, typename R=nothing>
class distribution: public object
{
public:
    distribution(unsigned int dim) : fdim(dim) {}
    unsigned int dim() const { return fdim; }
    using X_t = X;
    using R_t = R;

    static condition<X> c(const scenario<X>& s) { R r; return r(s); }
private:
    unsigned int fdim;
};


template <typename X, typename R, typename P>
class parametricdistribution: virtual public distribution<X,R>
{
public:
    parametricdistribution(const P& p): fp(p) {}
protected:
    const P& p() { return fp; }
private:
    P fp;
};


/// \addtogroup mcdisc MC distributions
/// \ingroup Distributions
/// @{

template <typename X, typename R>
class mcdistribution: virtual public distribution<X, R>
{
public:
    rvector<X> draw(const condition<X>& c) const
    {
        rvector<X> r(this->dim());
        draw(c,r);
        assert(r.size()==this->dim());
        return r;
    }

//    rvector<X> draw(const scenario<X>& s) const
//    {
//        return draw(distribution<X,R>::c(s));
//    }
private:
    virtual void draw(const condition<X>& c, rvector<X>& r ) const = 0;
};

template <typename X>
class imcdistribution: public mcdistribution<X, nothing>
{
public:
    rvector<X> draw() const
    {
        rvector<X> r(this->dim());
        draw(r);
        assert(r.size()==this->dim());
        return r;
    }

private:
    virtual void draw(const condition<X>& c, rvector<X>& r) const
    {
        draw(r);
    }

    virtual void draw(rvector<X> &r ) const = 0;
};

/// @}

/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{

template <typename X>
using rvatom=atom<rvector<X>>;

template <typename X, typename R>
class ddistribution: virtual public distribution<X,R>
{
public:
    ddistribution() {}

    void atoms(vector<rvatom<X>>& a, const condition<X>& c) const
    {
        a.resize(0);
        atoms_are(a,c);
        probability p=0;
        for(unsigned int i=0; i< a.size(); i++)
        {
            assert(a[i].x.size() == this->dim());
            p += a[i].p;
        }

        assert(fabs(p-1.0) < 0.000000000001);
    }

    unsigned int natoms(const condition<X>& c) const
    {
        return natoms_is(c);
    }

    void atoms(vector<rvatom<X>>& a, const scenario<X>& s) const
    {
        atoms(a,distribution<X,R>::c(s));
    }

    rvatom<X> atom(unsigned int i, const condition<X>& c) const
    {
        assert(i<natoms_is(c));
        return atom_is(i,c);
    }

protected:
    virtual void atoms_are(vector<rvatom<X>>& a, const condition<X>& c) const
    {
        a.resize(0);
        unsigned int s= natoms_is(c);
        for(unsigned int i=0; i<s; i++)
            a.push_back(atom_is(i,c));
    }
    virtual rvatom<X> atom_is(unsigned int i, const condition<X>& c) const  = 0;
    virtual unsigned int natoms_is(const condition<X>& c) const = 0;
};

template <typename X>
class iddistribution: virtual public ddistribution<X,nothing>
{
public:
    iddistribution() {}

    rvatom<X> operator [] (unsigned int i) const
    {
        assert(i < natoms_is());
        return atom_is(i);
    }

    void atoms(vector<rvatom<X>>& a) const
    {
        a.clear();
        atoms_are(a);
    }

    unsigned int natoms() const
    {
        return natoms_is();
    }

protected:
    virtual void atoms_are(vector<rvatom<X>>& a, const condition<X>&) const
      { atoms_are(a); }
    virtual void atoms_are(vector<rvatom<X>>& a) const
    {
        a.resize(0);
        unsigned int s= natoms_is();
        for(unsigned int i=0; i<s; i++)
            a.push_back(atom_is(i));
    }

    virtual unsigned int natoms_is(const condition<X>&) const
      { return natoms_is(); }
    virtual unsigned int natoms_is() const = 0;
    virtual rvatom<X> atom_is(unsigned int i, const condition<X>&) const
    {
        return atom_is(i);
    }
    virtual rvatom<X> atom_is(unsigned int i) const = 0;
};


template <typename X=double>
class gddistribution: public iddistribution<X>
{   
public:
    gddistribution(const vector<rvatom<X>>& atoms) :
        distribution<X>(atoms[0].x.size()), fatoms(atoms)
    {
    }

    gddistribution(const iddistribution<X>& d) :
       fatoms(d.dim())
    {
       d.atoms(fatoms);
    }

    gddistribution(const vector<X>& values) :
        distribution<X>(1), fatoms(values.size())
    {
        assert(fatoms.size());
        probability p=1.0/fatoms.size();
        for(unsigned int i=0; i<fatoms.size(); i++)
        {
            fatoms[i].x = {values[i]};
            fatoms[i].p = p;
        }
    }

    gddistribution(const vector<rvector<X>>& values) :
        distribution<X>(values[0].size()), fatoms(values.size())
    {
        assert(fatoms.size());
        probability p=1.0/fatoms.size();
        for(unsigned int i=0; i<fatoms.size(); i++)
        {
            fatoms[i].x = values[i];
            fatoms[i].p = p;
            assert(fatoms[i].x.size() == fatoms[0].x.size());
        }
    }

protected:
    virtual void atoms_are(vector<rvatom<X>>& a) const
    {
        a = fatoms;
    };

    unsigned int natoms_is() const
    {
        return fatoms.size();
    }

    virtual rvatom<X> atom_is(unsigned int i) const
    {
        return fatoms[i];
    }

private:
    vector<rvatom<X>> fatoms;
};

template <typename X>
gddistribution<X> operator *(const iddistribution<X>& x,
                             const iddistribution<X>& y)
{
    assert(x.natoms());
    assert(y.natoms());
    vector<rvatom<X>> a(x.natoms()*y.natoms());

    unsigned int d=x.dim()+y.dim();
    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(); i++)
        for(unsigned int j=0; j<y.natoms(); j++)
        {
            rvatom<X> xa = x[i];
            rvatom<X> ya = y[j];
            assert(xa.x.size()+ya.x.size()==d);
            rvector<X> r = xa.x;
            for(unsigned int k=0; k<ya.x.size(); k++)
                r.push_back(ya.x[k]);
            assert(r.size()==d);
            a[dst++]= { r, xa.p*ya.p};
        }
    return gddistribution<X>(a);
}



/// @} - discrete distributions


/// @} - distribution

/// \addtogroup processes Proceses
/// @{

template <typename D>
class processdistribution : public object
{
public:
    processdistribution(const rvector<typename D::X_t>& xi0,
                        const D& d, unsigned int T)
      : fxi0(xi0), fd(T,d)
    {
    }
    processdistribution(const D& d, unsigned int T)
      : fxi0(0), fd(T,d)
    {
    }

    processdistribution(const rvector<typename D::X_t>& xi0, const vector<D>& d)
      : fxi0(xi0), fd(d)
    {
    }

    rvector<typename D::X_t> xi0() const { return fxi0; }
    D d(unsigned int k) const { assert(k<=fd.size()); return fd[k-1]; }
    unsigned int T() const { return fd.size(); }
private:
    rvector<typename D::X_t> fxi0;
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
class   indexedpath : public path<indexedatom<I>>
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


template<typename X, typename S=void>
using sctreecallback=ptreecallback<rvector<X>,S>;


template <typename X>
class scenariotree: public ptree<rvector<X>>
{
public:
    using X_t = X;
};

template <typename D>
class distrscenariotree : public scenariotree<typename D::X_t>
{
    using X=typename D::X_t;
public:
    distrscenariotree(const processdistribution<D>& p):fp(p)
    {}
    distrscenariotree(const D& d, unsigned int T) :
        fp(processdistribution<D>(d,T))
    {}
    distrscenariotree(const rvector<X>& xi0, const D& d, unsigned int T) :
        fp(processdistribution<D>(xi0,d,T))
    {}
private:

    virtual void branches_are(vector<rvatom<X>>& bchs,
                              const indexedpath<rvector<X>>& s) const
    {
        unsigned int k=s.size();
        assert(k);
        assert(fp.T());
        fp.d(k).ddistribution<typename D::X_t,typename D::R_t>::atoms(bchs,s.pth());
    }
    virtual rvector<X> root_is() const
    {
        return fp.xi0();
    }
    virtual unsigned int  T_is() const { return fp.T(); }
private:
    processdistribution<D> fp;
};

template <typename X>
using gdscenariotree=distrscenariotree<gddistribution<X>>;

/// @} - processes


} // namespace


#endif // RANDOM_H
