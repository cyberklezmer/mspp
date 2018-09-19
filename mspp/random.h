#ifndef RANDOM_H
#define RANDOM_H

#include "mspp/commons.h"
#include <cmath>

namespace mspp
{

/// \addtogroup Distributions
/// @{

using probability = double;

template <typename I=double, typename C=nocondition>
class distribution: public object
{
public:
    distribution() {}
    using I_t = I;
    using C_t = C;
};

template<typename X=double>
using rvector = vector<X>;


template <typename X=double, typename C=nocondition>
class mddistribution : virtual public distribution<rvector<X>>
{
public:
    mddistribution(unsigned int adim) : dim(adim) {}
    const unsigned dim;
};

template <typename P, typename I=double, typename C=nocondition>
class parametricdistribution: virtual public distribution<I,C>
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

template <typename I=double, typename C=nocondition>
class mcdistribution: virtual public distribution<I, C>
{
public:
    I draw(const C& c) const
    {
        return do_draw(c);
    }

private:
    virtual I do_draw(const C& c) const = 0;
};


template <typename I=double>
class imcdistribution: public mcdistribution<I, nocondition>
{
public:
    I draw() const
    {
        return draw();
    }

private:
    virtual I do_draw(const nocondition& c) const
    {
        return do_draw();
    }


    virtual I do_draw() const = 0;
};


template <typename D, typename E>
class iterativedistribution:
    virtual public distribution
        <pair<typename D::I_t, typename E::I_t>, typename D::C_t>
{
public:
    using I_t = pair<typename D::I_t, typename E::I_t>;
    using C_t = typename D::C_t;
    iterativedistribution(const D& d, const E& e) : fd(d), fe(e)
    {
        static_assert(std::is_same<typename E::C_t, typename D::I_t>::value);
    }
    const D& d() const { return fd; }
    const E& e() const { return fe; }
private:
    D fd;
    E fe;
};

template <typename D, typename E>
class mciterativedistribution: public iterativedistribution<D,E>
{
public:
    using I_t = typename iterativedistribution<D,E>::I_t;
    using C_t = typename iterativedistribution<D,E>::C_t;
    mciterativedistribution(const D& d, const E& e) :
      iterativedistribution<D,E>(d,e)
    {
    }
    virtual I_t do_draw(const C_t& c) const
    {
        I_t r;
        r.i = this->d().draw(c);
        r.j = this->e().draw(r.i);
        return r;
    }
};


/// @}

/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{

template <typename I=double>
struct atom
{
    using I_t = I;
    I x;
    probability p;
};

template <typename I=double, typename C=nocondition>
class ddistribution: virtual public distribution<I,C>,
                     virtual public mcdistribution<I,C>
{
public:
    ddistribution() {}

    void atoms(vector<atom<I>>& a, const C& c) const
    {
        a.resize(0);
        atoms_are(a,c);
        probability p=0;
        for(unsigned int i=0; i< a.size(); i++)
        {
            p += a[i].p;
        }

        assert(fabs(p-1.0) < probabilitytolerance);
    }

    unsigned int natoms(const C& c) const
    {
        return natoms_is(c);
    }


    atom<I> operator () (unsigned int i, const C& c) const
    {
        assert(i<natoms_is(c));

        return atom_is(i,c);
    }

private:
    virtual void atoms_are(vector<atom<I>>& a, const C& c) const
    {
        a.resize(0);
        unsigned int s= natoms_is(c);
        for(unsigned int i=0; i<s; i++)
            a.push_back(atom_is(i,c));
    }
    virtual atom<I> atom_is(unsigned int i, const C& c) const  = 0;
    virtual unsigned int natoms_is(const C& c) const = 0;

    virtual I do_draw(const C& c) const
    {
        double u = sys::uniform();
        unsigned int n = natoms(c);
        double sum = 0;
        for(unsigned int i=0; i<n; i++)
        {
            atom<I> a = (*this)(i,c) ;
            sum += a.p;
            if(sum >= u)
                return a.x;
        }
        assert(0);
    }
};

template <typename I>
class iddistribution: virtual public ddistribution<I,nocondition>
{
public:
    iddistribution() {}

    atom<I> operator [] (unsigned int i) const
    {
        assert(i < natoms_is());
        return atom_is(i);
    }

    void atoms(vector<atom<I>>& a) const
    {
        a.clear();
        atoms_are(a);
    }

    unsigned int natoms() const
    {
        return natoms_is();
    }

protected:
    virtual void atoms_are(vector<atom<I>>& a, const nocondition&) const
      { atoms_are(a); }
    virtual void atoms_are(vector<atom<I>>& a) const
    {
        a.resize(0);
        unsigned int s= natoms_is();
        for(unsigned int i=0; i<s; i++)
            a.push_back(atom_is(i));
    }

    virtual unsigned int natoms_is(const nocondition&) const
      { return natoms_is(); }
    virtual unsigned int natoms_is() const = 0;
    virtual atom<I> atom_is(unsigned int i, const nocondition&) const
    {
        return atom_is(i);
    }
    virtual atom<I> atom_is(unsigned int i) const = 0;
};


template <typename I=double>
class gddistribution: public iddistribution<I>
{   
public:
    gddistribution(const vector<atom<I>>& atoms) :
        fatoms(atoms)
    {
    }

    gddistribution(const iddistribution<I>& d)
    {
       d.atoms(fatoms);
    }

    gddistribution(const vector<I>& values) :
        fatoms(values.size())
    {
        assert(fatoms.size());
        probability p=1.0/fatoms.size();
        for(unsigned int i=0; i<fatoms.size(); i++)
        {
            fatoms[i].x = values[i];
            fatoms[i].p = p;
        }
    }


protected:
    virtual void atoms_are(vector<atom<I>>& a) const
    {
        a = fatoms;
    };

    unsigned int natoms_is() const
    {
        return fatoms.size();
    }

    virtual atom<I> atom_is(unsigned int i) const
    {
        return fatoms[i];
    }

private:
    vector<atom<I>> fatoms;
};

template <typename X=double>
class gmdddistribution: public gddistribution<rvector<X>>,
     public mddistribution<X,nocondition>
{
public:
    gmdddistribution(const vector<rvector<X>>& values) :
        gddistribution<rvector<X>>(values),
        mddistribution<X>(values[0].size())
    {
    }
    gmdddistribution(const vector<atom<rvector<X>>>& atoms) :
        gddistribution<rvector<X>>(atoms),
        mddistribution<X>(atoms[0].x.size())
     {}


};


template <typename I>
gmdddistribution<I> operator *(const iddistribution<I>& x,
                             const iddistribution<I>& y)
{
    assert(x.natoms());
    assert(y.natoms());
    vector<atom<vector<I>>> a(x.natoms()*y.natoms());

    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(); i++)
        for(unsigned int j=0; j<y.natoms(); j++)
        {
            atom<I> xa = x[i];
            atom<I> ya = y[j];
            rvector<I> r = { xa.x, ya.x };
            a[dst++]= { r, xa.p*ya.p};
        }
    return gmdddistribution<I>(a);
}


template <typename D, typename E>
class diterativedistribution:
        public iterativedistribution<D,E>,
        public ddistribution<pair<typename D::I_t,typename E::I_t>,
                                               typename D::C_t>

{
public:
    using I_t = typename iterativedistribution<D,E>::I_t;
    using C_t = typename iterativedistribution<D,E>::C_t;
    diterativedistribution(const D& d, const E& e) :
        iterativedistribution<D,E>(d,e)
    {
    }
private:
    virtual void atoms_are(vector<atom<I_t>>& a, const C_t& c) const
    {
        a.clear();
        vector<atom<typename D::I_t>> da;
        this->d().atoms(da,c);
        unsigned int nd=da.size();
        for(unsigned int i=0; i<nd; i++)
        {
            atom<typename D::I_t> x = da[i];

            vector<atom<typename E::I_t>> ea;
            this->e().atoms(ea,x.x);
            unsigned int ne=ea.size();

            for(unsigned int j=0; j<ne; j++)
            {
                atom<typename E::I_t> y = ea[j];
                a.push_back({{x.x,y.x},x.p*y.p});
            }
        }
    }


    virtual atom<I_t> atom_is(unsigned int i, const C_t& c) const
    {
        vector<atom<I_t>> a;
        atoms_are(a,c);
        assert(i<a.size());
        return a[i];
    }
    virtual unsigned int natoms_is(const C_t& c) const
    {
        vector<atom<I_t>> a;
        atoms_are(a,c);
        return a.size();
    }
};


/// @} - discrete distributions


/// @} - distribution


} // namespace


#endif // RANDOM_H
