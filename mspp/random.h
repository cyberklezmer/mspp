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

template <typename I, typename J, typename C=nocondition>
class jointdistribution: virtual public distribution
        <pair<I, J>, C>
{
};


template <typename X=double, typename C=nocondition>
class vdistribution : virtual public distribution<vector<X>>
{
public:
    vdistribution(unsigned int adim) : dim(adim) {}
    const unsigned dim;
};

/*template <typename P, typename I=double, typename C=nocondition>
class parametricdistribution: virtual public distribution<I,C>
{
public:
    parametricdistribution(const P& p): fp(p) {}
protected:
    const P& p() { return fp; }
private:
    P fp;
};*/


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
class ddistribution: public distribution<I,C>,
                     public mcdistribution<I,C>
{
public:
    atom<I> operator () (unsigned int i, const C& c=C()) const
    {
        return atom_is(i,c);
    }
private:
    virtual atom<I> atom_is(unsigned int i, const C& c) const  = 0;
    virtual I do_draw(const C& c) const
    {
        double u = sys::uniform();
        double sum = 0;
        for(unsigned int i=0; ; i++)
        {
            atom<I> a = (*this)(i,c) ;
            sum += a.p;
            if(sum >= u)
                return a.x;
        }
        assert(0);
    }
};

template <typename I=double, typename C=nocondition>
class fdistribution: public ddistribution<I,C>
{
public:
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

private:
    virtual void atoms_are(vector<atom<I>>& a, const C& c) const
    {
        a.resize(0);
        unsigned int s= natoms_is(c);
        for(unsigned int i=0; i<s; i++)
            a.push_back((*this)(i,c));
    }
    virtual unsigned int natoms_is(const C& c) const = 0;
};



/*
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
*/

template <typename I=double>
class ldistribution: public fdistribution<I,nocondition>
{   
public:
    ldistribution(const vector<atom<I>>& atoms) :
        fatoms(atoms)
    {
    }

    ldistribution(const vector<I>& values) :
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
private:
    virtual void atoms_are(vector<atom<I>>& a, const nocondition&) const
    {
        a = fatoms;
    };

    virtual unsigned int natoms_is(const nocondition&) const
    {
        return fatoms.size();
    }

    virtual atom<I> atom_is(unsigned int i, const nocondition&) const
    {
        return fatoms[i];
    }
private:
    vector<atom<I>> fatoms;
};

template <typename I>
class diracdistribution: public ldistribution<I>
{
public:
    diracdistribution(const I& a)
        : ldistribution<I>({a}) {}
};


/*
template <typename X=double>
class gmdddistribution: public ldistribution<rvector<X>>,
     public mddistribution<X,nocondition>
{
public:
    gmdddistribution(const vector<rvector<X>>& values) :
        ldistribution<rvector<X>>(values),
        mddistribution<X>(values[0].size())
    {
    }
    gmdddistribution(const vector<atom<rvector<X>>>& atoms) :
        ldistribution<rvector<X>>(atoms),
        mddistribution<X>(atoms[0].x.size())
     {}
};*/


template <typename I>
ldistribution<vector<I>> operator *(const fdistribution<I,nocondition>& x,
                             const fdistribution<I,nocondition>& y)
{
    assert(x.natoms(nocondition()));
    assert(y.natoms(nocondition()));
    vector<atom<vector<I>>> a(x.natoms(nocondition())*y.natoms(nocondition()));

    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(nocondition()); i++)
        for(unsigned int j=0; j<y.natoms(nocondition()); j++)
        {
            atom<I> xa = x(i,nocondition());
            atom<I> ya = y(j,nocondition());
            vector<I> r = { xa.x, ya.x };
            a[dst++]= { r, xa.p*ya.p};
        }
    return ldistribution<vector<I>>(a);
}



/// @} - discrete distributions

/// \addtogroup iterativedists Iteratively defined distributions
/// \ingroup Distributions
/// @{



template <typename D, typename E,
          typename Z=nomapping<pair<typename D::C_t,typename D::I_t>>>
class ipdistribution:
     public jointdistribution<typename D::I_t, typename E::I_t, typename D::C_t>
{
public:
//    using I_t = pair<typename D::I_t, typename E::I_t>;
//    using C_t = typename D::C_t;
    ipdistribution(const D& d, const E& e) : fd(d), fe(e)
    {
        static_assert(std::is_convertible<typename Z::R_t&,typename E::C_t&>::value);
    }
    const D& d() const { return fd; }
    const E& e() const { return fe; }
private:
    D fd;
    E fe;
};

template <typename D, typename E, typename Z>
class mcipdistribution: public ipdistribution<D,E,Z>
{
public:
    using I_t = typename ipdistribution<D,E>::I_t;
    using C_t = typename ipdistribution<D,E>::C_t;
    mcipdistribution(const D& d, const E& e) :
      ipdistribution<D,E,Z>(d,e)
    {
    }
    virtual I_t do_draw(const C_t& c) const
    {
        typename D::I_t i = this->d().draw(c);
        pair<C_t,typename D::I_t> p = { c, i };
        typename E::I_t j = this->e().draw(Z()(p));
        return {i,j};
    }
};



template <typename D,
          typename Z=nomapping<vector<typename D::I_t>>>
class ivdistribution:
     public vdistribution<typename D::I_t,typename Z::R_t>
{
public:
    ivdistribution(const vector<D>& d)
       : vdistribution<typename D::I_t,typename Z::R_t>(d.size())
    {
        assert(d.size());
        static_assert(std::is_same<typename Z::D_t,
                           std::vector<typename D::I_t>>::value);
    }
    const vector<D>& d() const { return fd; }
private:
    vector<D> fd;
};


/*
template <typename D, typename E, typename Z>
class dipdistribution:
        public ipdistribution<D,E,Z>,
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

    virtual I_t do_draw(const C_t& c) const
    {
        I_t r;
        r.i = this->d().draw(c);
        r.j = this->e().draw(r.i);
        return r;
    }
};

*/

/// @} - iterative distributions



/// @} - distribution


} // namespace


#endif // RANDOM_H
