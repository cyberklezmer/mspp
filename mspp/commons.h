#ifndef COMMONS_H
#define COMMONS_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <assert.h>
#include <limits>

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

static_assert(__cplusplus >= 201703, "C++17 language standard required");

static_assert(RAND_MAX >=2147483647, "RAND_MAX too small");


using namespace std;

namespace mspp
{

/// \addtogroup general General Definitions
///  @{

const double probabilitytolerance = 1e-15;


template <typename V>
const V inf() { return std::numeric_limits<V>::infinity(); }

template <typename V>
const V minf() { return -inf<V>(); }

// TBD
template <typename V>
const V max() { return 1e10; }

template <typename V> const V min() {
    return -1e10;
              }


template <typename T>
using ptr=std::shared_ptr<T>;

//template <typename T>
//using vector=std::vector<T>;


template <typename T>
class vectors : public vector<vector<T>>
{
public:
    vectors() {}
    vectors(unsigned int dim) : vector<vector<T>>(dim) {}
    vectors(const vector<vector<T>>& v): vector<vector<T>>(v) {}

    void converttoonedim(vector<T>& d)
    {
        d.clear();
        for(typename vectors<T>::iterator i = this->begin();
             i != this->end();
             i++)
            for(typename vector<T>::iterator j = i->begin();
                 j != i->end();
                 j++)
                d.push_back(*j);
    }

    operator vector<T>() const
    {
        vector<T> d;
        converttoonedim(d);
        return d;
    }
    unsigned int totaldim() const
    {
        unsigned int r=0;
        for(typename vectors<T>::iterator i = this->begin();
             i != this->end();
             i++)
            r += (*this)[i].size();
        return r;
    }
};


/// @}

/** \addtogroup sysclasses System Classes
 *  @{
 */


class object
{
public:
    virtual ~object() {}
};


class exception : public std::exception, public object
{
public:
    exception(const std::string &error = "unknown mspp error")
        : fmsg(error)
    {
    }

    const std::string& msg() { return fmsg; }
private:
    std::string fmsg;
};


class sys : public object
{
    static sys& self()
    {
       static sys s;
       return s;
    }
public:
    sys() : fout(0), flog(0), ferr(0)
    {
    }
    ~sys()
    {
    }

    static void setout(std::ostream& o)
    {
        self().fout = &o;
    }
    static void resetout()
    {
        self().fout = 0;
    }

    static void setlog(std::ostream& l)
    {
        self().flog = &l;
    }
    static void resetlog()
    {
        self().flog=0;
    }

    static void seterr(std::ostream& e)
    {
        self().ferr = &e;
    }
    static void reseterr()
    {
        self().ferr = 0;
    }

    static std::ostream& out()
    {
        return self().fout ? *(self().fout) : std::cout;
    }
    static std::ostream& log()
    {
        return self().flog ? *(self().flog) : std::clog;
    }
    static std::ostream& err()
    {
        return self().ferr ? *(self().ferr) : std::cerr;
    }
    static void seed()
    {
        srand(time(0));
    }

    static void seed(unsigned int seed)
    {
        srand(seed);
    }

    static double uniform()
    {
        return ((double) rand() + 0.5) / ((double) RAND_MAX + 1);
    }

private:
    std::ostream* fout;
    std::ostream* flog;
    std::ostream* ferr;
};

/// @}




/// \addtogroup general General Definitions
///  @{

/// \addtogroup pvar Variables
/// @{

using realvar=double;
using intvar=long;
using binvar=bool;

union mixedvar
{
    realvar x;
    intvar n;
    binvar b;
    mixedvar(const realvar& ax) : x(ax) {}
    mixedvar(const intvar& an) : n(an) {}
    mixedvar(const binvar& ab) : b(ab) {}
    operator realvar() const { return x; }
    operator intvar() const { return n; }
    operator binvar() const { return b; }
//    mixedvar& operator =(realvar& v) { x=v; return *this; }
//    mixedvar& operator =(intvar& v) { n=v; return *this; }
};

template <typename V>
using variables=vector<V>;

template <typename V>
class range : public object
{
public:
    using V_t = V;

    enum type { none, realt, intt, bint };

    range()
    {
        if constexpr( std::is_same<V,realvar>::value)
           ft = realt;
        else if constexpr( std::is_same<V,intvar>::value)
           ft = intt;
        else if constexpr( std::is_same<V,binvar>::value)
           ft = bint;
        else
           ft = none;

        setlimits();
    }


    range(type t, V l=minf<V>(), V h=inf<V>())
    {
        if constexpr( std::is_same<V,realvar>::value)
        {
           assert(t==realt);
           ft = realt;
        }
        else if constexpr( std::is_same<V,intvar>::value)
        {
            assert(t==intt);
            ft = intt;
        }
        else if constexpr( std::is_same<V,binvar>::value)
        {
            assert(t==bint);
            ft = bint;
        }
        else
        {
           assert(t==realt || t==intt || t==bint);
           ft = t;
        }
        setlimits(l,h);
    }

    void setreal(realvar l=minf<realvar>(), realvar h=inf<realvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = realt;
        setlimits(l,h);
    }

    void setint(intvar l=minf<intvar>(), intvar h=inf<intvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = intt;
        setlimits(l,h);
    }

    void setlimits(V l=minf<V>(), V h=inf<V>())
    {
        assert(ft==realt || ft==intt);
        fl=l; fh=h;
    }

    void setpositive()
    {
       setlimits(0);
    }
    V l() const
    {;
        assert(ft==realt || ft==intt);
        return fl;
    }
    V h() const
    {
        assert(ft==realt || ft==intt);
        return fh;
    }
    type t() const
    {
        assert(ft != none);
        return ft;
    }
    bool ishinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar h = fh;
            return h==inf<realvar>();
        }
        else
        {
            intvar h = fh;
            return h==inf<int>();
        }
    }
    bool islinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar l = fl;
            return l==-inf<realvar>();
        }
        else
        {
            intvar l = fl;
            return l==-inf<int>();
        }
    }
private:
    V fl;
    V fh;
    type ft;
};


///@}



/// \addtogroup cns Vector restrictions
/// @{

template <typename T>
using subvectors = vectors<T>;

class tilde: public object
{
public:
    bool included(unsigned int i, unsigned int s) const
    {
        assert(i<s);
        return is_included(i,s);
    }
    bool included(unsigned int i, unsigned int j, unsigned int s) const
    {
        assert(i<s);
        return is_included(i,j,s);
    }

    template <typename T>
    subvectors<T> operator () (const vectors<T>& v, unsigned int s)
    {
        vectors<T> d(v.size());
        unsigned int i=0;
        for( ;i<s; i++)
        {
            for(unsigned int j=0; j<v[i].size(); j++)
                if(included(i,j,s) )
                    d[i].push_back(v[i][j]);
        }
        for( ;i<v.size(); i++)
            d[i]=v[i];
        return d;
    }

    template <typename T>
    subvectors<T> operator () (const vectors<T>& v)
    {
        return operator ()(v,v.size());
    }

private:    
    virtual bool is_included(unsigned int i, unsigned int s) const = 0;
    virtual bool is_included(unsigned int i, unsigned int j, unsigned int s) const = 0;
};


class allx: public tilde
{
private:
    virtual bool is_included(unsigned int, unsigned int ) const { return true; }
    virtual bool is_included(unsigned int, unsigned int, unsigned int) const { return true; }
};

class lastx: public tilde
{
private:
    virtual bool is_included(unsigned int i, unsigned int s) const
       { return i==s-1; }
    virtual bool is_included(unsigned int i, unsigned int j, unsigned int s) const
       { return i==s-1; }
};


/// @}


///@}


/// \addtogroup problems Problems
/// @{

class constraint: public object
{
public:
    enum type {eq, geq, leq};
    constraint(unsigned int xdim, type t=eq)
        : fxdim(xdim), ft(t) {}

    type t() const { return ft; }
    void settype(type t) { ft = t; }
    unsigned int xdim() const { return fxdim; }
private:
    type ft;
    unsigned int fxdim;
};

const constraint::type eq = constraint::eq;
const constraint::type leq = constraint::leq;
const constraint::type geq = constraint::geq;

class criterion: public object
{
};

///@}


/// \addtogroup fns Functions
/// @{

class function: public object
{
public:
    function(unsigned int xdim) : fxdim(xdim) {}
    unsigned int xdim() const { return fxdim; }
    double operator() (const vector<double>& x) const
    {
        assert(x.size() == fxdim);
        return value_is(x);
    }
private:
    unsigned int fxdim;
    virtual double value_is(const vector<double>& x) const = 0;
};

class convexfunction: public function
{
public:
    convexfunction(unsigned int xdim) : function(xdim) {}
};

class subdiffunction: public convexfunction
{
public:
    subdiffunction(unsigned int dim) : convexfunction(dim) {}

    vector<double> sg(const vector<double> x) const
    {
        return sg_is(x);
    }
private:
   virtual vector<double> sg_is(const vector<double>& x) const = 0;
};

class linearfunction: public subdiffunction
{
public:
    linearfunction(unsigned int dim):
        subdiffunction(dim), fc(dim) {}

    linearfunction(vector<double> c):
        subdiffunction(c.size()), fc(c) {}

    void setc(const vector<double>& c)
    {
        assert(c.size()==xdim());
        fc = c;
    }
    void setc(unsigned int i, double c)
    {
        assert(i<fc.size());
        fc[i] = c;
    }
    double c(unsigned int i) const
    {
        assert(i < fc.size());
        return fc[i];
    }
    vector<double> c() const
    {
        return fc;
    }

private:
    virtual vector<double> sg_is(const vector<double>& x) const
    {
        assert(x.size()==xdim());
        return fc;
    }
    virtual double value_is(const vector<double>& x) const
    {
        assert(x.size()==xdim());
        double s=0;
        for(unsigned int i=0; i<fc.size(); i++)
            s += fc[i] * x[i];
        return s;
    }
private:
    vector<double> fc;
};

///@}




/// \ingroup Distributions
struct nocondition {};



/// \ingroup Processes
template <typename T>
using path=std::vector<T>;


/// \ingroup Processes
template <typename I=double>
using scenario = path<I>;

/// \ingroup Processes
template <typename C>
class zeta: public object
{
public:
    using C_t=C;
    virtual operator C () const = 0;
};

/// \ingroup Processes

template <typename I>
class allxi: public zeta <scenario<I>>, public scenario<I>
{
public:
    allxi(const scenario<I>& s) : scenario<I>(s)
    {
    }
    operator scenario<I>() const { return *this; }
};

/// \ingroup Processes

template <typename I>
class lastxi: public zeta<I>
{
public:
    lastxi(const scenario<I>& s)
    {
        assert(s.size());
        fxi = s[s.size()-1];
    }
    operator I() const { return fxi; }
private:
    I fxi;
};

/// \ingroup Processes

template <typename I>
class noxi: public zeta<nocondition>
{
public:
    noxi(const scenario<I>& s)
    {}
    operator nocondition() const { return nocondition(); }
};

} // namespace


/// \defgroup sms Solution Methods


#endif // COMMONS_H
