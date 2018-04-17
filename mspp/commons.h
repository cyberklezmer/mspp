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

namespace mspp
{

/// \addtogroup general General Definitions
///  @{


const double inf = std::numeric_limits<double>::infinity();
const double minf = -inf;

template <typename T>
using ptr=std::shared_ptr<T>;

using rn=std::vector<double>;

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


/*class sys : public object
{
    static sys* fself;
public:
    sys() { fself=this; }
    ~sys() { fself = 0; }
    static void set_log(const ptr<std::ostream>& l)
    {
        assert(fself);
        fself->flog = l;
    }
    static void reset_log()
    {
        assert(fself);
        fself->flog.reset();
    }

    static void set_err(const ptr<std::ostream>& e)
    {
        assert(fself);
        fself->ferr = e;
    }
    static void reset_err()
    {
        assert(fself);
        fself->ferr.reset();
    }

    static std::ostream& log()
    {
        assert(fself);
        return fself->flog ? *(fself->flog) : std::cout;
    }
    static std::ostream& err()
    {
        assert(fself);
        return fself->ferr ? *(fself->ferr) : std::cerr;
    }
private:
    ptr<std::ostream> flog;
    ptr<std::ostream> ferr;
};*/

/// @}



/// \addtogroup general General Definitions
///  @{


template <typename X>
using scenario = std::vector<X>;

template <typename X>
class condition: public object
{
public:
    using X_t = X;
};

/// \addtogroup cns Conditions
/// \ingroup general
/// @{



template <typename X>
class emptycondition: public condition<X>
{
public:
    emptycondition(const scenario<X>&) {}
};

template <typename X>
class fullhistory: public condition<X>
{
public:
    fullhistory(const scenario<X>& s): fs(s)
    {}
    const scenario<X>& s() const { return fs; }
    const X& operator[](unsigned int k) const
                    { assert(k<fs.size()); return fs[k]; }
private:
    scenario<X> fs;
};

template <typename X>
class lastvalue: public condition<X>
{
    unsigned int checksize(unsigned int s) const
    {
        assert(s);
        return s;
    }
public:
    lastvalue(const scenario<X>& s):
        fx(s[checksize(s.size())-1])
    {
    }
    X x() const { return fx; }
    operator X() const {return x();}
private:
    X fx;
};

template <typename X, typename M>
class markovcondition: public condition<X>
{
public:
    struct result { X x; M m; };
    markovcondition(const scenario<X>& s):
        fs(s)
    {
        assert(s.size());
    }
    X x() const { return fs[fs.size()-1]; }
    M m() const { assert(fs.size() > 1); return m_is(fs); }
    result r() const { return { x(), m() }; }
    operator result() const {return r();}
private:
    virtual M m_is(const scenario<X> fs ) const = 0;
    scenario<X> fs;
};


/// @}

/// \addtogroup pvar Variables
/// @{



using variable=double;

using variables=std::vector<variable>;

class vardef : public object
{
};

template <typename V>
using vardefs = std::vector<V>;

class realvar : public vardef
{
public:
    enum type { R, Rplus, Rminus };

    realvar(type t=R)
    {
        set(t);
    }

    void set(double l=minf, double h=inf)
    { fl=l; fh=h; }

    void set(type at)
    {
        switch(at)
        {
        case R: fl=minf; fh=inf; break;
        case Rplus: fl=0; fh=inf; break;
        case Rminus: fl=minf; fh=0; break;
        }
    }
    double l() const { return fl; }
    double h() const { return fh; }
private:
    double fl;
    double fh;
};

///@}

///@}


/// \addtogroup problems Problems
/// @{

class constraint: public object
{
public:
    enum type {eq, geq, leq};
    constraint(type t=eq) : ft(t) {}

    type t() const { return ft; }
    void settype(type t) { ft = t; }
private:
    type ft;
};

///@}


/// \addtogroup fns Functions
/// @{

class function: public object
{
public:
    function(unsigned int dim) : fdim(dim) {}
    unsigned int dim() const { return fdim; }
    double operator() (const std::vector<variable>& x) const
    {
        assert(x.size() == fdim);
    }
private:
    unsigned int fdim;
    virtual double value_is(const std::vector<variable>& x) const = 0;
};

class convexfunction: public function
{
public:
    convexfunction(unsigned int dim) : function(dim) {}
};

class sgcfunction: public convexfunction
{
public:
    sgcfunction(unsigned int dim) : convexfunction(dim) {}

    std::vector<double> sg(const std::vector<variable> x) const
    {
        return sg_is(x);
    }
private:
   virtual std::vector<double> sg_is(const std::vector<variable>& x) const = 0;
};

class linearfunction: public sgcfunction
{
public:
    linearfunction(unsigned int dim):
        sgcfunction(dim), fc(dim) {}

    linearfunction(std::vector<double> c):
        sgcfunction(c.size()), fc(c) {}

    void setc(const std::vector<double>& c)
    {
        assert(c.size()==dim());
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
    std::vector<double> c() const
    {
        return fc;
    }

private:
    virtual std::vector<double> sg_is(const std::vector<variable>& x) const
    {
        assert(x.size()==dim());
        return fc;
    }
    virtual double value_is(const std::vector<variable>& x) const
    {
        assert(x.size()==dim());
        double s=0;
        for(unsigned int i=0; i<fc.size(); i++)
            s += fc[i] * x[i];
        return s;
    }
private:
    std::vector<double> fc;
};

///@}



} // namespace

/// \defgroup sms Solution Methods


#endif // COMMONS_H
