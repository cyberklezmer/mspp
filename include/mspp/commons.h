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
    static sys* fself;
public:
    sys() { fself=this; }
    ~sys() { fself = 0; }
    static void set_log(const ptr<std::ostream>& l)
    {
        fself->flog = l;
    }
    static void reset_log()
    {
        fself->flog.reset();
    }

    static void set_err(const ptr<std::ostream>& e)
    {
        fself->ferr = e;
    }
    static void reset_err()
    {
        fself->ferr.reset();
    }

    static std::ostream& log()
    {  return fself->flog ? *(fself->flog) : std::cout; }
    static std::ostream& err()
    {  return fself->ferr ? *(fself->ferr) : std::cerr; }
private:
    ptr<std::ostream> flog;
    ptr<std::ostream> ferr;
};

/// @}



/// \addtogroup general General Definitions
///  @{


template <typename X>
using scenario = std::vector<X>;

template <typename X>
class condition: public object
{
};

/// \addtogroup cns Conditions
/// \ingroup general
/// @{



template <typename X>
class nocondition: public condition<X>
{
public:
    nocondition(const scenario<X>&) {}
};

template <typename X>
class fullpath: public condition<X>
{
public:
    fullpath(const scenario<X>& s): fs(s)
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

/// @}

/// @}

} // namespace


#endif // COMMONS_H
