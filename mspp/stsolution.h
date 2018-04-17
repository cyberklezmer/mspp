#ifndef STSOLUTION_H
#define STSOLUTION_H

#include "mspp/msproblem.h"
#include "mspp/distribution.h"

namespace mspp
{

template <typename P, typename S>
class stsolution: public object, public sctreecallback<variables>
{
public:
    stsolution(const P& p, const S& s) :
        fs(ptr<S>(new S(s))), fps(p.d)
       {}
    stsolution(const P& p, const ptr<S> s) :
        fs(s), fps(p.d)
       {}
    void set(const variables& x) { fx=ptr<variables>(new variables(x)); }
    void set(const ptr<variables> x) { fx=x; }
    void foreach()
    {
        fs->foreach(this);
    }
    virtual void callback(const indexedhistory<variables>&) {};
    variable &x(unsigned int i) const { return (*fx)[i]; }
private:
    ptr<variables> fx;
    ptr<S> fs;
    msproblemstructure fps;
};

} // namespace


#endif // STSOLUTION_H
