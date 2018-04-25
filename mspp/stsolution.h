#ifndef STSOLUTION_H
#define STSOLUTION_H

#include "mspp/msproblem.h"
#include "mspp/distribution.h"

namespace mspp
{

template<typename X>
struct stsolitem { X x; variables v; };

template<typename X>
class stsolutioncallback
{
public:
    virtual void callback(const indexedhistory<stsolitem<X>>& s) = 0;
};

template <typename P, typename S>
class stsolution: public object, public sctreecallback<typename S::X_t>
{
public:
    using P_t = P;
    using S_t = S;
    using X_t = typename S::X_t;

    stsolution(const P& p, const S& s) :
        fs(ptr<S>(new S(s))), fps(p.d), fcallee(0)
       {}
    stsolution(const P& p, const ptr<S> s) :
        fs(s), fps(p.d), fcallee(0)
       {}
    void set(const variables& x) { fx=ptr<variables>(new variables(x)); }
    void set(const ptr<variables> x) { fx=x; }
    const msproblemstructure& ps() const { return fps; }
private: // state variables
    stsolutioncallback<typename S::X_t> *fcallee;
    indexedhistory<stsolitem<typename S::X_t>> fh;
    unsigned int fi;
public:
    void foreachnode(stsolutioncallback<typename S::X_t> *callee)
    {
        assert(fcallee == 0);
        assert(fh.size()==0);
        fcallee = callee;
        fi=0;
        fs->foreachnode(this);
        assert(fi==fx->size());
        fh.clear();
        fcallee = 0;
    }
    virtual void callback(const indexedhistory<typename S::X_t>& h)
    {
        assert(h.size() <= fps.size());
        assert(h.size());
        unsigned int k = h.size()-1;
        variables n;
        for(unsigned int i=0; i<fps[k]; i++)
        {
            assert(fi<fx->size());
            n.push_back((*fx)[fi++]);
        }
        fh.resize(h.size());
        const indexedatom<typename S::X_t>& a = h[k];
        fh[k] = { {a.x,n }, a.p, a.i };
        for(unsigned int i=0; i<h.size(); i++)
            assert(fh[i].i==h[i].i);
        fcallee->callback(fh);
    }
//    virtual void callback(const indexedhistory<stsolitem<typename S::X_t>>&) {};
    variable &x(unsigned int i) const { return (*fx)[i]; }
    const variables x() const { return *fx; }
private:
    ptr<variables> fx;
    ptr<S> fs;
    msproblemstructure fps;
};

template <typename S, typename D>
class stsolreducer : public object, public stsolutioncallback<typename S::X_t>
{
public:
    void convert(S& s, D& d)
    {
        static_assert(std::is_same<typename S::S_t, typename D::S_t>::value);
        assert(!fv);
        fv.reset(new variables);
        fdps = d.ps();
        s.foreachnode(this);
        d.set(fv);
        fv.reset();
    }

private:
// state variables
    ptr<variables> fv;
    msproblemstructure fdps;
    virtual void callback(const indexedhistory<stsolitem<typename S::X_t>>& h)
    {
        assert(h.size());
        assert(h.size() <= fdps.size());
        unsigned int k=h.size()-1;
        const variables& v=h[k].x.v;
        assert(v.size() >= fdps[k]);
        for(unsigned int i=0; i<fdps[k]; i++)
            fv->push_back(v[i]);

    }
};


} // namespace


#endif // STSOLUTION_H
