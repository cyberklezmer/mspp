#ifndef STSOLUTION_H
#define STSOLUTION_H

#include "mspp/msproblem.h"
#include "mspp/random.h"

namespace mspp
{

template<typename P, typename Z>
struct stsolitem
{
    rvector<typename Z::X_t> xi;
    variables<typename P::V_t> x;
};


template <typename P, typename Z>
class stsolcallback
{
public:
    virtual void callback(const indexedpath<stsolitem<P,Z>>& path, void* state=0) const;
};

template<typename P, typename Z>
struct stsolutionstate
{
public:
    stsolcallback<P,Z> *callee;
    indexedpath<stsolitem<P,Z>> h;
    unsigned int i;
    void* callees;
};


template <typename P, typename Z>
class stsolution: public object,
        public sctreecallback<typename Z::X_t, stsolutionstate<P,Z>>
{
public:
    using P_t = P;
    using Z_t = Z;
    using X_t = typename Z::X_t;
    using V_t = typename P::V_t;

    stsolution(const P& p, const Z& z) :
        fz(ptr<Z>(new Z(z))), fps(p.d)
       {}
    stsolution(const P& p, const ptr<Z> z) :
        fz(z), fps(p.d)
       {}
    void set(const variables<V_t>& x) { fx=ptr<variables<V_t>>(new variables<V_t>(x)); }
    void set(const ptr<variables<V_t>> x) { fx=x; }
    const msproblemstructure& ps() const { return fps; }

    void foreachnode(stsolcallback<P,Z> *callee, void* cs=0) const
    {
        stsolutionstate<P,Z> s;
        s.callee = callee;
        s.i = 0;
        s.callees = cs;
        fz->foreachnode(this, &s);
    }


    virtual void callback(const indexedpath<rvector<X_t>>& h,
                                 stsolutionstate<P,Z>* s) const
    {
        assert(h.size() <= fps.size());
        assert(h.size());
        unsigned int k = h.size()-1;
        variables<V_t> n;
        for(unsigned int i=0; i<fps[k]; i++)
        {
            assert(s->i<fx->size());
            n.push_back((*fx)[s->i++]);
        }
        s->h.resize(h.size());
        auto a=h[k];
        s->h[k] = { stsolitem<P_t,Z_t>({a.x,n }), a.p,  a.i };

        for(unsigned int i=0; i<h.size(); i++)
            assert(s->h[i].i==h[i].i);
        s->callee->callback(s->h, s->callees);
    }
    V_t &x(unsigned int i) const { return (*fx)[i]; }
    const variables<V_t> x() const { return *fx; }
private:
    ptr<variables<V_t>> fx;
    ptr<Z> fz;
    msproblemstructure fps;
};

template <typename V>
struct  stsolreducerstate
{
    ptr<variables<V>> fv;
    msproblemstructure fdps;
};

template <typename S, typename D>
class stsolreducer : public object,
        public stsolcallback
          <typename S::P_t, typename S::Z_t>
{
public:
    using V_t = typename S::P_t::V_t;
    void convert(S& s, D& d) //const
    {
        static_assert(std::is_same<typename S::Z_t, typename D::Z_t>::value);
        stsolreducerstate<V_t> st;
        st.fv.reset(new variables<V_t>);
        st.fdps = d.ps();
        s.foreachnode(this, &st);
        d.set(*(st.fv));
    }

private:
// state variables
    virtual void callback
      (const indexedpath<stsolitem<typename S::P_t, typename S::Z_t>>& h,
                           void* sptr) const
    {
        assert(h.size());
        stsolreducerstate<V_t>* state
                = static_cast<stsolreducerstate<V_t>*>(sptr);
        assert(h.size() <= state->fdps.size());
        unsigned int k=h.size()-1;
        const variables<V_t>& v=h[k].x.x;
        assert(v.size() >= state->fdps[k]);

        for(unsigned int i=v.size()-state->fdps[k]; i<v.size(); i++)
            state->fv->push_back(v[i]);
    }
};



} // namespace


#endif // STSOLUTION_H
