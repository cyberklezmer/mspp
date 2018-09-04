#ifndef STSOLUTION_H
#define STSOLUTION_H

#include "mspp/msproblem.h"
#include "mspp/process.h"

namespace mspp
{

template<typename P, typename X>
struct stsolitem
{
    typename X::I_t xi;
    variables<typename P::V_t> x;
};


template <typename P, typename X>
class stsolcallback
{
public:
    virtual void callback(const indexedpath<stsolitem<P,X>>& path, void* state=0) const;
};

template<typename P, typename X>
struct stsolutionstate
{
public:
    stsolcallback<P,X> *callee;
    indexedpath<stsolitem<P,X>> h;
    unsigned int i;
    void* callees;
};


template <typename P, typename X>
class stsolution: public object,
        public sctreecallback<typename X::I_t, stsolutionstate<P,X>>
{
public:
    using P_t = P;
    using X_t = X;
    using I_t = typename X::I_t;
    using V_t = typename P::V_t;

    stsolution(const P& p, const X& x) :
        fxi(ptr<X>(new X(x))), fps(p.d)
       {}
    stsolution(const P& p, const ptr<X> x) :
        fxi(x), fps(p.d)
       {}
    void set(const variables<V_t>& x)
    {
        fx=ptr<variables<V_t>>(new variables<V_t>(x));
    }
    void set(const ptr<variables<V_t>> x)
    {
        fx=x;
    }
    const msproblemstructure& ps() const { return fps; }

    void foreachnode(stsolcallback<P,X> *callee, void* cs=0) const
    {
        stsolutionstate<P,X> s;
        s.callee = callee;
        s.i = 0;
        s.callees = cs;
        fxi->foreachnode(this, &s);
    }


    virtual void callback(const indexedpath<I_t>& h,
                                 stsolutionstate<P,X>* s) const
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
        s->h[k] = { stsolitem<P_t,X_t>({a.x,n }), a.p,  a.i };

        for(unsigned int i=0; i<h.size(); i++)
            assert(s->h[i].i==h[i].i);
        s->callee->callback(s->h, s->callees);
    }
    V_t &x(unsigned int i) const { return (*fx)[i]; }
    const variables<V_t> x() const { return *fx; }
private:
    ptr<variables<V_t>> fx;
    ptr<X> fxi;
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
          <typename S::P_t, typename S::X_t>
{
public:
    using V_t = typename S::P_t::V_t;
    void convert(S& s, D& d) //const
    {
        static_assert(std::is_same<typename S::X_t, typename D::X_t>::value);
        stsolreducerstate<V_t> st;
        st.fv.reset(new variables<V_t>);
        st.fdps = d.ps();
        s.foreachnode(this, &st);
        d.set(*(st.fv));
    }

private:
// state variables
    virtual void callback
      (const indexedpath<stsolitem<typename S::P_t, typename S::X_t>>& h,
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
