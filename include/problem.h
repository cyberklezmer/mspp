#ifndef PROBLEM_H
#define PROBLEM_H

#include "commons.h"
#include "probability.h"
#include <limits>
#include <string>
#include <sstream>

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

const double inf = std::numeric_limits<double>::infinity();
const double minf = -inf;

struct varinfo
{
public:
    enum type { R, Rplus, Rminus };
    double l;
    double h;

    varinfo(double al=minf, double ah=inf, const char* an=0)
      : l(al),h(ah),n(an ? an : "")
    {}
    varinfo(type at,const char* an=0) : n(an ? an : "")
    {
        switch(at)
        {
        case R: l=minf; h=inf; break;
        case Rplus: l=0; h=inf; break;
        case Rminus: l=minf; h=0; break;
        }
    }
    std::string n;
};

using varinfo_list = std::vector<varinfo>;
using varinfo_list_ptr = std::shared_ptr<varinfo_list>;

template<typename C>
using constraint_list=std::vector<C>;

template<typename C>
using constraint_list_ptr = std::shared_ptr<constraint_list<C>>;

template<typename O>
using objective_ptr = std::shared_ptr<O>;

template<class O, class C, class Xi>
class problem : public object
{
    public:
        problem(const std::vector<unsigned int>& stagedims):fstagedims(stagedims)
        {
            assert(stagedims.size());
        }

        unsigned int dimupto(unsigned int k) const
        {
            assert(k < fstagedims.size());
            unsigned int s=0;
            for(int i=0; i<=k; i++)
                s+=fstagedims[i];
            return s;
        }

        unsigned int stageoffset(unsigned int k) const
        {
            assert(k < fstagedims.size());
            return k==0 ? 0 : dimupto(k-1);
        }

        unsigned int T() const { return fstagedims.size()-1; }

        unsigned int totaldim() const { return dimupto(T()); }

        unsigned int stagedim(unsigned int k) const
        {
            assert(k<fstagedims.size());
            return fstagedims[k];
        }
        const std::vector<unsigned int>& stagedims() const
        {
            return fstagedims;
        }

        virtual std::string varname(unsigned int k)
        {
            assert(k<totaldim());
            unsigned int s;
            for(s=0; s<T();s++)
            {
                if(stageoffset(s+1) > k)
                    break;
            }
            std::ostringstream st;
            st << varname(s,k-stageoffset(s)) << "_" << s;
            return st.str();
        }

        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            std::ostringstream s;
            s << "x" << i;
            return s.str();
        }

protected:
        virtual void f(
                unsigned int stage,
                const scenario<Xi>& xi,
                objective_ptr<O>& f) const = 0;

        virtual void constraints(
                unsigned int stage,
                const scenario<Xi>& xi,
                varinfo_list_ptr& xs,
                constraint_list_ptr<C>& constraints
                ) const = 0;
public:
        void get_f(
                unsigned int stage,
                const scenario<Xi>& xi,
                objective_ptr<O>& f) const
        {
            assert(stage <= T());
            assert(xi.size()==stage+1);
            this->f(stage,xi,f);
            assert(f.get());
        }

        void get_constraints(
                unsigned int stage,
                const scenario<Xi>& xi,
                varinfo_list_ptr& vars,
                constraint_list_ptr<C>& constraints
                ) const
        {
            assert(stage <= T());
            assert(xi.size()==stage+1);
            this->constraints(stage,xi,vars,constraints);
            assert(vars.get());
            assert(vars->size()==this->fstagedims[stage]);
            if(constraints) // tbd presunout  do linearproblemu., udelat precheck a postcheck
            {
                for(unsigned int i=0; i<constraints->size(); i++)
                    assert((*constraints)[i].lhs.size()==this->dimupto(stage));
            }
        }
    private:
        std::vector<unsigned int> fstagedims;
};


#endif // PROBLEM_H
