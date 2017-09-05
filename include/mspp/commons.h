#ifndef COMMONS_H
#define COMMONS_H

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <assert.h>
#include <limits>

// tbd
// - shiftedscenario to solver

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");


namespace mspp
{

const double inf = std::numeric_limits<double>::infinity();
const double minf = -inf;


class object
{
public:
    virtual ~object() {};
};

// tbd use this

class exception : public std::exception
{

public:
    exception(std::string error = "unknown mspp error")
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
    sys() : flog(0) { fself=this; }
    ~sys() { delete flog; fself = 0; }
    static void set_log(const std::string& s)
    {
        delete fself->flog;
        fself->flog = new std::ofstream(s);
        if(!(fself->flog ))
        {
            std::cerr << "Cannot initialize log file " << s << std::endl;
            throw;
        }
    }
    static void reset_log()
    {
        delete fself->flog;
        fself->flog = 0;
    }

    static std::ostream& log()
    {  return fself->flog ? *(fself->flog) : std::cout; }
private:
    std::ostream* flog;
};


using path = std::vector<unsigned int>;

class tree;

class treecallback
{
public:
    virtual void callback(const path& p) = 0;
};

class tree : public object
{
public:
    virtual unsigned int depth() const = 0;
    virtual unsigned int numbranches(const path& p) const = 0;
    void foreachnode(treecallback *callee, path p=path(0))
    {
         if(depth() > p.size())
             doforeachnode(callee,p);
    }
private:
    void doforeachnode(treecallback *callee, path& ap)
    {
        unsigned int k=ap.size();
        unsigned int n=numbranches(ap);

        path p(ap);
        p.push_back(0);

        bool cb = k+1 < depth();
        for(p[k]=0; p[k]<n; p[k]++)
        {
            callee->callback(p);
            if(cb)
                doforeachnode(callee, p);
        }
    }
};

using tree_ptr = std::shared_ptr<tree>;



template <typename L>
using branch = std::vector<L>;

template <typename L>
class treemapping : public object
{
public:
    treemapping(const tree_ptr& t) : ft(t) {}

    // tbd possibly make a reference
    virtual L operator()(const path& p) const = 0;
    virtual branch<L> s(const path& p) const
    {
        branch<L> r;
        path np;
        for(unsigned int i=0; i<p.size(); i++)
        {
            np.push_back(p[i]);
            r.push_back((*this)(np));
        }
        return r;
    }
    const tree_ptr& t() const { return ft; }
protected:
    tree_ptr ft;
};

template <typename L>
using treemapping_ptr = std::shared_ptr<treemapping<L>>;


template <typename Xi>
using process=treemapping<Xi>;

template <typename Xi>
using scenario = branch<Xi>;


template <typename Xi>
using process_ptr = std::shared_ptr<process<Xi>>;


using prob = double;

class probability: public treemapping<prob>
{
public:
    probability(const tree_ptr& t) : treemapping(t) {}
    virtual prob up(const path& ap) const
    {
        prob pr=1.0;
        path p;
        for(unsigned int i=0; i<ap.size(); i++)
        {
            p.push_back(ap[i]);
            pr*=(*this)(p);
        }
        return pr;
    }
};

using probability_ptr = std::shared_ptr<probability>;




template <typename Xi>
class scenariotree : public object
{
public:
    scenariotree()
      {}

    virtual const process_ptr<Xi>& tm() const = 0;
    virtual const probability_ptr& tp() const = 0;

    const tree_ptr& t() const { return tp()->t(); }
    virtual unsigned int depth() const { return tp()->t()->depth();}
    Xi x(const path& p) const  { return (*tm())(p); }
    scenario<Xi> s(const path& p) const { return tm()->s(p); }
    virtual prob p(const path& p) const { return (*tp())(p); }
    virtual prob up(const path& p) const { return tp()->up(p); }
};

template<class Xi>
using scenariotree_ptr=std::shared_ptr<scenariotree<Xi>>;



template<class Xi>
class modularscenariotree : public scenariotree<Xi>
{
public:
    modularscenariotree(const process_ptr<Xi>& x, const probability_ptr& p)
         : fx(x), fp(p)
    {
    }

    virtual const process_ptr<Xi>& tm() const { return fx; }
    virtual const probability_ptr& tp() const { return fp; }

protected:
    process_ptr<Xi> fx;
    probability_ptr fp;
};


template<class Xi>
using modularscenariotree_ptr=std::shared_ptr<modularscenariotree<Xi>>;


struct varinfo
{
public:
    enum type { R, Rplus, Rminus };

    varinfo(double al=minf, double ah=inf)
      : l(al),h(ah)
    {}

    varinfo(type at)
    {
        switch(at)
        {
        case R: l=minf; h=inf; break;
        case Rplus: l=0; h=inf; break;
        case Rminus: l=minf; h=0; break;
        }
    }

    double l;
    double h;
};

using varinfo_list = std::vector<varinfo>;
using varinfo_list_ptr = std::shared_ptr<varinfo_list>;

template<typename C>
using constraint_list=std::vector<C>;

template<typename C>
using constraint_list_ptr = std::shared_ptr<constraint_list<C>>;

template<typename O>
using objective_ptr = std::shared_ptr<O>;

/// Base class for MSSPP implementations.
/**
 * Even though it is not abstract, it itself cannot be used
 * (if used, it raises an assertion error).
*/
template<class O, class C, class Xi>
class problem : public object
{
    public:
        problem(const std::vector<unsigned int>& stagedims)
            : fstagedims(stagedims)
        {
            assert(stagedims.size());
        }

        unsigned int T() const { return fstagedims.size()-1; }

        const std::vector<unsigned int>& stagedims() const
        {
            return fstagedims;
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

        unsigned int stagedim(unsigned int k) const
        {
            assert(k<fstagedims.size());
            return fstagedims[k];
        }

        unsigned int totaldim() const { return dimupto(T()); }


        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            std::ostringstream s;
            s << "x" << i;
            return s.str();
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

protected:
        /**
         * @brief f  should return the objective function
         * @param k  the stage
         * @param barxi  history of the random vector up to #k
         * @param f  return value
         */
        virtual void f(unsigned int k,const scenario<Xi>& barxi,O& f) const
        {
            assert(0);
        }
        /**
         * @brief constraints should return
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param xs return value - variable ranges
         * @param c return value - constraints
         */
        virtual void constraints(
                unsigned int k,
                const scenario<Xi>& barxi,
                varinfo_list& xs,
                constraint_list<C>& c
                ) const
        {
            assert(0);
        }

public:
        virtual void f(
                unsigned int stage,
                const scenario<Xi>& xi,
                objective_ptr<O>& af) const
        {
            af.reset(new O);
            f(stage,xi,*af);
        }

        virtual void constraints(
                unsigned int stage,
                const scenario<Xi>& xi,
                varinfo_list_ptr& xs,
                constraint_list_ptr<C>& csts
                ) const
        {
            xs.reset(new varinfo_list);
            csts.reset(new constraint_list<C>);
            this->constraints(stage,xi,*xs,*csts);
            assert(xs->size()==this->fstagedims[stage]);
        }

public:
    private:
        std::vector<unsigned int> fstagedims;
};

class treesolution: public object, public treecallback
{
    enum foreachnodemode {elist, estats};
public:
    treesolution(const std::vector<unsigned int>& stagedims,
                 const probability_ptr& tp,
                 std::vector<double>& x) :
             fsd(stagedims), fso(stagedims.size(),0), ftp(tp), fx(x)
    {
        unsigned int o=0;
        for(unsigned int i=0; i<fsd.size(); i++)
        {
            fso[i] = o;
            o += fsd[i];
        }
    }
    void list(std::ostream& os)
    {
        fimode = elist;
        fos = &os;
        fi = 0;
        ftp->t()->foreachnode(this);
        assert(fi==fx.size());

    }

    void stats(std::vector<double>& E, std::vector<double>& var)
    {
        unsigned int totaldim = 0;
        for(unsigned int i=0; i<fsd.size(); i++)
            totaldim += fsd[i];

        fs.resize(totaldim);
        std::fill(fs.begin(), fs.end(), 0);

        fs2.resize(totaldim);
        std::fill(fs2.begin(), fs2.end(), 0);

        fi=0;

        fimode = estats;
        ftp->t()->foreachnode(this);

        assert(fi==fx.size());

        E.resize(totaldim);
        var.resize(totaldim);
        for(unsigned int i=0; i<totaldim; i++)
        {
            E[i] = fs[i];
            var[i] = fs2[i] - E[i]*E[i];
        }
    }
    double x(unsigned int i) { return fx[i]; }
    double nx()
    {
        return fx.size();
    }
    const probability_ptr& tp() { return ftp; }

    unsigned int sd(unsigned int i) { return fsd[i];}
    unsigned int so(unsigned int i) { return fso[i];}
    unsigned int numstages() { return fsd.size();}

protected:

    std::vector<unsigned int> fsd;
    std::vector<unsigned int> fso;
    probability_ptr ftp;
    std::vector<double> fx;

// state variables

protected:
    unsigned int fi;
private:
    foreachnodemode fimode;
    std::ostream* fos;
    std::vector<double> fs;
    std::vector<double> fs2;

protected:
    virtual void callback(const path& p)
    {
        prob pr = ftp->up(p);
        unsigned int s = p.size()-1;
        if(fimode==elist)
        {
            for(unsigned int i=0; ; i++)
            {
                (*fos) << p[i];
                if(i== p.size()-1)
                    break;
                *fos << "-";
            }

            for(unsigned int i=0; i<s; i++)
                for(unsigned int j=0; j<fsd[i]; j++)
                    (*fos) << ",";
            for(unsigned int i=0; i<fsd[s] ; i++)
            {
                *fos << "," << fx[fi++];
            }
            *fos << std::endl;
        }
        else if(fimode==estats)
        {
            for(unsigned int i=0; i<fsd[s]; i++)
            {
                unsigned int k = fso[s]+i;
                fs[k] += pr * fx[fi];
                fs2[k] += pr * fx[fi]*fx[fi];
                fi++;
            }
        }
        assert(fi <= fx.size());
    }
};

using treesolution_ptr = std::shared_ptr<treesolution>;


template<typename P, typename D, typename S>
class solutionmethod : public object
{
public:
    solutionmethod(const std::shared_ptr<P>& problem,
                   const std::shared_ptr<D>& distribution) :
       fp(problem), fd(distribution)
    {

    }
    virtual void solve(std::shared_ptr<S>& sol, double& optimal) = 0;
protected:
    std::shared_ptr<P> fp;
    std::shared_ptr<D> fd;
};


}


#endif // COMMONS_H
