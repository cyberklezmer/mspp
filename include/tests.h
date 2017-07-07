#ifndef TESTS_H
#define TESTS_H

#include "linear.h"
#include "scenarios.h"
#include "mmpcvar.h"
#include <fstream>

class csvlpsolver : public lpsolver
{
public:

    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<sparselinearconstraint_ptr>& constraints,
            const std::vector<std::string>& varnames,
            std::vector<double>& x,
            double& ) const
    {
        using namespace std;
        ofstream f("problem.csv");
        f << "variables" << endl;
        for(int i=0;i<vars.size(); i++)
            f << varnames[i] << ",";
        f << endl;
        f << endl;

        f << "objective function" << endl;
        for(int i=0;i<objective.coefs.size(); i++)
           f << objective.coefs[i] << ",";
        f << endl;

        f << "lower consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].l == minf)
                f << ",";
            else
                f << vars[i].l << ",";
        f << endl;

        f << "upper consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].h == inf)
                f << ",";
            else
                f << vars[i].h << ",";

        f << endl;

        for(int y=0; y<3; y++)
        {
            linearconstraint::type t = (linearconstraint::type) y;
            f << "consttraints ";
            switch(t)
            {
                case linearconstraint::eq:
                    f << "=";
                break;
                case linearconstraint::geq:
                    f << ">=";
                break;
                case linearconstraint::leq:
                    f << "<=";
                break;
            }
            f << endl;
            for(int i=0; i<constraints.size(); i++)
            {
                const sparselinearconstraint& c = *constraints[i];
                if(c.t == t)
                {
                    shared_ptr<std::vector<double>> clhs(c.lhs());
                    int j=0;
                    for(; j<clhs->size(); j++)
                        f << (*clhs)[j] << ",";
                    for(; j<vars.size(); j++)
                        f << "0,";
                    f << c.rhs << endl;
                }
             }
        }
        std::cout << "Test solver only produces problem.csv file.";
    }
};


template<typename Xi>
class scenariolister : public treecallback
{
public:
    scenariolister(const scenariotree_ptr<Xi>& s, std::ostream& o)
        : fs(s), fo(o)
    {}
    void callback(const path& p)
    {
        for(unsigned int i=0; i<p.size(); i++)
        {
            fo << p[i] <<"-";
        }
        fo << "p=" << fs->p(p) << " ";
        fo << fs->x(p) << std::endl;
    }
    void list()
    {
        fs->t()->foreachnode(this);
    }

private:
    scenariotree_ptr<Xi> fs;
    std::ostream& fo;
};

template<typename T>
void output(const std::vector<T>& v, std::ostream& os)
{
    for(unsigned int i=0;;)
    {
        os << v[i];
        if(++i==v.size())
            break;
        os << ",";
    }
}


inline void printstats(treesolution& ts)
{
    std::vector<double> E,var;
    ts.stats(E,var);
    assert(E.size());
    assert(var.size());

    std::cout << "E=(";

    output(E,std::cout);
    std::cout << ")" << std::endl;

    std::cout << "var=(";

    output(var,std::cout);

    std::cout << ")" << std::endl;
}

template<typename O, typename C, typename Xi>
inline void printvarnames(problem<O,C,Xi>& p)
{
    for(unsigned int i=0;;)
    {
        std::cout << p.varname(i);
        if(++i==p.totaldim())
            break;
        std::cout << ",";
    }
    std::cout << std::endl;
}

#endif // TESTS_H
