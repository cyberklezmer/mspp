#ifndef TESTS_H
#define TESTS_H

#include "mspp/linear.h"
//#include "mspp/mmpcvar.h"
#include "mspp/lpsolver.h"
//#include "io.h"
#include <fstream>

namespace mspp
{

class csvlpsolver : public lpsolver
{
public:

    virtual void solve(const std::vector<realvar>& vars,
            const std::vector<double>& f,
            const std::vector<sparselinearconstraint_ptr>& constraints,
            const std::vector<std::string>& varnames,
            std::vector<double>& x,
            double& ) const
    {
        using namespace std;
        ofstream os("problem.csv");
        os << "variables" << endl;
        for(int i=0;i<vars.size(); i++)
            os << varnames[i] << ",";
        os << endl;
        os << endl;

        os << "f function" << endl;
        for(int i=0;i<f.size(); i++)
           os << f[i] << ",";
        os << endl;

        os << "lower consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].l() == minf)
                os << ",";
            else
                os << vars[i].l()<< ",";
        os << endl;

        os << "upper consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].h() == inf)
                os << ",";
            else
                os << vars[i].h() << ",";

        os << endl;

        for(int y=0; y<3; y++)
        {
            linearconstraint::type t = (linearconstraint::type) y;
            os << "consttraints ";
            switch(t)
            {
                case linearconstraint::eq:
                    os << "=";
                break;
                case linearconstraint::geq:
                    os << ">=";
                break;
                case linearconstraint::leq:
                    os << "<=";
                break;
            }
            os << endl;
            for(int i=0; i<constraints.size(); i++)
            {
                const sparselinearconstraint& c = *constraints[i];
                if(c.t == t)
                {
/*                    shared_ptr<std::vector<double>> clhs(c.lhs());
                    int j=0;
                    for(; j<clhs->size(); j++)
                        os << (*clhs)[j] << ",";*/

                    int j=0;
                    for(; j<c.lhs()->size(); j++)
                        os << (*c.lhs())[j] << ",";
                    for(; j<vars.size(); j++)
                        os << "0,";
                    os << c.rhs << endl;
                }
             }
        }
        std::cout << "Test solver only produces problem.csv osile.";
    }
};

/*

template<typename Xi>
class scenariolister : public treecallback
{
public:
    scenariolister(const scenariotree_ptr<Xi>& s, std::ostream& o)
        : fs(s), fo(o)
    {}
    void callback(const path& p)
    {
        for(unsigned int i=0; ; i++)
        {
            fo << p[i];
            if(i==p.size()-1)
                break;
            fo <<"-";
        }
        fo << ",p=" << fs->p(p) ;
        scenario<Xi> sc = fs->s(p);
        for(unsigned int i=0; i<p.size(); i++)
            fo << "," << sc[i];
        fo << std::endl;
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
        if(++i==p.sumd())
            break;
        std::cout << ",";
    }
    std::cout << std::endl;
}
*/

}
#endif // TESTS_H
