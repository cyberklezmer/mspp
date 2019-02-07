#ifndef CPLEX_H
#define CPLEX_H

#include "mspp/lpsolver.h"
#include <ilcplex/ilocplex.h>

namespace mspp
{

/// \addtogroup solvers Solvers
/// @{



template <typename V>
class cplex : public lpsolver<V>
{
public:
    virtual void solve(
            const lpproblem<V>& lp,
            vector<double>& sol,
            double& objvalue) const
    {
        assert(lp.check());

        IloEnv env;

        try
        {
            IloModel model(env);

            IloNumVarArray x(env);

            for(int i=0; i<lp.vars.size(); i++)
            {
                const range<V>& v = lp.vars[i];

                const IloNum l = v.islinf() ? -IloInfinity : v.l();
                const IloNum h = v.ishinf() ? IloInfinity : v.h();
                IloNumVar::Type t;
                switch(v.t())
                {
                    case range<V>::realt:
                       t = IloNumVar::Float;
                    break;
                    case range<V>::intt:
                       t = IloNumVar::Int;
                    break;
                    case range<V>::bint:
                       t = IloNumVar::Bool;
                    break;
                default: assert(0);
                }

                x.add(IloNumVar(env,l,h,t,lp.varnames[i].c_str()));
            }

            IloObjective obj = IloMinimize(env);

            for(int i=0; i<lp.f.size(); i++)
                obj.setLinearCoef(x[i], lp.f[i]);

            model.add(obj);

            for(int k=0; k< lp.constraints.size();
                k++)
            {
                IloExpr v(env);

                const sparselinearconstraint& c = lp.constraints[k];

    //            assert(c.lhs.size()==vars.size());

                for(int i=0; i<c.numnz(); i++)
                {
                    const sparselinearconstraint::iitem& it = c.nz(i);
                    if(it.value) // tbd ensure in lc
                        v+= it.value*x[it.index];
                }
                switch(c.t())
                {
                case constraint::eq:
                    model.add(v==c.rhs);
                    break;
                case constraint::leq:
                    model.add(v<=c.rhs);
                    break;
                case constraint::geq:
                    model.add(v>=c.rhs);
                    break;
                default:
                    assert(1);
                }
                v.end();
            }

            IloCplex cplex(model);
            cplex.setOut(sys::log());
            cplex.setWarning(sys::log());
            cplex.exportModel("model.lp"); //lp, mps, sav


           // Optimize the problem and obtain solution.
           int serr;
           if ( !(serr = cplex.solve()) )
           {
              sys::err() << "Error " << serr << " when optimizing LP" << std::endl;
              throw exception("Error  when optimizing LP", serr);
           }

           for(int i=0; i<lp.vars.size(); i++)
               sol[i] = cplex.getValue(x[i]);
           objvalue = cplex.getObjValue();
//std::cout << std::endl << "ov:" << cplex.getObjValue() <<std::endl<< std::endl;
        }
        catch (IloException& e)
        {
           throw mspp::exception(e.getMessage());
        }
        catch (...)
        {
           throw mspp::exception("Unknown concert exception caught" );
        }

        env.end();
    }

};

/// @}

} // namespace


#endif // CPLEX_H
