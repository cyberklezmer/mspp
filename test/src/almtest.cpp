/*
class almproblem: public linearproblem<double>
{

public:
    almproblem() : linearproblem<double>({1,1,1,1})
    {}
    virtual void f(unsigned int stage,
                   const scenario<double>& xi,
                   linearobjective& r) const
    {
       double c=xi[0];
       for(unsigned int i=1; i<=stage; i++)
          c *= xi[i];
       r.coefs[0] = c;
    }

    virtual void msconstraints(
                unsigned int stage,
                const scenario<double>& xi,
                varranges& vars,
                constraint_list<linearconstraint>& csts) const
    {
        if(stage==0)
           vars[0]=varrange(0,1);
        else if(stage==1)
        {
           vars[0]=varrange(varrange::Rplus);
           csts.add(linearconstraint({1.0,1.0},linearconstraint::geq, 0.0));
           csts.add(linearconstraint({1.0,1.0},linearconstraint::leq, 1.0));
        }
        else if(stage==2)
        {
           vars[0]=varrange(varrange::Rplus);
           csts.add(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::geq, 0.0));
           csts.add(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::leq, 1.0));
        }
        else
        {
            vars[0]=varrange(varrange::Rplus);
            csts.add(linearconstraint({1.0,1.0,1.0,1.0},
                                     linearconstraint::eq, 1.0));
        }
    }
};

void almtest(double alpha, double lambda, const lpsolver_ptr& cps)
{
    const int numleaves = 2;
    ntree_ptr tp(new ntree(3,numleaves));
    iidprobability_xxx_ptr pp(new iidprobability_xxx(tp,{0.5,0.5}));
    iidprocess_ptr<double> xp(new iidprocess<double>(tp,{0.8,1.1}));

    shiftedscenariotree_ptr<double> ss(new shiftedscenariotree<double>(xp,pp,1.0));

    linearproblem_ptr<double> prp(new almproblem());

    std::cout <<"ALMtest..." << std::endl;

    linearproblem_ptr<double> cvp;
    cvp.reset(new mmpcvarequivalent<double>(prp,alpha,lambda));

    printvarnames(*cvp);
    biglpmethod<double> b(cvp,ss,cps);

    treesolution_ptr s;
    double ov;
    b.solve(s,ov);

    printstats(*s);

    const double tol = 1e-5;

    double almsol[]={0,0,0,0,0.64,1,0,0,0,0,0,0,0.727273,0,0.264,0.272727,-0.036,0.272727,0,0,
                     0,0.88,1,0,0,0,0,0,0,0.727273,0,0.363,0.272727,-0.0495,0.272727,0};
    double almobj = 0.906063;

    if(fabs(ov-almobj) > tol)
    {
        std::cerr << "almtest: opt="
             << almobj << " expected, " << ov << " achieved." << std::endl;
        throw;
    }

    unsigned int k = sizeof(almsol)/sizeof(almsol[0]);
    for(unsigned int i=0; i<k; i++)
    {
        if(fabs(s->x(i)-almsol[i]) > tol)
        {
            std::cerr << "almtest: x[" << i << "]="
                 <<https://slovnik.seznam.cz/en/?q=z%C3%A1stupce almsol[i] << " expected, " << s->x(i)
                 << " achieved." << std::endl;
            std::cerr << "x=";

            for(unsigned int j=0; j<s->nx(); j++)
            {
                std::cerr << s->x(j) << ",";
                if((j+1)%20 == 0)
                     std::cerr << std::endl;
            }
            throw;
        }
    }
    std::cout <<  " ALMtest passed." << std::endl;
}

*/
