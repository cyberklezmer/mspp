#ifndef SOLUTION_H
#define SOLUTION_H

#include "probability.h"


class treesolution: public object, public treecallback
{
public:
    treesolution(const std::vector<unsigned int>& stagedims,
                 const treeprobability_ptr& tp,
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

private:
    std::vector<unsigned int> fsd;
    std::vector<unsigned int> fso;
    treeprobability_ptr ftp;
    std::vector<double> fx;

// state variables

    unsigned int fi;
    std::vector<double> fs;
    std::vector<double> fs2;

    virtual void callback(const path& p)
    {
        prob pr = ftp->up(p);
        unsigned int s = p.size()-1;
        for(unsigned int i=0; i<fsd[s]; i++)
        {
            unsigned int k = fso[s]+i;
            fs[k] += pr * fx[fi];
            fs2[k] += pr * fx[fi]*fx[fi];
            fi++;
        }
        assert(fi <= fx.size());
    }
};

using treesolution_ptr = std::shared_ptr<treesolution>;

#endif // SOLUTION_H
