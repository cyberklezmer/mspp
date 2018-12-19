#ifndef MARKOV_H
#define MARKOV_H

#include "mspp/random.h"

namespace mspp
{


/// \addtogroup Distributions
/// @{

/// \addtogroup markov Markov distributions
/// \ingroup Distributions
/// @{

class mcdistribution : virtual public ddistribution<unsigned int,unsigned int>
{
public:
    unsigned int nstates() const
    {
        return nstates_is();
    }
private:
    virtual probability transprob(unsigned int from, unsigned int to) const = 0;
    virtual unsigned int nstates_is() const = 0;
protected:
    virtual atom<unsigned int>
         atom_is(unsigned int i,const unsigned int& c) const
    {
        return {i,this->transprob(c,i)};
    }
};

class fmcdistribution :
        public mcdistribution,
        public fdistribution<unsigned int, unsigned int>
{
public:
    unsigned int nstates() const
    {
        return nstates_is();
    }
private:
    virtual unsigned int natoms_is(const unsigned int& ) const
    {
         return nstates();
    }
    virtual atom<unsigned int> atom_is(unsigned int i,const unsigned int& c) const
    {
        return mcdistribution::atom_is(i,c);
    }

    virtual unsigned int nstates_is() const = 0;
};

class mmcdistribution : public fmcdistribution
{
public:
    mmcdistribution(const vector<vector<double>> m) : fm(m)
    {
        assert(fm.size());
        assert(fm[0].size());
        for(unsigned int i=0; i<fm.size(); i++)
            for(unsigned int j=1; j<fm.size(); j++)
                assert(fm[i].size()==fm[0].size());
    }
private:
    virtual probability transprob(unsigned int from, unsigned int to) const
    {
        assert(from < fm.size());
        assert(to<fm[0].size());
        return fm[from][to];
    }
    virtual unsigned int nstates_is() const
    {
        return fm[0].size();
    }
    vector<vector<double>> fm;
};

/// @}

/// @}


} // namespace

#endif // MARKOV_H
