#ifndef COMMONS_H
#define COMMONS_H

#include <vector>
#include <memory>
#include <assert.h>

// stages index from zero to T-1


class object
{
public:
    virtual ~object() {};
};



/* TBD
template<class X>
class distribution
{
    public:
        distribution(unsigned int adim): dim(adim)
        {}

        distribution(const distribution& adist): dim(adist.dim) {}
        virtual ~distribution() = 0;

    protected:

    private:
        unsigned int dim;
};

class solver
{
    static void checkstage(tbd );
};
*/

typedef std::vector<double> stagesolution;

#endif // COMMONS_H
