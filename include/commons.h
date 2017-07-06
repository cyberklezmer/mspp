#ifndef COMMONS_H
#define COMMONS_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <assert.h>

// stages index from zero to T-1


class object
{
public:
    virtual ~object() {};
};

class sys : public object
{
    static sys* fself;
public:
    sys() : flog(0) { fself=this; }
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
    static std::ostream& log()
    {  return fself->flog ? *(fself->flog) : std::cout; }
private:
    std::ostream* flog;
};

#endif // COMMONS_H
