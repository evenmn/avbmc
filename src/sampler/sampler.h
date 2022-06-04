/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.
---------------------------------------------------------------------------- */

#pragma once
#include <vector>
#include <valarray>
#include <string>


class Sampler
{
public:
    Sampler(class System *);

    // declare pure virtual functions
    virtual double w(int) = 0;
    virtual ~Sampler() = default;
    
    // declare public functions
    void sample(int);

    std::string label;

protected:
    class System* system = nullptr;
    class RandomNumberGenerator* rng = nullptr;
};
