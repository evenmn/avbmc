/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Basic umbrella sampling with a single umbrella. First proposed by Torrie
  and Valleau, "Nonphysical sampling distributions in Monte Carlo free-energy
  estimation: Umbrella sampling" (1977)
---------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>
#include <functional>

#include "umbrella.h"
#include "../system.h"


/* ----------------------------------------------------------------------------
   Umbrella sampling initialization. This class takes only one weight function
   'f_in', and is therefore basic compared to more advanced Umbrella sampling 
   approaches.
---------------------------------------------------------------------------- */

Umbrella::Umbrella(System* system_in, std::function<double (int npar)> f_in,
    int maxpar_in)
    : Sampler(system_in), f(f_in)
{
    maxpar = maxpar_in;
    tabulated = tabulate(f_in, maxpar_in);
    label = "Umbrella";
}


Umbrella::Umbrella(System* system_in, std::valarray<double> tabulated_in)
    : Sampler(system_in), tabulated(tabulated_in)
{
    maxpar = tabulated.size() - 1;
    f = [](int n) -> double {
        std::cout << "Number of particles exceeded tabulated umbrellas! "
                  << "Aborting." << std::endl;
        exit(0);
        return 0.;
    };
    label = "Umbrella";
}


/* ----------------------------------------------------------------------------
   Tabulate the weight function 'f_in'. Since the weight function takes integer
   inputs, it is discrete and can be tabulated. 'maxpar_in' is the maximum int
   to tabulate.
---------------------------------------------------------------------------- */

std::valarray<double> Umbrella::tabulate(std::function<double (int npar)> f_in,
    int maxpar_in)
{
    std::valarray<double> tabulated(maxpar_in);
    for (int i=0; i<maxpar; i++) {
        tabulated[i] = f_in(i);
    }
    return tabulated;
}


/* ----------------------------------------------------------------------------
   Return value of weight function
---------------------------------------------------------------------------- */

double Umbrella::w(const int npar)
{
    if (npar < maxpar) {
        return tabulated[npar];
    }
    else {
        return f(npar);
    }
}
