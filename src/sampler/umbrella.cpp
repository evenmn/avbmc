#include <iostream>
#include <valarray>
#include <functional>

#include "umbrella.h"
#include "../box.h"
#include "../moves/moves.h"


/* --------------------------------------------------------
   Umbrella sampling initialization. This class takes
   only one weight function 'f_in', and is therefore
   basic compared to more advanced Umbrella sampling 
   approaches.
----------------------------------------------------------- */

Umbrella::Umbrella(Box* box_in, std::function<double (int npar)> f_in, const int maxpar_in)
    : Sampler(box_in) 
{
    f = f_in;
    maxpar = maxpar_in;
    tabulated = tabulate(f_in, maxpar_in);
}


/* ----------------------------------------------------------
   Tabulate the weight function 'f_in'. Since the weight 
   function takes integer inputs, it is discrete and can be 
   tabulated. 'maxpar_in' is the maximum int to tabulate
------------------------------------------------------------- */

std::valarray<double> Umbrella::tabulate(std::function<double (int npar)> f_in, const int maxpar_in)
{
    std::valarray<double> tabulated(maxpar_in);
    for(int i=0; i<maxpar; i++){
        tabulated[i] = f_in(i);
    }
    return tabulated;
}


/* ----------------------------------------------------------
   Return value of weight function
------------------------------------------------------------- */

double Umbrella::w(const int npar)
{
    if(npar < maxpar){
        return tabulated[npar];
    }
    else{
        return f(npar);
    }
}
