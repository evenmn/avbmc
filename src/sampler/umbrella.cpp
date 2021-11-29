#include "umbrella.h"
#include "../box.h"
#include "../moves/moves.h"


Umbrella::Umbrella(Box* box_in, std::function<double (int npar)> f_in, const int maxpar_in)
    : Sampler(box_in) 
{
    f = f_in;
    maxpar = maxpar_in;
    tabulated = tabulate(f_in, maxpar_in);
}

std::valarray<double> Umbrella::tabulate(std::function<double (int npar)> f_in, const int maxpar_in)
{
    /* Tabulate the f-function for system sizes up to maxpar
     */
    std::valarray<double> tabulated(maxpar_in);
    for(int i=0; i<maxpar; i++){
        tabulated[i] = f_in(i);
    }
    return tabulated;
}


double Umbrella::w(const int npar)
{
    /* Return weight value
     */
    if(npar < maxpar){
        return tabulated[npar];
    }
    else{
        return f(npar);
    }
}
