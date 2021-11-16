#include "umbrella.h"
#include "../box.h"
#include "../moves/moves.h"


Umbrella::Umbrella(class Box* box_in, std::function<double (int npar)> f_in, const int maxpar_in)
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


double Umbrella::get_acceptance_prob(class Moves* move, double temp, double chempot)
{
    /* Get acceptance probability of biased sampling
     */
    int npar = box->npar;
    if(npar < maxpar){
        return move->accept() * tabulated[npar] * std::exp(-(du + chempot)/temp);
    }
    else{
        return move->accept() * f(npar) * std::exp(-(du + chempot)/temp);
    }
}
