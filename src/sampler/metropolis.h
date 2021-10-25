#include <cmath>
#include "sampler.h"

class Metropolis : public Sampler
{
public:
    Metropolis(class Box* box_in);

    double get_acceptance_prob(class Moves* move, double temp, double chempot);
};
