#include "mersennetwister.h"

MersenneTwister::MersenneTwister()
    : RandomNumberGenerator()
{}

double MersenneTwister::next_gaussian(double mean, double variance)
{
    /* Returns a double drawn from a Gaussian 
     * distribution with mean 'mean' and variance
     * 'variance'
     */
    normal_distribution<double> dis(mean, variance);
    return dis(generator);
}

int MersenneTwister::next_int(int upper_limit)
{
    /* Returns an integer between 0 and 'upper_limit'
     */
    uniform_int_distribution<int> dis(0, upper_limit - 1);
    return dis(generator);
}

double MersenneTwister::next_double()
{
    /* Returns a double between 0 and 1 drawn 
     * from a uniform distribution
     */
    std::uniform_real_distribution<double> dis(0, 1);
    return dis(generator);
}

/*
auto MersenneTwister::choice(vector<typename iterator_traits<Iter>::value_type> samples, vector<double> probabilities)
{
    // Inspired by numpy.random.choice
    //
{
    // assert that sum of probabilities is 1
    double sum_probs = 0;
    for (double prob : probabilities){
        sum_probs += prob;
    }
    assert((sum_probs - 1.0) < 0.01);
    assert(samples.size() == probabilities.size());

    // pick sample based on probabilities
    double r = next_double();
    double p = 0.0;
    for (int i=0; i<samples.size(); i++) {
        p += probabilities[i];
        if (r>=p) return samples[i];
    }
}
*/

int MersenneTwister::choice(vector<double> probabilities)
{
    /* Inspired by numpy.random.choice, but simplified as we always
     * want to draw only one sample. The samples are
     * also always range(len(probabilities), so we don't
     * need to take them as an argument
     */

    // assert that sum of probabilities is 1
    double sum_probs = 0;
    for (double prob : probabilities){
        sum_probs += prob;
    }
    assert(abs(sum_probs - 1.0) < 0.01);

    double r = next_double();
    double p = 0.0;
    for (int i=0; i<probabilities.size(); i++) {
        p += probabilities[i];
        if (r<=p) return i;
    }
    cout << "outside" << endl;
    return 0;
}