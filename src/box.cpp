#include "box.h"

Box::Box(string working_dir_in, double temp_in, double chempot_in) //: integrator(Euler(this)), forcefield(LennardJones(this, ".in")
{
    working_dir = working_dir_in;
    temp = temp_in;
    chempot = chempot_in;

    npar = 0;

    //integrator = new Euler(this);
    //forcefield = new LennardJones(this, ".in");
}


void Box::set_forcefield(class ForceField* forcefield_in)
{
    forcefield = forcefield_in;
}

void Box::set_integrator(class Integrator* integrator_in)
{
    integrator = integrator_in;
}

/*
void Box::set_sampler(class Sampler* sampler_in)
{
    sampler = sampler_in;
}
*/

/*
void Box::add_move(class Moves move, double prob)
{
    moves.push_back(move);
    moves_prob.push_back(prob);
}
*/

void Box::add_particles(const int type, const double mass, const mat position, const mat velocity, const string chem)
{
    /* Add particles of type "type", mass "mass",
     * positions "position" and velocities "velocity"
     */
    int npar_added = position.n_rows;
    npar += npar_added;
    ndim = position.n_cols;
    assert (velocity.n_rows == npar_added);
    
    types = join_cols(types, type * ones(npar_added));
    masses = join_cols(masses, mass * ones(npar_added));
    positions = join_cols(positions, position);
    velocities = join_cols(velocities, velocity);

    vector<string> chem_symbol_vec(npar_added, chem);
    chem_symbol.insert(chem_symbol.end(), chem_symbol_vec.begin(), chem_symbol_vec.end());
}

void Box::write_xyz(string filename)
{
    /* Write an Armadillo matrix containing
     * positions at a given timestep to file with
     * xyz-format. The cube has to have particles on
     * first axis and dimensions on second axis
     */
    vector<string> dims = {"x", "y", "z"};

    ofstream f;

    f.open(filename);
    f << npar << endl;
    f << "symbol";
    for(int j=0; j<ndim; j++){
        f << " " + dims[j];
    }
    f << endl;
    for(int i=0; i<npar; i++){
        f << chem_symbol[i];
        for(int j=0; j<ndim; j++){
            f << " " << positions(i, j);
        }
        f << endl;
    }
    f.close();
}

void Box::snapshot(string filename){
    /* Dump snapshot of system using the
     * "write_xyz"-function. 
     */
    write_xyz(filename);
}


void Box::run_md(int steps)
{
    // Compute initial acceleration
    forcefield->eval_acc(positions, accelerations);

    // Run molecular dynamics simulation
    for(int step=0; step<steps; step++){
        integrator->next_step();
        // dump
        cout << step << endl;
    }
}
