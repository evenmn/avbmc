#include "parser.h"

vector<string> split(const string s)
{
    /* Split string by whitespace
     */
    stringstream ss(s);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> vstrings(begin, end);
    return vstrings;
}


void parser(int argc, char** argv)
{
    /* Parser main function
     */
    if(argc >= 2){
        string filename = argv[1];
        ifstream f(filename);

        if (f.is_open()){ 
            Box box("", 300., 0.);
            string line;
            while (getline(f, line))
            {
                if (line.empty()) {
                    /* Continue if line is blank */
                    continue;
                }
                else if(line.rfind("#", 0) == 0){
                    /* Continue if line is commented out */
                    continue;
                }
                else{
                    vector<string> splitted = split(line);
                    int argc = splitted.size();
                    string cmd_cat = splitted[0];
                    if(cmd_cat == "set"){
                        set(box, splitted, argc);
                    }
                    else if(cmd_cat == "add"){
                        add(box, splitted, argc);
                    }
                    else if(cmd_cat == "take"){
                        take(box, splitted, argc);
                    }
                    else if(cmd_cat == "run"){
                        run(box, splitted, argc);
                    }
                    else if(cmd_cat == "thermo"){
                        thermo(box, splitted, argc);
                    }
                    else if(cmd_cat == "dump"){
                        dump(box, splitted, argc);
                    }
                    else{
                        cout << "There is no category '" + cmd_cat + "'! Aborting." << endl;
                        exit(0);
                    }
                }
            }
        }
        else{
            cout << "Could not open the file!" << endl;
        }
    }
    else{
        cout << "avbmc version 0.0.1" << endl;
        cout << "Author: Even M. Nordhagen" << endl;
        cout << "url: github.com/evenmn/avbmc" << endl;
    }
}

void set(Box& box, const vector<string> splitted, const int argc)
{
    assert(argc > 2);
    string keyword = splitted[1];
    if(keyword == "temp"){
        box.set_temp(stod(splitted[2]));
    }
    else if(keyword == "chempot"){
        box.set_chempot(stod(splitted[2]));
    }
    else if(keyword == "mass"){
        assert(argc > 3);
        string chem_symbol = splitted[2];
        double mass = stod(splitted[3]);
        box.set_mass(chem_symbol, mass);
    }
    else if(keyword == "forcefield"){
        assert(argc > 2);
        string forcefield = splitted[2];
        string paramfile = splitted[3];
        if(forcefield == "lennardjones"){
            box.set_forcefield(new LennardJones(&box, paramfile));
        }
        //else if(forcefield == "vashishta"){
        //    box.set_forcefield(new Vashishta(&box, paramfile));
        //}
        else{
            cout << "Forcefield '" + forcefield + "' is not known! Aborting." << endl;
            exit(0);
        }
    }
    else if(keyword == "integrator"){
        string integrator = splitted[2];
        if(argc == 3){
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box));
            }
            else{
                cout << "Integrator '" + integrator + "' is not known! Aborting." << endl;
                exit(0);
            }
        }
        else if(argc >= 4){
            double dt = stod(splitted[3]);
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box, dt));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box, dt));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box, dt));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box, dt));
            }
            else{
                cout << "Integrator '" + integrator + "' is not known! Aborting." << endl;
                exit(0);
            }
        }
        else{
            cout << "Integrator does not have a sufficient number of arguments! Aborting." << endl;
            exit(0);
        }
    }
    else if(keyword == "sampler"){
        string sampler = splitted[2];
        if(sampler == "metropolis"){
            box.set_sampler(new Metropolis(&box));
        }
        //else if(sampler == "umbrella"){
        //    box.set_sampler(new Umbrella(&box));
        //}
        else {
            cout << "Sampler '" + sampler + "' is not known! Aborting." << endl;
            exit(0);
        }
    }
    else if(keyword == "rng"){
        string rng = splitted[2];
        if(rng == "mersennetwister"){
            box.set_rng(new MersenneTwister());
        }
        else{
            cout << "Random number generator '" + rng + "' is not known! Aborting." << endl;
        }
    }
    else if(keyword == "boundary"){
        assert(argc > 3);
        string boundary = splitted[2];
        if(boundary == "stillinger"){
            double rc = stod(splitted[3]);
            box.set_boundary(new Stillinger(&box, rc));
        }
        else if(boundary == "fixed"){
            if(argc == 4){
                // one box length, indicating one dimension
                double lx = stod(splitted[3]);
                box.set_boundary(new Fixed(&box, lx));
            }
            else if(argc == 5){
                // two box lengths, indicating two dimensions
                double lx = stod(splitted[3]);
                double ly = stod(splitted[4]);
                box.set_boundary(new Fixed(&box, lx, ly));
            }
            else if(argc == 6){
                double lx = stod(splitted[3]);
                double ly = stod(splitted[4]);
                double lz = stod(splitted[5]);
                box.set_boundary(new Fixed(&box, lx, ly, lz));
            }
        }
    }
    else{
        cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
        exit(0);
    }
}


void add(Box& box, const vector<string> splitted, const int argc)
{
    assert(argc > 1);
    string keyword = splitted[1];

    if(keyword == "move"){
        assert(argc > 3);
        string move = splitted[2];
        double prob = stod(splitted[3]);
        if(move == "trans"){
            if(argc == 4){
                box.add_move(new Trans(&box), prob);
            }
            else if(argc >= 5){
                double dx = stod(splitted[4]);
                box.add_move(new Trans(&box, dx), prob);
            }
        }
        else if(move == "transmh"){
            if(argc == 4){
                box.add_move(new TransMH(&box), prob);
            }
            else if(argc == 5){
                double dx = stod(splitted[4]);
                box.add_move(new TransMH(&box, dx), prob);
            }
            else if(argc >= 6){
                double dx = stod(splitted[4]);
                double Ddt = stod(splitted[5]);
                box.add_move(new TransMH(&box, dx, Ddt), prob);
            }
        }
        else{
            cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
            exit(0);
        }
    }
    else if(keyword == "particle"){
        assert(argc > 3);
        if(argc == 4){
            // assuming 1d
            string chem_symbol = splitted[2];
            double x = stod(splitted[3]);
            mat positions_in = {x};
            box.add_particles(chem_symbol, positions_in);
        }
        if(argc == 5){
            // assuming 2d
            string chem_symbol = splitted[2];
            double x = stod(splitted[3]);
            double y = stod(splitted[4]);
            mat positions_in = {x, y};
            box.add_particles(chem_symbol, positions_in);
        }
        else{
            // assuming 3d
            string chem_symbol = splitted[2];
            double x = stod(splitted[3]);
            double y = stod(splitted[4]);
            double z = stod(splitted[5]);
            mat positions_in = {x, y, z};
            box.add_particles(chem_symbol, positions_in);
        }
    }
    else if(keyword == "particles"){
        assert(argc > 3);
        string init_type = splitted[2];
        if(init_type == "fcc"){
            assert(argc > 5);
            string chem_symbol = splitted[3];
            int ncell = stoi(splitted[4]);
            double lenbulk = stod(splitted[5]);
            mat positions_in;
            if(argc == 6){
                positions_in = fcc(ncell, lenbulk);
            }
            else{
                int ndim = stoi(splitted[6]);
                positions_in = fcc(ncell, lenbulk, ndim);
            }
            box.add_particles(chem_symbol, positions_in);
        }
        //else if(init_type == "xyz"){
        //    string filename = splitted[3];
        //    vector<string> chem_symbols_in;
        //    mat positions_in = from_xyz(filename, chem_symbols_in);
        //    box.add_particles(chem_symbols_in, positions_in);
        //}
        else{
            cout << "Particle initialization type '" + init_type + "' is not known! Aborting." << endl;
        }
    }
    else{
        cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
        exit(0);
    }
}

void take(Box& box, const vector<string> splitted, const int argc)
{
    assert(argc > 1);
    string keyword = splitted[1];
    if(keyword == "snapshot"){
        assert(argc > 2);
        string filename = splitted[2];
        box.snapshot(filename);
    }
    else{
        cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
        exit(0);
    }
}

void run(Box& box, const vector<string> splitted, const int argc)
{
    assert(argc > 2);
    string keyword = splitted[1];
    int nsteps = stod(splitted[2]);
    if(keyword == "mc"){
        assert(argc > 3);
        int nmoves = stod(splitted[3]);
        box.run_mc(nsteps, nmoves);
    }
    else if(keyword == "md"){
        box.run_md(nsteps);
    }
    else{
        cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
        exit(0);
    }
}

void thermo(Box& box, const vector<string> splitted, const int argc)
{
    /* Thermo output
     */
    assert(argc > 3);
    int freq = stoi(splitted[1]);
    string filename = splitted[2];
    vector<string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    box.set_thermo(freq, filename, outputs);
}

void dump(Box& box, const vector<string> splitted, const int argc)
{
    /* Dump output
     */
    assert(argc > 3);
    int freq = stoi(splitted[1]);
    string filename = splitted[2];
    vector<string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    box.set_dump(freq, filename, outputs);
}
