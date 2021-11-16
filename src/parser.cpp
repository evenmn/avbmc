#include "parser.h"

void parser(int argc, char** argv)
{
    /* Parser main function
     */
    if(argc >= 2){
        std::string filename = argv[1];
        std::ifstream f(filename);

        if (f.is_open()){ 
            Box box("", 300., 0.);
            std::string line;
            while (std::getline(f, line))
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
                    std::vector<std::string> splitted = split(line);
                    int argc = splitted.size();
                    std::string cmd_cat = splitted[0];
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
                        std::cout << "There is no category '" + cmd_cat + "'! Aborting." << std::endl;
                        exit(0);
                    }
                }
            }
        }
        else{
            std::cout << "Could not open the file!" << std::endl;
            exit(0);
        }
    }
    else{
        std::cout << "avbmc version 0.0.1" << std::endl;
        std::cout << "Author: Even M. Nordhagen" << std::endl;
        std::cout << "url: github.com/evenmn/avbmc" << std::endl;
    }
}

void set(Box& box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 2);
    std::string keyword = splitted[1];
    if(keyword == "temp"){
        box.set_temp(std::stod(splitted[2]));
    }
    else if(keyword == "chempot"){
        box.set_chempot(std::stod(splitted[2]));
    }
    else if(keyword == "mass"){
        assert(argc > 3);
        std::string label = splitted[2];
        double mass = std::stod(splitted[3]);
        box.set_mass(label, mass);
    }
    else if(keyword == "forcefield"){
        assert(argc > 2);
        std::string forcefield = splitted[2];
        std::string paramfile = splitted[3];
        if(forcefield == "lennardjones"){
            box.set_forcefield(new LennardJones(&box, paramfile));
        }
        //else if(forcefield == "vashishta"){
        //    box.set_forcefield(new Vashishta(&box, paramfile));
        //}
        else{
            std::cout << "Forcefield '" + forcefield + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    /*
    else if(keyword == "integrator"){
        std::string integrator = splitted[2];
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
                std::cout << "Integrator '" + integrator + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else if(argc >= 4){
            double dt = std::stod(splitted[3]);
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
                std::cout << "Integrator '" + integrator + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else{
            std::cout << "Integrator does not have a sufficient number of arguments! Aborting." << std::endl;
            exit(0);
        }
    }
    */
    else if(keyword == "sampler"){
        std::string sampler = splitted[2];
        if(sampler == "metropolis"){
            box.set_sampler(new Metropolis(&box));
        }
        else if(sampler == "umbrella"){
            assert(argc > 3);
            std::string umbrella_type = splitted[3];
            if(umbrella_type == "poly"){
                assert(argc > 6);
                double a = std::stod(splitted[4]);
                double b = std::stod(splitted[5]);
                double c = std::stod(splitted[6]);
                auto f = [a, b, c] (const int x) { return (a * x * x + b * x + c); };
                box.set_sampler(new Umbrella(&box, f));
            }
            else if(umbrella_type == "square"){
                assert(argc > 5);
                double nc = std::stod(splitted[4]);
                double k = std::stod(splitted[5]);
                auto f = [nc, k] (const int n) { return (k * (n - nc) * (n - nc)); };
                box.set_sampler(new Umbrella(&box, f));
            }
            else{
                std::cout << "Umbrella type '" + umbrella_type + "' is not known! Aborting." << std::endl;
                exit(0);
            }
        }
        else {
            std::cout << "Sampler '" + sampler + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "rng"){
        std::string rng = splitted[2];
        if(rng == "mersennetwister"){
            box.set_rng(new MersenneTwister());
        }
        else{
            std::cout << "Random number generator '" + rng + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "boundary"){
        assert(argc > 3);
        std::string boundary = splitted[2];
        if(boundary == "stillinger"){
            double rc = std::stod(splitted[3]);
            box.set_boundary(new Stillinger(&box, rc));
        }
        /*
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
        */
        else{
            std::cout << "Boundary '" + boundary + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    /*
    else if(keyword == "velocity"){
        assert(argc > 3);
        std::string init_style = splitted[2];
        if(init_style == "temp"){
            double temp = std::stod(splitted[3]);
            box.velocity = new Temp(box.rng, temp);
            box.velocities = box.velocity->get_velocity(box.npar, box.ndim);
        }
        else if(init_style == "gauss"){
            assert(argc > 4);
            double mean = std::stod(splitted[3]);
            double var = std::stod(splitted[4]);
            box.velocity = new Gauss(box.rng, mean, var);
            box.velocities = box.velocity->get_velocity(box.npar, box.ndim);
        }
        else{
            std::cout << "Velocity style '" + init_style + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    */
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}


void add(Box& box, const vector<string> splitted, const int argc)
{
    assert(argc > 1);
    std::string keyword = splitted[1];

    if(keyword == "move"){
        assert(argc > 3);
        std::string move = splitted[2];
        double prob = std::stod(splitted[3]);
        if(move == "trans"){
            if(argc == 4){
                box.add_move(new Trans(&box), prob);
            }
            else{
                double dx = std::stod(splitted[4]);
                box.add_move(new Trans(&box, dx), prob);
            }
        }
        /*
        else if(move == "transmh"){
            if(argc == 4){
                box.add_move(new TransMH(&box), prob);
            }
            else if(argc == 5){
                double dx = std::stod(splitted[4]);
                box.add_move(new TransMH(&box, dx), prob);
            }
            else{
                double dx = std::stod(splitted[4]);
                double Ddt = std::stod(splitted[5]);
                box.add_move(new TransMH(&box, dx, Ddt), prob);
            }
        }
        */
        else if(move == "avbmc"){
            if(argc == 4){
                box.add_move(new AVBMC(&box), prob);
            }
            else if(argc == 5){
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMC(&box, r_below), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                double r_above = std::stod(splitted[5]);
                box.add_move(new AVBMC(&box, r_below, r_above), prob);
            }
        }
        else if(move == "avbmcin"){
            if(argc == 4){
                box.add_move(new AVBMCIn(&box), prob);
            }
            else if(argc == 5){
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMCIn(&box, r_below), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                double r_above = std::stod(splitted[5]);
                box.add_move(new AVBMCIn(&box, r_below, r_above), prob);
            }
        }
        else if(move == "avbmcout"){
            if(argc == 4){
                box.add_move(new AVBMCOut(&box), prob);
            }
            else{
                double r_below = std::stod(splitted[4]);
                box.add_move(new AVBMCOut(&box, r_below), prob);
            }
        }
        else{
            std::cout << "Move '" + move + "' is not known! Aborting." << std::endl;
            exit(0);
        }
    }
    else if(keyword == "particle"){
        assert(argc > 3);
        Particle* particle;
        particle->label = splitted[2];
        std::string label = splitted[2];
        double x = std::stod(splitted[3]);
        if(argc == 4){
            // assuming 1d
            std::valarray<double> r = {x};  // initializer list supported after C++11
            box.add_particle(label, r);
        }
        if(argc == 5){
            // assuming 2d
            double y = std::stod(splitted[4]);
            std::valarray<double> r = {x, y};  // initializer list supported after C++11
            box.add_particle(label, r);
        }
        else{
            // assuming 3d
            double y = std::stod(splitted[4]);
            double z = std::stod(splitted[5]);
            std::cout << "Here3" << std::endl;
            std::valarray<double> r = {x, y, z}; // initializer list supported after C++11
            std::cout << "Here32" << std::endl;
            std::cout << particle->r.size() << std::endl;
            particle->r.resize(3);
            std::cout << particle->r.size() << std::endl;
            std::cout << "Here33" << std::endl;
            //particle->r = r;
            std::cout << "Here4" << std::endl;
            box.add_particle(particle);
            std::cout << "Here5" << std::endl;
        }
    }
    else if(keyword == "particles"){
        assert(argc > 3);
        std::string init_type = splitted[2];
        if(init_type == "fcc"){
            assert(argc > 5);
            std::string label = splitted[3];
            int ncell = std::stoi(splitted[4]);
            double lenbulk = std::stod(splitted[5]);
            std::vector<Particle *> particles_in;
            if(argc == 6){
                particles_in = fcc(ncell, lenbulk);
            }
            else{
                int ndim = std::stoi(splitted[6]);
                particles_in = fcc(ncell, lenbulk, ndim);
            }
            for(Particle *particle : particles_in){
                particle->label = label;
            }
            box.add_particles(particles_in);
        }
        else if(init_type == "xyz"){
            std::string filename = splitted[3];
            std::vector<Particle *> particles_in = read_xyz(filename);
            box.add_particles(particles_in);
        }
        else{
            std::cout << "Particle initialization type '" + init_type + "' is not known! Aborting." << std::endl;
        }
    }
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}

void take(Box& box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 1);
    std::string keyword = splitted[1];
    if(keyword == "snapshot"){
        assert(argc > 2);
        std::string filename = splitted[2];
        box.snapshot(filename);
    }
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}

void run(Box& box, const std::vector<std::string> splitted, const int argc)
{
    assert(argc > 2);
    std::string keyword = splitted[1];
    int nsteps = std::stod(splitted[2]);
    if(keyword == "mc"){
        assert(argc > 3);
        int nmoves = std::stod(splitted[3]);
        box.run_mc(nsteps, nmoves);
    }
    /*
    else if(keyword == "md"){
        box.run_md(nsteps);
    }
    */
    else{
        std::cout << "Keyword '" + keyword + "' is not known! Aborting." << std::endl;
        exit(0);
    }
}

void thermo(Box& box, const std::vector<std::string> splitted, const int argc)
{
    /* Thermo output
     */
    assert(argc > 3);
    int freq = stoi(splitted[1]);
    std::string filename = splitted[2];
    std::vector<std::string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    box.set_thermo(freq, filename, outputs);
}

void dump(Box& box, const std::vector<std::string> splitted, const int argc)
{
    /* Dump output
     */
    assert(argc > 3);
    int freq = std::stoi(splitted[1]);
    std::string filename = splitted[2];
    std::vector<std::string> outputs;
    for(int i=3; i<argc; i++){
        outputs.push_back(splitted[i]);
    }
    box.set_dump(freq, filename, outputs);
}
