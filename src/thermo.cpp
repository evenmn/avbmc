#include "thermo.h"
#include "../box.h"



auto step = [] (Box* box) -> double {
    return box->step;
};

auto time_ = [] (Box* box) -> double {
    return box->time;
};

auto atoms = [] (Box* box) -> double {
    return box->npar;
};

auto types = [] (Box* box) -> double {
    return box->ntype;
};

auto poteng = [] (Box* box) -> double {
    return box->poteng;
};
/*
auto kineng = [] (Box* box) -> double {
    vec vel = box->velocities;
    return sum(sum(vel % vel));
};
*/
auto acceptance_ratio = [] (Box* box) -> double {
    return box->sampler->acceptance_ratio;
};

auto move_idx = [] (Box* box) -> double {
    return box->sampler->move_idx;
};


Thermo::Thermo(Box* box_in, const int freq_in, const std::string filename, const std::vector<std::string> outputs_in)
{
    // store box and outputs
    freq = freq_in;
    box = box_in;
    outputs = outputs_in;

    // open file
    f.open(filename);

    // fill vector with output functions
    for(string i : outputs_in){
        if(i == "step"){
            output_functions.push_back(step);
        }
        else if(i == "time"){
            output_functions.push_back(time_);
        }
        else if(i == "atoms"){
            output_functions.push_back(atoms);
        }
        else if(i == "types"){
            output_functions.push_back(types);
        }
        else if(i == "poteng"){
            output_functions.push_back(poteng);
        }
        /*
        else if(i == "kineng"){
            output_functions.push_back(kineng);
        }
        */
        else if(i == "acceptanceratio"){
            output_functions.push_back(acceptance_ratio);
        }
        else if(i == "move"){
            output_functions.push_back(move_idx);
        }
        else{
            std::cout << "No output style '" + i + "' exists! Aborting." << std::endl;
            exit(0);
        }
    }
}


void Thermo::print_header()
{
    /* Print header of thermo outputs
     */
    f << "# ";
    for(std::string i : outputs){
        std::cout << i << " ";
        f << i << " ";
    }
    std::cout << std::endl;
    f << std::endl;
}


void Thermo::print_line()
{
    /* Print thermo output at the current
     * step
     */
    if(box->step % freq == 0){
        for(auto func : output_functions){
            std::cout << func(box) << " ";
            f << func(box) << " ";
        }
        std::cout << std::endl;
        f << std::endl;
    }
}

Thermo::~Thermo()
{
    f.close();
}
