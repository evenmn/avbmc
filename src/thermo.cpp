#include "thermo.h"
#include "../box.h"



auto step = [] (class Box* box) -> double {
    return box->step;
};

auto time_ = [] (class Box* box) -> double {
    return box->time;
};

auto atoms = [] (class Box* box) -> double {
    return box->npar;
};

auto types = [] (class Box* box) -> double {
    return box->ntype;
};

auto poteng = [] (class Box* box) -> double {
    return box->poteng;
};

auto kineng = [] (class Box* box) -> double {
    vec vel = box->velocities;
    return sum(sum(vel % vel));
};

auto acceptance_ratio = [] (class Box* box) -> double {
    return box->sampler->acceptance_ratio;
};


Thermo::Thermo(class Box* box_in, const int freq_in, const string filename, const vector<string> outputs_in)
{
    // store box and outputs
    freq = freq_in;
    box = box_in;
    outputs = outputs_in;

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
        else if(i == "kineng"){
            output_functions.push_back(kineng);
        }
        else if(i == "acceptanceratio"){
            output_functions.push_back(acceptance_ratio);
        }
        else{
            cout << "No output style '" + i + "' exists!" << endl;
            exit(0);
        }
    }
}


void Thermo::print_header()
{
    /* Print header of thermo outputs
     */
    for(string i : outputs){
        cout << i << " ";
    }
    cout << endl;
}


void Thermo::print_line()
{
    /* Print thermo output at the current
     * step
     */
    if(box->step % freq == 0){
        for(auto f : output_functions){
            cout << f(box) << " ";
        }
        cout << endl;
    }
}

