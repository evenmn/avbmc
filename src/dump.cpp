#include "dump.h"
#include "../box.h"


auto x = [] (class Box* box) -> mat {
    return box->positions.col(0);
};

auto y = [] (class Box* box) -> mat {
    return box->positions.col(1);
};

auto z = [] (class Box* box) -> mat {
    return box->positions.col(2);
};

auto xy = [] (class Box* box) -> mat {
    return box->positions.cols(0, 1);
};

auto xyz = [] (class Box* box) -> mat {
    return box->positions;
};

auto vx = [] (class Box* box) -> mat {
    return box->velocities.col(0);
};

auto vy = [] (class Box* box) -> mat {
    return box->velocities.col(1);
};

auto vz = [] (class Box* box) -> mat {
    return box->velocities.col(2);
};

auto vxvy = [] (class Box* box) -> mat {
    return box->velocities.cols(0, 1);
};

auto vxvyvz = [] (class Box* box) -> mat {
    return box->velocities;
};

auto ax = [] (class Box* box) -> mat {
    return box->accelerations.col(0);
};

auto ay = [] (class Box* box) -> mat {
    return box->accelerations.col(1);
};

auto az = [] (class Box* box) -> mat {
    return box->accelerations.col(2);
};

auto axay = [] (class Box* box) -> mat {
    return box->accelerations.cols(0, 1);
};

auto axayaz = [] (class Box* box) -> mat {
    return box->accelerations;
};


Dump::Dump(class Box* box_in, const int freq_in, const string filename, const vector<string> outputs_in)
{
    // store box and outputs
    freq = freq_in;
    box = box_in;

    // sort outputs
    outputs = outputs_in;  //sort(outputs_in.begin(), outputs_in.end());

    // fill vector with output functions
    info_line = "";
    for(string i : outputs){
        info_line += i + " ";
        if(i == "ax"){
            output_functions.push_back(ax);
        }
        else if(i == "ay"){
            output_functions.push_back(ay);
        }
        else if(i == "az"){
            output_functions.push_back(az);
        }
        else if(i == "axay"){
            output_functions.push_back(axay);
        }
        else if(i == "axayaz"){
            output_functions.push_back(axayaz);
        }
        if(i == "vx"){
            output_functions.push_back(vx);
        }
        else if(i == "vy"){
            output_functions.push_back(vy);
        }
        else if(i == "vz"){
            output_functions.push_back(vz);
        }
        else if(i == "vxvy"){
            output_functions.push_back(vxvy);
        }
        else if(i == "vxvyvz"){
            output_functions.push_back(vxvyvz);
        }
        if(i == "x"){
            output_functions.push_back(x);
        }
        else if(i == "y"){
            output_functions.push_back(y);
        }
        else if(i == "z"){
            output_functions.push_back(z);
        }
        else if(i == "xy"){
            output_functions.push_back(xy);
        }
        else if(i == "xyz"){
            output_functions.push_back(xyz);
        }
        else{
            cout << "No output style '" + i + "' exists!" << endl;
            exit(0);
        }
    }

    // open file and write
    f.open(filename);
}


void Dump::print_frame()
{
    /* Print thermo output at the current
     * step
     */
    if(box->step % freq == 0){
        mat dump_data;
        for(auto func : output_functions){
            dump_data = join_rows(dump_data, func(box));
        }
        write_xyz(f, dump_data, box->chem_symbols, info_line);
    }
}


Dump::~Dump() { f.close(); }
