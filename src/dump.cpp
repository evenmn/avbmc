#include "dump.h"
#include "../box.h"
#include "../particle.h"


auto x = [] (class Box *box) -> double ** {
    double** data = new double*[box->npar];
    for(int i = 0; i<box->npar; i++){
        data[i] = new double[1];
        data[i][0] = box->particles[i]->r[0];
    }
    return data;
};

auto y = [] (class Box *box) -> double ** {
    double** data = new double*[box->npar];
    for(int i = 0; i<box->npar; i++){
        data[i] = new double[1];
        data[i][0] = box->particles[i]->r[1];
    }
    return data;
};

auto z = [] (class Box *box) -> double ** {
    double** data = new double*[box->npar];
    for(int i = 0; i<box->npar; i++){
        data[i] = new double[1];
        data[i][0] = box->particles[i]->r[2];
    }
    return data;
};

auto xy = [] (class Box* box) -> double ** {
    double** data = new double*[box->npar];
    for(int i = 0; i<box->npar; i++){
        data[i] = new double[2];
        data[i][0] = box->particles[i]->r[0];
        data[i][1] = box->particles[i]->r[1];
    }
    return data;
};

auto xyz = [] (class Box* box) -> double ** {
    double** data = new double*[box->npar];
    for(int i = 0; i<box->npar; i++){
        data[i] = new double[3];
        data[i][0] = box->particles[i]->r[0];
        data[i][1] = box->particles[i]->r[1];
        data[i][2] = box->particles[i]->r[2];
    }
    return data;
};

/*
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
*/

Dump::Dump(class Box* box_in, const int freq_in, const std::string filename, const std::vector<std::string> outputs_in)
{
    // store box and outputs
    freq = freq_in;
    box = box_in;

    if(freq != 0){
        // sort outputs
        outputs = outputs_in;  //sort(outputs_in.begin(), outputs_in.end());

        // fill vector with output functions
        info_line = "";
        nvar = 0;
        nvars.clear();
        for(std::string i : outputs){
            info_line += i + " ";
            /*
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
            else if(i == "vx"){
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
            */
            if(i == "x"){
                nvar ++;
                nvars.push_back(1);
                output_functions.push_back(x);
            }
            else if(i == "y"){
                nvar ++;
                nvars.push_back(1);
                output_functions.push_back(y);
            }
            else if(i == "z"){
                nvar ++;
                nvars.push_back(1);
                output_functions.push_back(z);
            }
            else if(i == "xy"){
                nvar += 2;
                nvars.push_back(2);
                output_functions.push_back(xy);
            }
            else if(i == "xyz"){
                nvar += 3;
                nvars.push_back(3);
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
}


void Dump::print_frame()
{
    /* Print thermo output at the current
     * step
     */
    if(freq != 0){
        if(box->step % freq == 0){
            // create array for all data to dump
            double** dump_data = new double*[box->npar];
            for(int i = 0; i<box->npar; i++){
                dump_data[i] = new double[nvar];
            }

            // fill up array 
            int o = 0;
            int cum_nvar = 0;
            for(auto func : output_functions){
                double **tmp_data = func(box);
                for(int i=0; i<box->npar; i++){
                    for(int j=0; j<nvars[o]; j++){
                        dump_data[i][cum_nvar + j] = tmp_data[i][j];
                    }
                    free(tmp_data[i]);
                }
                free(tmp_data);
                cum_nvar += nvars[o];
                o ++;
            }

            // get labels
            std::vector<std::string> labels;
            for(Particle* particle : box->particles)
                labels.push_back(particle->label);

            // write to file
            write_xyz(f, dump_data, box->npar, nvar, labels, info_line);

            // free memory
            for(int i = 0; i<box->npar; i++){
                free(dump_data[i]);
            }
            free(dump_data);
        }
    }
}


Dump::~Dump() { f.close(); }
