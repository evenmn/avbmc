#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "io.h"
#include "particle.h"


/* ----------------------------------------------------
   Split string by whitespace
------------------------------------------------------- */

std::vector<std::string> split(const std::string s)
{
    std::stringstream ss(s);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    return vstrings;
}


/* -----------------------------------------------------
   Read xyz-file consisting of one time frame into 
   vector of particle objects 
-------------------------------------------------------- */

std::vector<Particle *> read_xyz(const std::string filename)
{
    int npar, ndim;
    std::cout << "read_xyz" << std::endl;
    std::vector<Particle *> particles;
    std::ifstream f(filename);
    if(f.is_open()){
        std::string line;
        int line_num = 0;
        while(std::getline(f, line))
        {
            std::cout << line << std::endl;
            if(line_num == 0){
                // line containing number of particles
                npar = std::stoi(line);
            }
            else if(line_num == 1){
                // ignore info line
                line_num ++;
                continue;
            }
            else{ 
                // coordinate line
                std::vector<std::string> splitted = split(line);
                if(line_num == 2){
                    // count number of dimensions
                    ndim = splitted.size() - 1;
                }
                std::valarray<double> tmp_r(ndim);
                for(int i=0; i<ndim; i++){
                    tmp_r[i] = std::stod(splitted[i+1]);
                }
                Particle *particle = new Particle(splitted[0], tmp_r);
                particles.push_back(particle);
            }
            line_num ++;
        }
    }
    else{
        cout << "Could not open file '" + filename + "'! Aborting." << endl;
        exit(0);
    }
    std::cout << "end" << std::endl;
    return particles;
}


/* -----------------------------------------------------------------
  Write matrix to file, assuming that all valarrays inside vector 
  has the same length
-------------------------------------------------------------------- */

void write_xyz(std::ofstream &f, double **data, const int nrow, const int ncol, 
               const std::vector<std::string> labels, const std::string info)
{
    f << nrow << endl;
    f << info << endl;;
    for(int i=0; i<nrow; i++){
        f << labels[i];
        for(int j=0; j<ncol; j++){
            f << " " << data[i][j];
        }
        f << endl;
    }
}
