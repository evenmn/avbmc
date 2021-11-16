#include "io.h"
#include "particle.h"


std::vector<std::string> split(const std::string s)
{
    /* Split string by whitespace
     */
    std::stringstream ss(s);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    return vstrings;
}

std::vector<Particle *> read_xyz(const std::string filename)
{
    /* Read xyz-file consisting of one time frame into 
     * vector of particle objects
     */
    int npar, ndim;
    std::vector<Particle *> particles;
    std::ifstream f(filename);
    if(f.is_open()){
        std::string line;
        int line_num = 0;
        while(std::getline(f, line))
        {
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
                Particle *particle;
                particle->label = splitted[0];
                std::valarray<double> tmp_r(ndim);
                for(int i=0; i<ndim; i++){
                    tmp_r[i] = std::stod(splitted[i+1]);
                }
                particle->r = tmp_r;
            }
            line_num ++;
        }
    }
    else{
        cout << "Could not open file '" + filename + "'! Aborting." << endl;
        exit(0);
    }
    return particles;
}

/*
void write_xyz(const string filename, const mat positions, const vector<string> chem_symbols, const string info, const bool append)
{
    // 
    //
    int npar = positions.n_rows;
    int ndim = positions.n_cols;

    ofstream f;
    if(append){
        f.open(filename, ofstream::out | ofstream::app);
    }
    else{
        f.open(filename, ofstream::out);
    }
    f << npar << endl;
    f << info << endl;;
    for(int i=0; i<npar; i++){
        f << chem_symbols[i];
        for(int j=0; j<ndim; j++){
            f << " " << positions(i, j);
        }
        f << endl;
    }
    f.close();
}
*/

void write_xyz(std::ofstream &f, double **data, const int nrow, const int ncol, const std::vector<std::string> labels, const std::string info)
{
    /* Write matrix to file, assuming that all valarrays inside vector 
     * has the same length
     */

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
