/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  The tools needed for input/output is implemented in this file. This includes
  reading and writing xyz-files and writing vectors and arrays to file.
---------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

#include "io.h"
#include "particle.h"


/* ----------------------------------------------------------------------------
   Split string 's' by whitespace.
---------------------------------------------------------------------------- */

std::vector<std::string> split(const std::string &s)
{
    std::stringstream ss(s);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);
    return vstrings;
}


/* ----------------------------------------------------------------------------
   Read xyz-file 'filename' consisting of one time frame into vector of
   particle objects.

    TODO: Allow for files with multiple frames. Then only read last frame
---------------------------------------------------------------------------- */

std::vector<Particle> read_xyz(const std::string &filename)
{
    unsigned int i, ndim, line_num;  //npar;
    std::string line;
    std::vector<Particle> particles;
    std::ifstream f(filename);

    ndim = line_num = 0;  //npar = 0;
    if (f.is_open()) {
        while (std::getline(f, line))
        {
            if (line_num == 0) {
                // line containing number of particles
                //npar = std::stoi(line);
            }
            else if (line_num == 1) {
                // ignore info line
                line_num ++;
                continue;
            }
            else { 
                // coordinate line
                std::vector<std::string> splitted = split(line);
                if (line_num == 2) {
                    // count number of dimensions
                    ndim = splitted.size() - 1;
                }
                std::valarray<double> tmp_r(ndim);
                for (i=0; i<ndim; i++) {
                    tmp_r[i] = std::stod(splitted[i+1]);
                }
                Particle particle(splitted[0], tmp_r);
                particles.push_back(particle);
            }
            line_num ++;
        }
    }
    else {
        std::cout << "Could not open file '" + filename + "'! Aborting." << std::endl;
        exit(0);
    }
    return particles;
}


/* ----------------------------------------------------------------------------
  Write particle properties (i.e. positions) to xyz-file. This is mainly used 
  by the dump class.

  Arguments:
      f : out file stream to write to
      data : matrix of size (npar, noutputs) containing data to write
      nrow : number of particles to write out to file
      ncol : number of outputs
      labels : vector containing labels of all involved particles
      info : information to put in info line, empty by default

---------------------------------------------------------------------------- */

void write_xyz(std::ofstream &f, double **data, const unsigned int nrow,
               const unsigned int ncol, const std::vector<std::string> &labels,
               const std::string &info = "")
{
    unsigned int i, j;

    f << nrow << std::endl;
    f << info << std::endl;;
    for(i=0; i<nrow; i++){
        f << labels[i];
        for(j=0; j<ncol; j++){
            f << " " << data[i][j];
        }
        f << std::endl;
    }
}


/* ----------------------------------------------------------------------------
   Write vector 'vec' to file 'filename' with instances separated by 'delim'.
---------------------------------------------------------------------------- */

void write_vector(std::vector<int> vec, const std::string &filename,
    const std::string &delim)
{
    std::ofstream f(filename);
    for (int element : vec)
    {
        f << element << delim;
    }
    f.close();
}


void write_vector(std::vector<double> vec, const std::string &filename,
    const std::string &delim)
{
    std::ofstream f(filename);
    for (double element : vec)
    {
        f << element << delim;
    }
    f.close();
}


void write_array(int* arr, int length, const std::string &filename,
    const std::string &delim)
{
    int i;
    std::ofstream f(filename);
    for (i=0; i < length; i++)
    {
        f << arr[i] << delim;
    }
    f.close();
}


void write_array(unsigned int* arr, unsigned int length,
    const std::string &filename, const std::string &delim)
{
    unsigned int i;

    std::ofstream f(filename);
    for (i=0; i < length; i++)
    {
        f << arr[i] << delim;
    }
    f.close();
}
