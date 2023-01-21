/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
  Output system and box properties to a file. This includes number of atoms,
  steps, potential energy and so on. Properties are written to a file format 
  that can be easily read by numpy.loadtxt or pandas.read_csv
---------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <string>
#include <functional>


#include "box.h"
#include "system.h"
#include "thermo.h"


/* ----------------------------------------------------------------------------
   Get current step
---------------------------------------------------------------------------- */

auto step = [] (Box* box) -> double {
    return box->system->step;
};


/* ----------------------------------------------------------------------------
   Get current time (for molecular dynamics)
---------------------------------------------------------------------------- */

auto time_ = [] (Box* box) -> double {
    return box->time;
};


/* ----------------------------------------------------------------------------
   Get current number of atoms
---------------------------------------------------------------------------- */

auto atoms = [] (Box* box) -> double {
    return box->npar;
};


/* ----------------------------------------------------------------------------
   Get current number of atom types
---------------------------------------------------------------------------- */

auto types = [] (Box* box) -> double {
    return box->ntype;
};

/* ----------------------------------------------------------------------------
   Get current potential energy
---------------------------------------------------------------------------- */

auto poteng = [] (Box* box) -> double {
    return box->poteng;
};


/* ----------------------------------------------------------------------------
   Get current kinetic energy (for molecular dynamics simulations)
---------------------------------------------------------------------------- */
/*
auto kineng = [] (Box* box) -> double {
    vec vel = box->velocities;
    return sum(sum(vel % vel));
};
*/

/* ----------------------------------------------------------------------------
   Get acceptance ratio of current cycle (for Monte Carlo simulations)
---------------------------------------------------------------------------- */
/*
auto acceptance_ratio = [] (Box* box) -> double {
    return box->sampler->acceptance_ratio;
};
*/

/* ----------------------------------------------------------------------------
   Get index of current move
---------------------------------------------------------------------------- */
/*
auto move_idx = [] (Box* box) -> double {
    return box->sampler->move_idx;
};
*/

/* ----------------------------------------------------------------------------
   Constructor of thermo class. Takes dump frequency, 'freq_in', thermo output
   filename 'filename' and vector of quantities to output 'outputs_in' as
   arguments.
---------------------------------------------------------------------------- */

Thermo::Thermo(Box* box_in, const int freq_in, const std::string &filename,
               const std::vector<std::string> &outputs_in) : outputs(outputs_in)
{
    // store box and outputs
    freq = freq_in;
    box = box_in;

    if (freq > 0) {
        // open file
        f.open(filename);
        print_header();

        // fill vector with output functions
        for (std::string i : outputs_in) {
            if (i == "step") {
                output_functions.push_back(step);
            }
            else if (i == "time") {
                output_functions.push_back(time_);
            }
            else if (i == "atoms") {
                output_functions.push_back(atoms);
            }
            else if (i == "types") {
                output_functions.push_back(types);
            }
            else if (i == "poteng") {
                output_functions.push_back(poteng);
            }
            /*
            else if (i == "kineng") {
                output_functions.push_back(kineng);
            }
            else if (i == "acceptanceratio") {
                output_functions.push_back(acceptance_ratio);
            }
            else if (i == "move") {
                output_functions.push_back(move_idx);
            }
            */
            else {
                std::cout << "No output style '" + i + "' exists! Aborting." << std::endl;
                exit(0);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
   Print header of thermo outputs
---------------------------------------------------------------------------- */

void Thermo::print_header()
{
    f << "# ";
    for (std::string i : outputs) {
        f << i << " ";
    }
    f << std::endl;
}


/* ----------------------------------------------------------------------------
   Print thermo outputs of current step
---------------------------------------------------------------------------- */

void Thermo::print_line (const std::size_t &step)
{
    if (freq > 0 && step % freq == 0) {
        //std::cout << std::fixed; // << std::setprecision(6);
        f << std::scientific;
        for (auto func : output_functions) {
            f << func(box) << "\t ";
        }
        f << std::endl;
    }
}


/* ----------------------------------------------------------------------------
  Thermo destructor, closing output file.
---------------------------------------------------------------------------- */

Thermo::~Thermo()
{
    f.close();
}
