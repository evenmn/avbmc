#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>        // automatic conversion between vector/valarray and numpy array
#include <pybind11/functional.h> // automatic conversion between functional and python function

#include "../src/io.h"
#include "../src/box.h"
#include "../src/system.h"
#include "../src/particle.h"
#include "../src/thermo.h"
#include "../src/dump.h"

#include "../src/rng/rng.h"
#include "../src/rng/mersennetwister.h"

#include "../src/sampler/sampler.h"
#include "../src/sampler/metropolis.h"
#include "../src/sampler/umbrella.h"

#include "../src/forcefield/idealgas.h"
#include "../src/forcefield/lennardjones.h"
#include "../src/forcefield/vashishta.h"

#include "../src/moves/moves.h"
#include "../src/moves/trans.h"
#include "../src/moves/transmh.h"
#include "../src/moves/avbmc.h"
#include "../src/moves/avbmcin.h"
#include "../src/moves/avbmcout.h"
#include "../src/moves/avbmcmol.h"
#include "../src/moves/avbmcmolin.h"
#include "../src/moves/avbmcmolout.h"

#include "../src/boundary/boundary.h"
#include "../src/boundary/open.h"
#include "../src/boundary/periodic.h"

#include "../src/constraint/constraint.h"
#include "../src/constraint/minneigh.h"
#include "../src/constraint/maxneigh.h"
#include "../src/constraint/mindistance.h"
#include "../src/constraint/maxdistance.h"
#include "../src/constraint/stillinger.h"

#include "../src/integrator/integrator.h"
#include "../src/integrator/euler.h"
#include "../src/integrator/eulercromer.h"
#include "../src/integrator/velocityverlet.h"
#include "../src/integrator/rungekutta4.h"


namespace py = pybind11;
using namespace pybind11::literals;


/* ----------------------------------------------------------------------------
   Python wrapper
------------------------------------------------------------------------------- */

PYBIND11_MODULE(avbmc, m) {

    // Particle
    py::class_<Particle>(m, "Particle", py::dynamic_attr())
        .def(py::init<std::string, std::valarray<double> >())
        .def_readonly("label", &Particle::label)
        .def_readwrite("r", &Particle::r)
        .def_readwrite("v", &Particle::v)
        .def_readwrite("f", &Particle::f);

    // Initialize position
    m.def("read_xyz", &read_xyz);

    // Random Number Generator
    py::class_<RandomNumberGenerator>(m, "RNG")
        .def_readonly("label", &RandomNumberGenerator::label);
    py::class_<MersenneTwister, RandomNumberGenerator>(m, "MersenneTwister")
        .def(py::init<>())
        .def("next_int", &MersenneTwister::next_int)
        .def("next_double", &MersenneTwister::next_double)
        .def("next_gaussian", &MersenneTwister::next_gaussian)
        .def("choice", &MersenneTwister::choice);

    // System
    py::class_<System>(m, "System")
        .def(py::init<std::string, bool>(), py::arg("working_dir") = "", py::arg("initialize") = true)
        .def("set_sampler", py::overload_cast<Sampler *>(&System::set_sampler))
        .def("set_sampler", py::overload_cast<const std::string &, std::function<double(int)> >(&System::set_sampler))
        .def("set_rng", py::overload_cast<RandomNumberGenerator *>(&System::set_rng))
        .def("set_rng", py::overload_cast<const std::string &>(&System::set_rng))
        .def("set_temp", &System::set_temp)
        .def("set_chempot", &System::set_chempot)
        .def("add_move", py::overload_cast<Moves *, double>(&System::add_move), "Add move to the list of moves",
            py::arg("system"), py::arg("prob") = 1.0)
        .def("add_move", py::overload_cast<const std::string &, double, double, double, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"), py::arg("prob") = 1.0, py::arg("dx") = 0.1, py::arg("Ddt") = 0.1, py::arg("box_id") = -1)
        .def("add_move", py::overload_cast<const std::string &, double, const std::string &, double, double, bool, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"), py::arg("prob") = 1.0, py::arg("particle"), py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0, py::arg("energy_bias") = false, py::arg("box_id") = -1)
        .def("add_move", py::overload_cast<const std::string &, double, std::vector<Particle>, double, double, double, bool, bool, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"), py::arg("prob") = 1.0, py::arg("molecule"), py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0, py::arg("r_inner") = 1.3, py::arg("energy_bias") = false, py::arg("target_mol") = false, py::arg("box_id") = -1)
        .def("add_box", py::overload_cast<>(&System::add_box))
        .def("add_box", py::overload_cast<Box *>(&System::add_box))
        .def("set_forcefield", py::overload_cast<ForceField *, int>(&System::set_forcefield),
            py::arg("forcefield"), py::arg("box_id") = -1)
        .def("set_forcefield", py::overload_cast<const std::string &, const std::vector<std::string> &, int>(&System::set_forcefield),
            py::arg("forcefield"), py::arg("elements"), py::arg("box_id") = -1)
        .def("set_forcefield", py::overload_cast<const std::string &, const std::string &, int>(&System::set_forcefield),
            py::arg("forcefield"), py::arg("paramfile"), py::arg("box_id") = -1)
        .def("set_boundary", py::overload_cast<Boundary *, int>(&System::set_boundary),
            py::arg("boundary"), py::arg("box_id") = -1)
        .def("set_boundary", py::overload_cast<const std::string &, std::valarray<double>, int>(&System::set_boundary),
            py::arg("boundary"), py::arg("length") = {}, py::arg("box_id") = -1)
        .def("add_constraint", &System::add_constraint, py::arg("constraint"), py::arg("box_id") = -1)
        .def("snapshot", &System::snapshot, py::arg("filename"), py::arg("box_id") = 0)
        .def("set_dump", &System::set_dump, py::arg("freq") = 1, py::arg("filename") = "mc.xyz", py::arg("outputs") = {}, py::arg("box_id") = 0)
        .def("set_thermo", &System::set_thermo, py::arg("freq") = 1, py::arg("filename") = "mc.log", py::arg("outputs") = {}, py::arg("box_id") = 0)
        .def("add_particle", py::overload_cast<Particle, int>(&System::add_particle), "Add a particle by object", py::arg("particle"), py::arg("box_id") = 0)
        .def("add_particle", py::overload_cast<const std::string &, std::valarray<double>, int >(&System::add_particle),
            "Add a particle by element and position", py::arg("element"), py::arg("position"), py::arg("box_id") = 0)
        .def("add_particles", py::overload_cast<std::vector<Particle>, int>(&System::add_particles), "Add particles by a list of objects",
            py::arg("particles"), py::arg("box_id") = 0)
        .def("add_particles", py::overload_cast<const std::string &, std::vector<std::valarray<double> >, int>(&System::add_particles),
            "Add particles of the same type", py::arg("element"), py::arg("positions"), py::arg("box_id") = 0)
        .def("read_particles", &System::read_particles, "Read particles from xyz-file", py::arg("filename"), py::arg("box_id") = 0)

        //.def("run_md", &System::run_md, py::arg("steps"))
        .def("run_mc", &System::run_mc, py::arg("cycles"), py::arg("steps") = 1)
        .def_readonly("nbox", &System::nbox)
        .def_readonly("ndim", &System::ndim)
        .def_readonly("temp", &System::temp)
        .def_readonly("chempot", &System::chempot)
        .def_readonly("working_dir", &System::working_dir);

    // Box
    py::class_<Box>(m, "Box")
        .def(py::init<System *>()) //, int>()) //, py::arg("memory_intensity") = 2)
        .def("set_forcefield", &Box::set_forcefield)
        .def("set_boundary", &Box::set_boundary)
        //.def("set_integrator", &Box::set_integrator)
        .def("add_particle", py::overload_cast<Particle>(&Box::add_particle))
        .def("add_particle", py::overload_cast<const std::string &, std::valarray<double> >(&Box::add_particle))
        .def("add_particles", py::overload_cast<std::vector<Particle> >(&Box::add_particles))
        .def("add_particles", py::overload_cast<const std::string &, std::vector<std::valarray<double> > >(&Box::add_particles))
        .def("read_particles", &Box::read_particles)
        .def("add_constraint", &Box::add_constraint)
        .def("snapshot", &Box::snapshot)
        .def("set_dump", &Box::set_dump)
        .def("set_thermo", &Box::set_thermo)
        .def_readonly("npar", &Box::npar)
        .def_readonly("step", &Box::step)
        .def_readonly("box_id", &Box::box_id)
        .def_readonly("nconstraint", &Box::nconstraint);

    // Sampler
    py::class_<Sampler>(m, "Sampler")
        .def("sample", &Sampler::sample)
        .def_readonly("label", &Sampler::label);
    py::class_<Metropolis, Sampler>(m, "Metropolis")
        .def(py::init<System *>());
    py::class_<Umbrella, Sampler>(m, "Umbrella")
        .def(py::init<System *, std::function<double (int)>, int>());

    // Moves
    py::class_<Moves>(m, "Moves")
        .def_readonly("label", &Moves::label)
        .def_readonly("ndrawn", &Moves::ndrawn)
        .def_readonly("naccept", &Moves::naccept)
        .def_readonly("cum_time", &Moves::cum_time);
    py::class_<Trans, Moves>(m, "Trans")
        .def(py::init<System *, Box *, double>());
    py::class_<TransMH, Moves>(m, "TransMH")
        .def(py::init<System *, Box *, double, double>());

    // ForceField
    py::class_<ForceField>(m, "ForceField")
        .def_readonly("label2type", &ForceField::label2type)
        .def_readonly("unique_labels", &ForceField::unique_labels)
        .def_readonly("ntype", &ForceField::ntype)
        .def_readonly("temp_scale", &ForceField::temp_scale)
        .def_readonly("label", &ForceField::label)
        .def_readonly("paramfile", &ForceField::paramfile);
    py::class_<LennardJones, ForceField>(m, "LennardJones")
        .def(py::init<Box *, std::string>());
    py::class_<Vashishta, ForceField>(m, "Vashishta")
        .def(py::init<Box *, std::string>());
    py::class_<IdealGas, ForceField>(m, "IdealGas")
        .def(py::init<Box *, std::vector<std::string> >());

    // Boundary
    py::class_<Boundary>(m, "Boundary")
        .def_readwrite("label", &Boundary::label);
    py::class_<Open, Boundary>(m, "Open")
        .def(py::init<Box *>());
    py::class_<Periodic, Boundary>(m, "Periodic")
        .def(py::init<Box *, std::valarray<double> >());

    // Constraint
    py::class_<Constraint>(m, "Constraint")
        .def_readwrite("label", &Constraint::label);
    py::class_<Stillinger, Constraint>(m, "Stillinger")
        .def(py::init<Box *, double>())
        .def("set_criterion", &Stillinger::set_criterion);
    py::class_<MinNeigh, Constraint>(m, "MinNeigh")
        .def(py::init<Box *, std::string, std::string, double, int>());
    py::class_<MaxNeigh, Constraint>(m, "MaxNeigh")
        .def(py::init<Box *, std::string, std::string, double, int>());
    py::class_<MinDistance, Constraint>(m, "MinDistance")
        .def(py::init<Box *, std::string, std::string, double>());
    py::class_<MaxDistance, Constraint>(m, "MaxDistance")
        .def(py::init<Box *, std::string, std::string, double>());

    // Integrator
    py::class_<Integrator>(m, "Integrator")
        .def_readwrite("dt", &Integrator::dt);
    py::class_<Euler, Integrator>(m, "Euler")
        .def(py::init<Box *, double>());
    py::class_<EulerCromer, Integrator>(m, "EulerCromer")
        .def(py::init<Box *, double>());
    py::class_<VelocityVerlet, Integrator>(m, "VelocityVerlet")
        .def(py::init<Box *, double>());
    py::class_<RungeKutta4, Integrator>(m, "RungeKutta4")
        .def(py::init<Box *, double>());

    // Thermo
    py::class_<Thermo>(m, "Thermo")
        .def(py::init<Box *, int, std::string, std::vector<std::string> >());

    // Dump
    py::class_<Dump>(m, "Dump")
        .def(py::init<Box *, int, std::string, std::vector<std::string> >());
}

int main()
{
    std::cout << "hey" << std::endl;
    return 0;
}
