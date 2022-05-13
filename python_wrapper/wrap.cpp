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
#include "../src/init_position.h"

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
#include "../src/moves/avbmcswapright.h"

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

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
using namespace pybind11::literals;


/* ----------------------------------------------------------------------------
   Python wrapper based on pybind11. Please see the pybind11 documentation
   for all function calls starting with "py::".
------------------------------------------------------------------------------- */

PYBIND11_MODULE(avbmc, m) {

    // Particle
    py::class_<Particle>(m, "Particle", py::dynamic_attr())
        .def(py::init<std::string, std::valarray<double> >())
        .def("__repr__",
            [](const Particle &particle) {
                return particle.label;
            }
        )
        .def_readonly("label", &Particle::label)
        .def_readwrite("r", &Particle::r)
        .def_readwrite("v", &Particle::v)
        .def_readwrite("f", &Particle::f);

    // Initialize position
    m.def("read_xyz", &read_xyz,
        "Read an XYZ-file and return list of Particle objects",
        py::arg("filename")
    );
    m.def("fcc", &fcc,
        "Initialize in face-centered cube",
        py::arg("n"),
        py::arg("L"),
        py::arg("ndim") = 3
    );

    // Random Number Generator
    py::class_<RandomNumberGenerator>(m, "RNG")
        .def("__repr__",
            [](const RandomNumberGenerator &rng) {
                return rng.label;
            }
        )
        .def_readonly("label", &RandomNumberGenerator::label);
    py::class_<MersenneTwister, RandomNumberGenerator>(m, "MersenneTwister")
        .def(py::init<>())
        .def("set_seed",
            &MersenneTwister::set_seed,
            "Set seed to be used by the random number generator",
            py::arg("seed")
        )
        .def("next_int",
            &MersenneTwister::next_int,
            "Get next int in the sequence",
            py::arg("upper_limit")
        )
        .def("next_double",
            &MersenneTwister::next_double,
            "Get next double (between 0 and 1) in the sequence"
        )
        .def("next_gaussian",
            &MersenneTwister::next_gaussian,
            "Get next gaussian in the sequence",
            py::arg("mean") = 0.,
            py::arg("var") = 1.
        )
        .def("choice",
            &MersenneTwister::choice,
            "Index is picked according to a list of probabilities",
            py::arg("probabilities")
        );

    // System
    py::class_<System>(m, "System")
        .def(py::init<std::string, bool>(), 
            "System base class constructor",
            py::arg("working_dir") = "",
            py::arg("initialize") = true
        )
        .def("set_sampler",
            py::overload_cast<Sampler *>(&System::set_sampler),
            "Set sampler used in Monte Carlo simulations",
            py::arg("sampler_object")
        )
        .def("set_sampler",
            py::overload_cast<const std::string &, std::function<double(int)> >(&System::set_sampler),
            "Set sampler used in Monte Carlo simulations",
            py::arg("sampler_label"),
            py::arg("weight_function") = py::none()
        )
        .def("set_rng",
            py::overload_cast<RandomNumberGenerator *>(&System::set_rng),
            "Set random number generator used by system",
            py::arg("rng_object")
        )
        .def("set_rng",
            py::overload_cast<const std::string &>(&System::set_rng),
            "Set random number generator used by system",
            py::arg("rng_label")
        )
        .def("set_temp",
            &System::set_temp,
            "Set system temperature to be used in canonical and grand-canonical Monte Carlo simulations",
            py::arg("temperature")
        )
        .def("set_seed",
            &System::set_seed,
            "Set system seed to be used by random number generator",
            py::arg("seed")
        )
        .def("set_chempot",
            &System::set_chempot,
            "Set chemical potential to be used in grand-canonical Monte Carlo simulations",
            py::arg("chempot")
        )
        .def("add_move",
            py::overload_cast<Moves *, double>(&System::add_move),
            "Add move to the list of moves",
            py::arg("system"),
            py::arg("prob") = 1.0
        )
        .def("add_move",
            py::overload_cast<const std::string &, double, double, double, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"),
            py::arg("prob") = 1.0,
            py::arg("dx") = 0.1,
            py::arg("Ddt") = 0.1,
            py::arg("box_id") = -1
        )
        .def("add_move",
            py::overload_cast<const std::string &, double, const std::string &, double, double, bool, int, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"),
            py::arg("prob") = 1.0,
            py::arg("particle"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("energy_bias") = false,
            py::arg("box_id") = 0,
            py::arg("box_id2") = 1
        )
        .def("add_move",
            py::overload_cast<const std::string &, double, std::vector<Particle>, double, double, double, bool, bool, int, int>(&System::add_move),
            "Add move to the list of moves",
            py::arg("move"),
            py::arg("prob") = 1.0,
            py::arg("molecule"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("r_inner") = 1.3,
            py::arg("energy_bias") = false,
            py::arg("target_mol") = false,
            py::arg("box_id") = 0,
            py::arg("box_id2") = 1
        )
        .def("add_box",
            py::overload_cast<>(&System::add_box),
            "Add system box automatically"
        )
        .def("add_box",
            py::overload_cast<Box *>(&System::add_box),
            "Add system box manually with a box object",
            py::arg("box_object")
        )
        .def("set_forcefield",
            py::overload_cast<ForceField *, int>(&System::set_forcefield),
            "Set forcefield of box by forcefield object. All boxes by default",
            py::arg("forcefield_object"),
            py::arg("box_id") = -1
        )
        .def("set_forcefield",
            py::overload_cast<const std::string &, const std::vector<std::string> &, int>(&System::set_forcefield),
            "Set forcefield of box by forcefield label. Used for ideal gas only. All boxes by default",
            py::arg("forcefield_label"),
            py::arg("elements"),
            py::arg("box_id") = -1
        )
        .def("set_forcefield",
            py::overload_cast<const std::string &, const std::string &, int>(&System::set_forcefield),
            "Set forcefield of box by forcefield label and parameter file. All boxes by default",
            py::arg("forcefield_label"),
            py::arg("paramfile"),
            py::arg("box_id") = -1
        )
        .def("set_boundary",
            py::overload_cast<Boundary *, int>(&System::set_boundary),
            "Set boundary of box by boundary object. All boxes by default",
            py::arg("boundary_object"),
            py::arg("box_id") = -1
        )
        .def("set_boundary",
            py::overload_cast<const std::string &, std::valarray<double>, int>(&System::set_boundary),
            "Set boundary of box by boundary label. All boxes by default",
            py::arg("boundary_label"),
            py::arg("length") = py::none(),
            py::arg("box_id") = -1
        )
        .def("add_constraint",
            py::overload_cast<Constraint *, int>(&System::add_constraint),
            "Add box constraint by constraint object. box_id = 0 by default",
            py::arg("constraint_object"),
            py::arg("box_id") = 0
        )
        .def("add_constraint",
            py::overload_cast<const std::string &, const std::string &, const std::string &, double, int, int>(&System::add_constraint),
            "Add box constraint by constraint label. box_id = 0 by default",
            py::arg("constraint"),
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance"),
            py::arg("nneigh") = 1,
            py::arg("box_id") = 0
        )
        .def("snapshot",
            &System::snapshot,
            "Take snapshot of box and write to xyz-file. box_id = 0 by default",
            py::arg("filename"),
            py::arg("box_id") = 0
        )
        .def("set_dump",
            &System::set_dump, 
            "Set particle output to be written to xyz-file. box_id = 0 by default",
            py::arg("freq") = 1,
            py::arg("filename") = "mc.xyz",
            py::arg("outputs") = py::none(),
            py::arg("box_id") = 0
        )
        .def("set_thermo",
            &System::set_thermo,
            "Set system output to be written to txt-file. box_id = 0 by default",
            py::arg("freq") = 1,
            py::arg("filename") = "mc.log",
            py::arg("outputs") = py::none(),
            py::arg("box_id") = 0
        )
        .def("add_particle",
            py::overload_cast<Particle, int>(&System::add_particle),
            "Add a particle by object. box_id = 0 by default",
            py::arg("particle_object"),
            py::arg("box_id") = 0
        )
        .def("add_particle",
            py::overload_cast<const std::string &, std::valarray<double>, int >(&System::add_particle),
            "Add a particle by element and position. box_id = 0 by default",
            py::arg("element"),
            py::arg("position"),
            py::arg("box_id") = 0
        )
        .def("add_particles",
            py::overload_cast<std::vector<Particle>, int>(&System::add_particles),
            "Add particles by a list of objects. box_id = 0 by default",
            py::arg("particle_objects"),
            py::arg("box_id") = 0
        )
        .def("add_particles",
            py::overload_cast<const std::string &, std::vector<std::valarray<double> >, int>(&System::add_particles),
            "Add a list of particles by position array. They have to be of the same type. box_id = 0 by default",
            py::arg("element"),
            py::arg("positions"),
            py::arg("box_id") = 0
        )
        .def("read_particles",
            &System::read_particles,
            "Read and add particles from an xyz-file. box_id = 0 by default",
            py::arg("filename"),
            py::arg("box_id") = 0
        )
        .def("get_size_histogram",
            &System::get_size_histogram,
            "Get histogram of all system sizes during simulation. box_id = 0 by default",
            py::arg("box_id") = 0
        )
        .def("write_size_histogram",
            &System::write_size_histogram,
            "Write size histogram out to txt-file. box_id = 0 by default",
            py::arg("filename"),
            py::arg("box_id") = 0
        )
        //.def("run_md",
        //  &System::run_md, 
        //  "Perform molecular dynamics run",
        //  py::arg("steps")
        //)
        .def("run_mc",
            &System::run_mc,
            "Perform Monte Carlo run",
            py::arg("cycles"),
            py::arg("steps") = 1
        )
        .def("run_mc_cycle",
            &System::run_mc_cycle,
            "Run a Monte Carlo cycle",
            py::arg("steps") = 1
        )
        .def("initialize_mc_run",
            &System::initialize_mc_run,
            "Initialize Monte Carlo run before running"
        )
        .def_readonly("nbox", &System::nbox)
        .def_readonly("ndim", &System::ndim)
        .def_readonly("temp", &System::temp)
        .def_readonly("chempot", &System::chempot)
        .def_readonly("working_dir", &System::working_dir)
        .def_readonly("boxes", &System::boxes)
        .def_readonly("moves", &System::moves);

    // Box
    py::class_<Box>(m, "Box")
        .def(py::init<System *>(), //, int>(),
            "Box class constructor",
            py::arg("system_object")
            //, py::arg("memory_intensity") = 2
        )
        .def("set_forcefield",
            &Box::set_forcefield,
            "Set box forcefield by forcefield object",
            py::arg("forcefield_object")
        )
        .def("set_boundary",
            &Box::set_boundary,
            "Set box boundary by boundary object",
            py::arg("boundary_object")
        )
        //.def("set_integrator",
        //  &Box::set_integrator,
        //  "Set box integrator by integrator object",
        //  py::arg("integrator_object")
        //)
        .def("add_particle",
            py::overload_cast<Particle>(&Box::add_particle),
            "Add a particle to the box py particle object",
            py::arg("particle_object")
        )
        .def("add_particle",
            py::overload_cast<const std::string &, std::valarray<double> >(&Box::add_particle),
            "Add a particle to the box by particle label and position",
            py::arg("particle_label"),
            py::arg("position")
        )
        .def("add_particles",
            py::overload_cast<std::vector<Particle> >(&Box::add_particles),
            "Add a list of particles to the box by particle objects",
            py::arg("particle_objects")
        )
        .def("add_particles",
            py::overload_cast<const std::string &, std::vector<std::valarray<double> > >(&Box::add_particles),
            "Add a list of particles to the box by particle label and positions",
            py::arg("particle_label"),
            py::arg("positions")
        )
        .def("read_particles",
            &Box::read_particles,
            "Add a set of particles to the box read from an xyz-file",
            py::arg("filename")
        )
        .def("add_constraint",
            &Box::add_constraint,
            "Add box constraint by object",
            py::arg("constraint_object")
        )
        .def("snapshot",
            &Box::snapshot,
            "Take snapshot of box and write to xyz-file",
            py::arg("filename")
        )
        .def("set_dump",
            &Box::set_dump,
            "Set particle output written to xyz-file",
            py::arg("freq"),
            py::arg("filename"),
            py::arg("outputs")
        )
        .def("set_thermo",
            &Box::set_thermo,
            "Set box output written to txt-file",
            py::arg("freq"),
            py::arg("filename"),
            py::arg("outputs")
        )
        .def_readonly("npar", &Box::npar)
        .def_readonly("step", &Box::step)
        .def_readonly("box_id", &Box::box_id)
        .def_readonly("nconstraint", &Box::nconstraint);

    // Sampler
    py::class_<Sampler>(m, "Sampler")
        .def("sample",
            &Sampler::sample,
            "Perform sampling moves",
            py::arg("moves") = 1
        )
        .def("__repr__",
            [](const Sampler &sampler) {
                return sampler.label;
            }
        )
        .def_readonly("label", &Sampler::label);
    py::class_<Metropolis, Sampler>(m, "Metropolis")
        .def(py::init<System *>(),
            "Metropolis class constructor",
            py::arg("system_object")
        );
    py::class_<Umbrella, Sampler>(m, "Umbrella")
        .def(py::init<System *, std::function<double (int)>, int>(),
            "Umbrella sampling class constructor",
            py::arg("system_object"),
            py::arg("weight_function"),
            py::arg("ntabulated") = 100
        );

    // Moves
    py::class_<Moves>(m, "Moves")
        .def("__repr__",
            [](const Moves &move) {
                return move.label;
            }
        )
        .def_readonly("label", &Moves::label)
        .def_readonly("ndrawn", &Moves::ndrawn)
        .def_readonly("naccept", &Moves::naccept)
        .def_readonly("cum_time", &Moves::cum_time);
    py::class_<Trans, Moves>(m, "Trans")
        .def(py::init<System *, Box *, double>(),
            "Naive translational move class constructor",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("dx") = 0.1
        );
    py::class_<TransMH, Moves>(m, "TransMH")
        .def(py::init<System *, Box *, double, double>(),
            "Biased translational move class constructor",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("dx") = 0.1,
            py::arg("Ddt") = 0.1
        );
    py::class_<AVBMC, Moves>(m, "AVBMC")
        .def(py::init<System *, Box *, const std::string &, double, double, bool>(),
            "Aggregate volume-biased insertation and deletion moves",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("particle"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("energy_bias") = false
        );
    py::class_<AVBMCIn, Moves>(m, "AVBMCIn")
        .def(py::init<System *, Box *, const std::string &, double, double, bool>(),
            "Aggregate volume-biased insertation move",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("particle"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("energy_bias") = false
        );
    py::class_<AVBMCOut, Moves>(m, "AVBMCOut")
        .def(py::init<System *, Box *, const std::string &, double, bool>(),
            "Aggregate volume-biased deletion move",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("particle"),
            py::arg("r_above") = 3.0,
            py::arg("energy_bias") = false
        );
    py::class_<AVBMCMol, Moves>(m, "AVBMCMol")
        .def(py::init<System *, Box *, std::vector<Particle>, double, double, double, bool, bool>(),
            "Aggregate volume-biased molecule insertation and deletion moves",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("molecule"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("r_inner") = 1.3,
            py::arg("energy_bias") = false,
            py::arg("target_molecule") = false
        );
    py::class_<AVBMCMolIn, Moves>(m, "AVBMCMolIn")
        .def(py::init<System *, Box *, std::vector<Particle>, double, double, double, bool, bool>(),
            "Aggregate volume-biased molecule insertation move",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("molecule"),
            py::arg("r_below") = 0.95,
            py::arg("r_above") = 3.0,
            py::arg("r_inner") = 1.3,
            py::arg("energy_bias") = false,
            py::arg("target_molecule") = false
        );
    py::class_<AVBMCMolOut, Moves>(m, "AVBMCMolOut")
        .def(py::init<System *, Box *, std::vector<Particle>, double, double, bool, bool>(),
            "Aggregate volume-biased molecule insertation move",
            py::arg("system_object"),
            py::arg("box_object"),
            py::arg("molecule"),
            py::arg("r_above") = 3.0,
            py::arg("r_inner") = 1.3,
            py::arg("energy_bias") = false,
            py::arg("target_molecule") = false
        );

    // ForceField
    py::class_<ForceField>(m, "ForceField")
        .def("__repr__",
            [](const ForceField &forcefield) {
                return forcefield.label;
            }
        )
        .def_readonly("label2type", &ForceField::label2type)
        .def_readonly("unique_labels", &ForceField::unique_labels)
        .def_readonly("ntype", &ForceField::ntype)
        .def_readonly("temp_scale", &ForceField::temp_scale)
        .def_readonly("label", &ForceField::label)
        .def_readonly("paramfile", &ForceField::paramfile);
    py::class_<LennardJones, ForceField>(m, "LennardJones")
        .def(py::init<Box *, std::string>(),
            "Lennard-Jones class constructor",
            py::arg("box_object"),
            py::arg("parameter_file")
        );
    py::class_<Vashishta, ForceField>(m, "Vashishta")
        .def(py::init<Box *, std::string>(),
            "Vashishta class constructor",
            py::arg("box_object"),
            py::arg("parameter_file")
        );
    py::class_<IdealGas, ForceField>(m, "IdealGas")
        .def(py::init<Box *, std::vector<std::string> >(),
            "Ideal gas class constructor",
            py::arg("box_object"),
            py::arg("elements")
        );

    // Boundary
    py::class_<Boundary>(m, "Boundary")
        .def("__repr__",
            [](const Boundary &boundary) {
                return boundary.label;
            }
        )
        .def_readonly("label", &Boundary::label);
    py::class_<Open, Boundary>(m, "Open")
        .def(py::init<Box *>(),
            "Open class constructor",
            py::arg("box_object")
        );
    py::class_<Periodic, Boundary>(m, "Periodic")
        .def(py::init<Box *, std::valarray<double> >(),
            "Periodic class constructor",
            py::arg("box_object"),
            py::arg("length")
        );

    // Constraint
    py::class_<Constraint>(m, "Constraint")
        .def("__repr__",
            [](const Constraint &constraint) {
                return constraint.label;
            }
        )
        .def_readonly("label", &Constraint::label);
    py::class_<Stillinger, Constraint>(m, "Stillinger")
        .def(py::init<Box *, double>(),
            "Stillinger class constructor",
            py::arg("box"),
            py::arg("distance") = 1.0
        )
        .def("set_criterion",
            &Stillinger::set_criterion,
            "Set distance criterion between two elements",
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance")
        );
    py::class_<MinNeigh, Constraint>(m, "MinNeigh")
        .def(py::init<Box *, std::string, std::string, double, int>(),
            "Minimum number of neighbors class constructor",
            py::arg("box_object"),
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance"),
            py::arg("nneigh")
        );
    py::class_<MaxNeigh, Constraint>(m, "MaxNeigh")
        .def(py::init<Box *, std::string, std::string, double, int>(),
            "Maximum number of neighbors class constructor",
            py::arg("box_object"),
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance"),
            py::arg("nneigh")
        );
    py::class_<MinDistance, Constraint>(m, "MinDistance")
        .def(py::init<Box *, std::string, std::string, double>(),
            "Minimum distance class constructor",
            py::arg("box_object"),
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance")
        );
    py::class_<MaxDistance, Constraint>(m, "MaxDistance")
        .def(py::init<Box *, std::string, std::string, double>(),
            "Maximum distance class constructor",
            py::arg("box_object"),
            py::arg("element1"),
            py::arg("element2"),
            py::arg("distance")
        );

    // Integrator
    py::class_<Integrator>(m, "Integrator")
        .def_readonly("dt", &Integrator::dt);
    py::class_<Euler, Integrator>(m, "Euler")
        .def(py::init<Box *, double>(),
            "Euler class constructor",
            py::arg("box_object"),
            py::arg("dt") = 0.001
        );
    py::class_<EulerCromer, Integrator>(m, "EulerCromer")
        .def(py::init<Box *, double>(),
            "Euler-Cromer class constructor",
            py::arg("box_object"),
            py::arg("dt") = 0.01
        );
    py::class_<VelocityVerlet, Integrator>(m, "VelocityVerlet")
        .def(py::init<Box *, double>(),
            "Velocity verlet class constructor",
            py::arg("box_object"),
            py::arg("dt") = 0.01
        );
    py::class_<RungeKutta4, Integrator>(m, "RungeKutta4")
        .def(py::init<Box *, double>(),
            py::arg("box_object"),
            py::arg("dt") = 0.01
        );

    // Thermo
    py::class_<Thermo>(m, "Thermo")
        .def(py::init<Box *, int, std::string, std::vector<std::string> >(),
            "Thermo class constructor",
            py::arg("box_object"),
            py::arg("freq"),
            py::arg("filename"),
            py::arg("outputs")
        );

    // Dump
    py::class_<Dump>(m, "Dump")
        .def(py::init<Box *, int, std::string, std::vector<std::string> >(),
            "Dump class constructor",
            py::arg("box_object"),
            py::arg("freq"),
            py::arg("filename"),
            py::arg("outputs")
        );
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
