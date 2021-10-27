#include "parser.h"


void parser(int argc, char** argv[])
{
    /* Parser main function
     */
    if(argc >= 2){
        Box box();
        string cmd_cat = argv[1];
        if(cmd_cat == "set"){
            set(box, argc, argv[]);
        }
        else if(cmd_cat == "add"){
            add(box, argc, argv[]);
        }
        else if(cmd_cat == "take"){
            take(box, argc, argv[]);
        }
        else if(cmd_cat == "run"){
            run(box, argc, argv[]);
        }
        else{
            cout << "There is no category '" + cmd_cat + "'! Aborting." << endl;
            exit(0);
        }
    }
    else{
        cout << "avbmc version 0.0.1" << endl;
        cout << "Author: Even M. Nordhagen" << endl;
        cout << "url: github.com/evenmn/avbmc" << endl;
    }
}


void set(Box& box, int argc, char** argv[])
{
    /*
     */
    
    assert(argc > 2);
    string keyword = argv[2];
    if(keyword == "temp"){
        box.set_temp(stod(argv[3]));
    }
    else if(keyword == "chempot"){
        box.set_chempot(stod(argv[3]));
    }
    else if(keyword == "forcefield"){
        assert(argc >= 5);
        string forcefield = argv[3];
        string paramfile = argv[4];
        if(forcefield == "lennardjones"){
            box.set_forcefield(new LennardJones(&box, paramfile));
        }
        else if(forcefield == "vashishta"){
            box.set_forcefield(new Vashishta(&box, paramfile));
        }
        else{
            cout << "Forcefield '" + forcefield + "' is not known! Aborting." << end;
            exit(0);
        }
    }
    else if(keyword == "integrator"){
        string integrator = argv[3];
        if(argc == 4){
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box));
            }
            else{
                cout << "Integrator '" + integrator + "' is not known! Aborting." << endl;
                exit(0);
            }
        else if(argc >= 5){
            double dt = stod(argv[4]);
            if(integrator == "euler"){
                box.set_integrator(new Euler(&box, dt));
            }
            else if(integrator == "eulercromer"){
                box.set_integrator(new EulerCromer(&box, dt));
            }
            else if(integrator == "velocityverlet"){
                box.set_integrator(new VelocityVerlet(&box, dt));
            }
            else if(integrator == "rungekutta4"){
                box.set_integrator(new RungeKutta4(&box, dt));
            }
            else{
                cout << "Integrator '" + integrator + "' is not known! Aborting." << endl;
                exit(0);
            }
        }
        else{
            cout << "Integrator does not have a sufficient number of arguments! Aborting." << endl;
            exit(0);
        }
    }
    else{
        cout << "Keyword '" + keyword + "' is not known! Aborting." << endl;
        exit(0);
    }
}
