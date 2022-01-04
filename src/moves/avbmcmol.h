#pragma once
#include <string>
#include <memory>

#include "avbmcinmol.h"
#include "avbmcoutmol.h"


class AVBMCMol : public AVBMCInMol, public AVBMCOutMol
{
public:
    //AVBMCMol(class System *, class Box *, double = 0.95, double = 3.0);
    AVBMCMol(class System *, std::shared_ptr<class Box>, double = 0.95, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();
    std::string repr();

private:
    double r_above, r_below;
    bool move_in;
    //class Box* box = nullptr;
    std::shared_ptr<class Box> box = nullptr;
};
