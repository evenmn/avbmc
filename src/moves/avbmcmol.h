#pragma once
#include "avbmcinmol.h"
#include "avbmcoutmol.h"


class AVBMCMol : public AVBMCInMol, public AVBMCOutMol
{
public:
    AVBMCMol(class Box*, double = 0.95, double = 3.0);
    void perform_move();
    double accept(double, double);
    void reset();

private:
    bool move_in;
};
