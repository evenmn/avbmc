#include "moves.h"
#include "../box.h"

Moves::Moves(class Box* box_in)
{
    box = box_in;
    rng = box->rng;
}
