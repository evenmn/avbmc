#include "moves.h"
#include "../box.h"

Moves::Moves(Box* box_in)
{
    box = box_in;
    rng = box->rng;
}
