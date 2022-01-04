#include <iostream>
#include <memory>

#include "boundary.h"
#include "../box.h"

/*
Boundary::Boundary(Box* box_in)
{
    box = box_in;
}
*/

Boundary::Boundary(std::shared_ptr<Box> box_in)
{
    box = box_in;
}
