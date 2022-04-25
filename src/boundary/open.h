#pragma once
#include <valarray>

#include "boundary.h"


class Open : public Boundary
{
public:
    Open(class Box *);
};
