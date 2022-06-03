/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

#include <iostream>

#include "boundary.h"
#include "../box.h"


/* ----------------------------------------------------------------------------
   Boundary base class constructor, overwriting box pointer
------------------------------------------------------------------------------- */

Boundary::Boundary(Box* box_in)
{
    box = box_in;
}
