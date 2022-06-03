/* ----------------------------------------------------------------------------
  This file is a part of the AVBMC library, which follows the GPL-3.0 License.
  For license information, see LICENSE file in the top directory, 
  https://github.com/evenmn/avbmc/LICENSE.

  Author(s): Even M. Nordhagen
  Email(s): evenmn@mn.uio.no
  Date: 2022-06-03 (last changed 2022-06-03)
---------------------------------------------------------------------------- */

#include <iostream>
#include <string>
#include <valarray>

#include "particle.h"

/* ----------------------------------------------------------------------------
   Particle constructor. Takes particle label 'label_in' and particle initial
   position 'r_in' as arguments.
-------------------------------------------------------------------------------
*/

Particle::Particle(const std::string label_in, const std::valarray<double> r_in)
{
    label = label_in;
    r = r_in;
}


/* ----------------------------------------------------------------------------
   Particle copy constructor.
------------------------------------------------------------------------------- */

Particle::Particle(const Particle &other)
{
    label = other.label;
    r = other.r;
    v = other.v;
    f = other.f;
    type = other.type;
    mass = other.mass;
    poteng = other.poteng;
}
