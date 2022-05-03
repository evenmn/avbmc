# AVBMC Code Style
The AVBMC C++ source code has to meet some requirements to be neet and readable.
Usually, every cpp source file contains a class, which might be a subclass. The
opposite is also true: All new classes are written to new files. 

Every cpp source file has to contain the following:

1. A doc string with the standard licence information, author and date last edited
2. A doc string about the tasks of the file/class
3. External includes (standard library, for instance)
4. Local includes (usually from AVBMC header files)
5. A doc string before every function, including the constructor and destructor

## Doc string rules
Every doc string starts with a separation line constructed with `/*` + 77 * `-` 
to make it 80 characters long:
``` c++
/* ----------------------------------------------------------------------------
```
Thereafter, the text is written. A text line should never exceed the length of
the separation line (80 characters).

The doc string is ended with another separation line constructed with 80 * `-`
+ `*/`.

## Spacing rules
There is one line gap from doc strings to code and two lines gap from code to
doc strings. 

## Code rules
We follow the clang code format when writing code. This can be done by running
`git-clang-format` before committing. Function arguments have always names 
ending with '`_in`', while class variables do not have any particular ending.
Code lines should not exceed 80 characters.

## An example
``` c++
/* ----------------------------------------------------------------------------
   Licence information ...

   Author(s): Even M. Nordhagen
   Email(s): evenmn@mn.uio.no
   Date: 2022-04-24 (last changed 2022-04-25)
------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
   This is the particle base class, used to create particle objects. A particle
   object consists of a label and position at the minimum. Might also contain
   information about mass (this is only needed by molecular dynamics).
------------------------------------------------------------------------------- */

#include <iostream>
#include <valarray>
#include <string>

#include "particle.h"


/* ----------------------------------------------------------------------------
   Particle constructor, taking in the particle label, 'label_in', and the
   initial particle position, 'position_in', as arguments.
------------------------------------------------------------------------------- */

Particle::Particle(const std::string &label_in, std::valarray<double> position_in)
{
    label = label_in;
    position = position_in;
}


/* ----------------------------------------------------------------------------
   Assign particle mass, 'mass_in'. Warning: Mass is overwritten every time
   this method is called, which is usually not the desired behavior.
------------------------------------------------------------------------------- */

void Particle::assign_mass(double mass_in)
{
    mass = mass_in;
}
```
