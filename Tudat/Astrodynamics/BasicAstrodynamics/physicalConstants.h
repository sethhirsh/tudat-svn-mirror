/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      131207    S. Hirsh         	File Created.
 *
 *    References
 *      NIST. NIST reference on constants, units and uncertainty.
 *          http://physics.nist.gov/cgi-bin/cuu/Value?mu0, last accessed: 7th December, 2013.
 *
 *    Notes
 *
 */

#ifndef TUDAT_PHYSICAL_CONSTANTS_H
#define TUDAT_PHYSICAL_CONSTANTS_H

#include <cmath>

#include "TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h"


namespace tudat
{
namespace basic_astrodynamics
{
namespace physical_constants
{

//! Magnetic constant
/*!
 * Magnetic constant in V*s/(A*m) (NIST, 2013).
 */
const static double MAGNETIC_CONSTANT = tudat::basic_mathematics::mathematical_constants::PI * 4.0e-7;


} // namespace physical_constants
} // namespace basic_astrodynamics
} // namespace tudat


#endif // TUDAT_PHYSICAL_CONSTANTS_H
