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
 *      131109	  S. Hirsh			File Created. 
 *
 *    References
 *		
 *
 *    Notes
 *
 */

#ifndef TUDAT_LOCAL_MAGNETIC_DIPOLE_FIELD_H
#define TUDAT_LOCAL_MAGNETIC_DIPOLE_FIELD_H

#include <Eigen/Core>
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute local magnetic field using static dipole model
/*!
 * Computes local magnetic field using static dipole model in an inertial reference frame.	
 * \param positionOfBodySubjectToAcceleration position of body which is being
 *			accelerated by the Lorentz force 													[m]
 * \param originOfMagneticField origin of magnetic field of body that is creating the
 *			magnetic dipole moment 																[m]
 * \param axisOfMagneticField unitvector that points in the direction of the symmetry axis of
 *			magnetic dipole field.																[-]
 * \param dipoleMomentOfMagneticField magnetic dipole moment from body exerting acceleration 	[T·m]
 * \param permeabilityOfSpaceBetweenBodies permeability of space between bodies. The default for 
 *			this value is the magnetic constant which is defined as 4*pi*10^-7 V·s/(A·m)  [V·s/(A·m)]
 * \return local magnetic field using static dipole model			                            [T·m]
 */

//! Compute local magnetic field using static dipole model.
Eigen::Vector3d computeLocalStaticMagneticDipoleField(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const Eigen::Vector3d& originOfMagneticField,
        const Eigen::Vector3d& axisOfMagneticField,
        const double dipoleMomentOfMagneticField,
        const double permeabilityOfSpaceBetweenBodies = basic_astrodynamics::physical_constants::MAGNETIC_CONSTANT );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_LOCAL_MAGNETIC_DIPOLE_FIELD_H
