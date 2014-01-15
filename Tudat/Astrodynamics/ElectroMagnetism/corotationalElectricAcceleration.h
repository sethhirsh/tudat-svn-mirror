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

#ifndef TUDAT_COROTATIONAL_ELECTRIC_ACCELERATION_H
#define TUDAT_COROTATIONAL_ELECTRIC_ACCELERATION_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute Lorentz force due to source with corotating electric field.
/*!
 * Computes Lorentz force due to corotation electric field.
 * Assumes plasma about rotating body has angular velocity equal to the rotating body.
 * \f[
 *		\bar{f} = q * \bar{\Omega} x \bar{r} x \bar{B}
 * \f] 
 * where \f$f\f$ is the force on the particle, \f$q\f$ is the charge of the accelerating
 * particle, \f$\bar{\Omega}\f$ is the angular velocity of the rotating body, \f$\bar{r}\f$
 * is the position of the body subject to the corotation electric field, and \f$\bar{v}\f$
 * is the magnetic field in the local region of the particle. 
 * \param angularVelocityOfRotatingBody angular velocity of body which is producing the 
 * 			corotation electric field                                    						[rad/s]
 * \param positionOfBodySubjectToAcceleration position of body in an inertial reference frame
 * 			subject to acceleration by corotation electric field.								[m]
 * \param originOfMagneticField origin of magnetic field in an inertial reference frame 		[m]
 * \param localMagneticField local magnetic field at position of body subject to acceleration	[T·m]
 * \param chargeOfBodySubjectToAcceleration charge of body which is being accelerated by
 *			Lorentz force 																		[C]
 * \return Lorentz force due to corotating electric field			                            [N]
 */

Eigen::Vector3d computeCorotationalElectricAcceleration(
        const Eigen::Vector3d& angularVelocityOfRotatingBody,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const Eigen::Vector3d& originOfRotation,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration,
        const double massOfBodySubjectToAcceleration );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_COROTATIONAL_ELECTRIC_ACCELERATION_H
