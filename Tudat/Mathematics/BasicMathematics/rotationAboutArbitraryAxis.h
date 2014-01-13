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

#ifndef TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H
#define TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

//! Compute rotation of point about arbitrary axis. 
/*!
 * Computes the rotation of a point about an arbitrary axis. Units for originOfRotation,
 * and initialPositionOfPoint must be the same. 
 * \param originOfRotation position origin of rotation axis.								
 * \param angleOfRotation angle which the point rotates with respect to the axis of rotation   [rad]
 * \param axisOfRotation unit vector point in the direction of rotation obtained using 
 *			right-hand rule																
 * \param initialPositionOfPoint position of point before rotation.					
 * \return position of point after rotation about axisOfRotation.
 */

//! Compute rotation of point about arbitrary axis
Eigen::Vector3d computeRotationOfPointAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfPoint );


//! Compute rotation of vector about an arbitrary axis.
/*!
 * Computes rotation of vector about an arbitrary axis. Uses 
 * computeRotationOfPointAboutArbitraryAxis to rotate the position of the head and tail of the 
 * vector. Units for originOfRotation, initialPositionOfVectorTail and initialVector must be the
 * same. 						
 * \param originOfRotation position of axis over which the body is being rotated.										
 * \param angleOfRotation angle over which the vector is rotated. A positive angle is determined
 * 			using the right hand rule.														  [rad]                                                      
 * \param axisOfRotation unit vector over which the rotated vector will be rotated.
 * \param initialPositionOfVectorTail initial position of tail of vector which will be rotated. 
 * \param initialVector vector which will be rotated. 
 * \return vector after rotation about an arbitrary axis.
 * \sa computeRotationOfPointAboutArbitraryAxis.
 */
//! Compute rotation of vector about arbitrary axis
Eigen::Vector3d computeRotationOfVectorAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfVectorTail,
        const Eigen::Vector3d& initialVector );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H
