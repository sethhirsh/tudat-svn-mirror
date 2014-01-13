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
 *      131204    S. Hirsh          File created.
 *
 *    References
 *    Wikipedia. Quaternions and spatial rotation, http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation,
 *         last accessed: 4 December, 2013.
 *
 *    Notes
 *
 */

#include <cmath>
#include <Eigen/Geometry>
#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h"

namespace tudat
{
namespace basic_mathematics
{

//! Compute rotation of point about arbitrary axis
Eigen::Vector3d computeRotationOfPointAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfPoint )
{

    // Compute rotation matrix using AngleAxis object.
    Eigen::Matrix3d rotationMatrix;

    rotationMatrix = Eigen::AngleAxisd(
                angleOfRotation,
                axisOfRotation.normalized( ) );


    // Compute displacement of position of point from the origin of rotation.
    const Eigen::Vector3d positionDisplacementFromOriginOfRotation = 
                initialPositionOfPoint 
                - originOfRotation;

    // Compute rotation of point about axis of rotation.
    const Eigen::Vector3d rotatedPositionWithRespectToOrigin = 
                rotationMatrix 
                * positionDisplacementFromOriginOfRotation;

    //Return position after rotation about arbitrary axis.
    return rotatedPositionWithRespectToOrigin + originOfRotation;


}

//! Compute rotation of vector about arbitrary axis
Eigen::Vector3d computeRotationOfVectorAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfVectorTail,
        const Eigen::Vector3d& initialVector )
{

// Compute rotation of the tail of vector.
Eigen::Vector3d rotatedPositionOfVectorTail = computeRotationOfPointAboutArbitraryAxis(
            originOfRotation,
            angleOfRotation,
            axisOfRotation,
            initialPositionOfVectorTail );

// Compute rotation of the head of vector.
Eigen::Vector3d rotatedPositionOfVectorHead = computeRotationOfPointAboutArbitraryAxis(
            originOfRotation,
            angleOfRotation,
            axisOfRotation,
            initialPositionOfVectorTail + initialVector);

// Return rotated vector
return rotatedPositionOfVectorHead - rotatedPositionOfVectorTail;

}

} // namespace basic_mathematics
} // namespace tudat
