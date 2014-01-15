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
#include "Tudat/Astrodynamics/ElectroMagnetism/localMagneticDipoleField.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/corotationalElectricAcceleration.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/lorentzStaticMagneticAcceleration.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/rotatingMagneticDipoleAcceleration.h"
#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h" 

namespace tudat
{
namespace electro_magnetism
{

//! Compute acceleration due to rotating magnetic dipole field
Eigen::Vector3d computeRotatingMagneticDipoleFieldAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& originOfMagneticField,
        const Eigen::Vector3d& originOfRotation,
        const Eigen::Vector3d& axisOfMagneticField,
        const Eigen::Vector3d& angularVelocityOfRotatingBody,
        const double chargeOfBodySubjectToAcceleration,
        const double permeabilityOfSpaceBetweenBodies,
        const double dipoleMomentOfMagneticField,
        const double massOfBodySubjectToAcceleration,
        const double timeElapsed )
{

// Compute angle of rotation.
const double angleOfRotation = angularVelocityOfRotatingBody.norm( ) 
                               * timeElapsed;

// Compute axis of rotation.
const Eigen::Vector3d axisOfRotation = angularVelocityOfRotatingBody.normalized( );


// Rotate position of body to find equivalent position in static dipole field
const Eigen::Vector3d equivalentPositionInStaticField = 
                tudat::basic_mathematics::computeRotationOfPointAboutArbitraryAxis(
                originOfMagneticField,
                -angleOfRotation,
                axisOfRotation,
                positionOfBodySubjectToAcceleration );

// Compute local magnetic field in equivalent static dipole field
const Eigen::Vector3d equivalentLocalStaticeMagneticField = 
                computeLocalStaticMagneticDipoleField(
                equivalentPositionInStaticField,
                originOfMagneticField,
                axisOfMagneticField,
                dipoleMomentOfMagneticField,
                permeabilityOfSpaceBetweenBodies );

// Compute rotation of equivalent local magnetic field to actual position
const Eigen::Vector3d rotatedLocalStaticMagneticField = 
                tudat::basic_mathematics::computeRotationOfVectorAboutArbitraryAxis(
                originOfRotation,
                angleOfRotation,
                axisOfRotation,
                equivalentPositionInStaticField,
                equivalentLocalStaticeMagneticField );

// Compute acceleration due to local static magnetic field
const Eigen::Vector3d accelerationDueToStaticMagneticField = 
                computeLorentzAccelerationDueToStaticMagneticField(
                velocityOfBodySubjectToAcceleration,
                rotatedLocalStaticMagneticField,
                chargeOfBodySubjectToAcceleration,
                massOfBodySubjectToAcceleration );

// Compute acceleration due to corotating electric field
const Eigen::Vector3d accelerationDueToCorotatingElectricField = 
                computeCorotationalElectricAcceleration(
                angularVelocityOfRotatingBody,
                positionOfBodySubjectToAcceleration,
                originOfRotation,
                rotatedLocalStaticMagneticField,
                chargeOfBodySubjectToAcceleration,
                massOfBodySubjectToAcceleration );

// Return acceleration due to static magnetic field and corotation electric field
return accelerationDueToStaticMagneticField + accelerationDueToCorotatingElectricField;

}

//! Get rotating magnetic dipole field acceleration.
Eigen::Vector3d RotatingMagneticDipoleField::getAcceleration( )
{
    return computeRotatingMagneticDipoleFieldAcceleration(
                currentAcceleratedBodyPosition_,
                currentAcceleratedBodyVelocity_,
                currentOriginOfMagneticField_,
                currentOriginOfRotation_,
                currentAxisOfMagneticField_,
                currentAngularVelocityOfRotatingBody_,
                currentChargeOfAcceleratedBody_,
                currentPermeabilityOfSpaceBetweenBodies_,
                currentDipoleMomentOfMagneticField_,
                currentMass_,
                currentTimeElapsed_ );
}

//! Update member variables used by the rotating magnetic dipole field acceleration model.
void RotatingMagneticDipoleField::updateMembers( )
{
    currentAcceleratedBodyPosition_ = acceleratedBodyPositionFunction_( );
    currentAcceleratedBodyVelocity_ = acceleratedBodyVelocityFunction_( );
    currentOriginOfMagneticField_ = originOfMagneticFieldFunction_( );
    currentAxisOfMagneticField_ = axisOfMagneticFieldFunction_( );
    currentAngularVelocityOfRotatingBody_ = angularVelocityOfRotatingBodyFunction_( );
    currentChargeOfAcceleratedBody_ = chargeOfAcceleratedBodyFunction_( );
    currentPermeabilityOfSpaceBetweenBodies_ = permeabilityOfSpaceBetweenBodiesFunction_( );
    currentDipoleMomentOfMagneticField_ = dipoleMomentOfMagneticFieldFunction_( );
    currentMass_ = massFunction_( );
    currentTimeElapsed_ = timeElapsedFunction_( );
}




} // namespace electro_magnetism
} // namespace tudat
