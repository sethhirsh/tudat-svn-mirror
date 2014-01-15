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
 *      131110    S. Hirsh          File created.
 *
 *    References
 *    Davis, J., Mathematical Modeling of Earth’s Magnetic Field,
 *         http://hanspeterschaub.info/Papers/UnderGradStudents/MagneticField.pdf
 *    Fowler, C., The Solid Earth: An Introduction to Global Geophysics, 2nd Edition,
 *         Cambridge University Press, Cambridge, UK, 2005.
 *
 *    Notes
 *
 */

#include <cmath>
#include <Eigen/Geometry>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/ElectroMagnetism/localMagneticDipoleField.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute local magnetic field using static dipole model.
Eigen::Vector3d computeLocalStaticMagneticDipoleField(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const Eigen::Vector3d& originOfMagneticField,
        const Eigen::Vector3d& axisOfMagneticField,
        const double dipoleMomentOfMagneticField,
        const double permeabilityOfSpaceBetweenBodies )
{
    using std::acos;

    // Set constant values reused for optimal computation of acceleration components.
    const Eigen::Vector3d displacementOfBodyFromMagneticField = 
            positionOfBodySubjectToAcceleration 
            - originOfMagneticField;

    const double distanceBetweenBodies = 
            displacementOfBodyFromMagneticField.norm( );

    const double preMultiplier = permeabilityOfSpaceBetweenBodies
    		* dipoleMomentOfMagneticField 
    		/ std::pow( distanceBetweenBodies, 5.0 )
    		/ 4.0 / tudat::basic_mathematics::mathematical_constants::PI;

    // Compute angle between axis of magnetic field and z axis
    const double angleOfMagneticAxisDeclination = std::acos(
                Eigen::Vector3d( 0,0,1 )
                .dot( axisOfMagneticField ) );

    // Compute vector perpendicular to axis of magnetic field and z axis
    const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( 0,0,1 ).cross( axisOfMagneticField );

    // Compute rotated position of body subject to acceleration
    const Eigen::Vector3d rotatedPosition = 
            tudat::basic_mathematics::computeRotationOfPointAboutArbitraryAxis( 
                                        originOfMagneticField,
                                        -angleOfMagneticAxisDeclination,
                                        axisOfRotation,
                                        displacementOfBodyFromMagneticField );

    Eigen::Vector3d rotatedMagneticFieldDueToDipoleMoment = 
            Eigen::Vector3d::Constant( preMultiplier );

    //Compute x component of acceleration due to dipole moment of body exerting acceleration
    rotatedMagneticFieldDueToDipoleMoment( basic_astrodynamics::xCartesianPositionIndex )
    		*= ( rotatedPosition.x( )
    			- originOfMagneticField.x( ) )
    		* ( rotatedPosition.z( )
    			- originOfMagneticField.z( ) )
    		* 3.0;

    //Compute y component of acceleration due to dipole moment of body exerting acceleration
    rotatedMagneticFieldDueToDipoleMoment( basic_astrodynamics::yCartesianPositionIndex )
    		*= ( rotatedPosition.y( )
    			- originOfMagneticField.y( ) )
    		* ( rotatedPosition.z( )
    			- originOfMagneticField.z( ) )
    		* 3.0;

    //Compute z component of acceleration due to dipole moment of body exerting acceleration
    rotatedMagneticFieldDueToDipoleMoment( basic_astrodynamics::zCartesianPositionIndex )
    		*= ( rotatedPosition.z( )
    			- originOfMagneticField.z( ) )
    		* ( rotatedPosition.z( )
    			- originOfMagneticField.z( ) )
    		* 3.0
    		- distanceBetweenBodies * distanceBetweenBodies;


    // Undo rotation of magnetic dipole moment vector back to "unrotated position"
    const Eigen::Vector3d MagneticFieldDueToDipoleMoment = 
            tudat::basic_mathematics::computeRotationOfVectorAboutArbitraryAxis(
                originOfMagneticField,
                angleOfMagneticAxisDeclination,
                axisOfRotation,
                rotatedPosition,
                rotatedMagneticFieldDueToDipoleMoment );

    //Return the local magnetic field using static dipole model
    return MagneticFieldDueToDipoleMoment;


}

} // namespace electro_magnetism
} // namespace tudat
