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
 *      140516    S. Hirsh          File created.
 *
 *    References
 *
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_local_magnetic_dipole_field_model )

// Set radiation pressure at 1 AU [N/m^2].
const double radiationPressureAtOneAU = 4.56e-6;

// Set 1 AU in metres [m].
const double astronomicalUnitInMeters = 1.49598e11;

//! Test implementation of radiation force model.
BOOST_AUTO_TEST_CASE( testRadiationPressureForceModelGaneffData )
{
    // Benchmark data is obtained from MATLAB script (Ganeff, 2012).

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedRadiationPressureForce =
            Eigen::Vector3d( -0.975383093968723e-6, -0.975383093968723e-6, 0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.21;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.5;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( astronomicalUnitInMeters, astronomicalUnitInMeters, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Compute radiation pressure force [N].
    const Eigen::Vector3d computedRadiationPressureForce
                = electro_magnetism::computeCannonBallRadiationPressureForce(
                    radiationPressureAtTarget, positionVectorToSource,
                    areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Compare computed and expected radiation pressure force vectors.
    BOOST_CHECK_EQUAL( computedRadiationPressureForce.z( ),
                       expectedRadiationPressureForce.z( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce.segment( 0, 2 ),
                                       expectedRadiationPressureForce.segment( 0, 2 ),
                                       1.0e-15 );
}

//! Test acceleration model on sample satellite at approximately 1 AU away from the Sun.
BOOST_AUTO_TEST_CASE( testRadiationAccelerationModelEarth )
{

    // Benchmark data is obtained using the General Astrodynamics Library (Willmott, 2011).

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedRadiationPressureAcceleration =
            Eigen::Vector3d( -2.964e-06, 0.0, 0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.3;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 2.0;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 4.0;

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration
            = electro_magnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Compare computed and expected radiation pressure acceleration vectors
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration,
                                       expectedRadiationPressureAcceleration,
                                       1.0e-15 );
}

//! Test acceleration model on sample satellite at approximately the position of Venus with respect
//! to the Sun.
BOOST_AUTO_TEST_CASE( testRadiationAccelerationModelVenus )
{
    // Benchmark data is obtained using the General Astrodynamics Library (Willmott, 2011).
    // Tests satellite at approximately the distance of Venus away from the Sun.

    // Set the distance from the Sun to Venus.
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedRadiationPressureAcceleration =
            Eigen::Vector3d( -2.05147517201883e-05, 
                             -2.05147517201883e-05, 
                             0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.5;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.005;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( distanceSunToVenus / std::sqrt( 2.0 ), 
                               distanceSunToVenus / std::sqrt( 2.0 ),
                               0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 0.0022;

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration
            = electro_magnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Compare computed and expected radiation pressure acceleration vectors.
    BOOST_CHECK_SMALL( computedRadiationPressureAcceleration.z( ),
                       std::numeric_limits< double >::min( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration.segment( 0, 2 ),
                                       expectedRadiationPressureAcceleration.segment( 0, 2 ),
                                       1.0e-14 );
}

//! Test force model on sample planet at approximately the position of Uranus with respect to the 
//! Sun.
BOOST_AUTO_TEST_CASE( testRadiationForceModelUranus )
{
    // Benchmark data is obtained using the General Astrodynamics Library (Willmott, 2011).
    // Tests satellite at approximately the distance of Uranus away from the Sun. 

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedRadiationPressureForce =
            Eigen::Vector3d( -4470411.61112176, 
                             -4470411.61112176, 
                             0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.8;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 69939064094327.4;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( 9.529 * astronomicalUnitInMeters / std::sqrt( 2.0 ), 
                               9.529 * astronomicalUnitInMeters / std::sqrt( 2.0 ), 
                               0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Compute radiation pressure force [N].
    const Eigen::Vector3d computedRadiationPressureForce
            = electro_magnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Compare computed and expected radiation pressure force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce,
                                       expectedRadiationPressureForce,
                                       1.0e-14 );
}

//! Test force model on random satellite at a random distance with respect to the 
//! Sun.
BOOST_AUTO_TEST_CASE( testRadiationForceModelRandom )
{
    // Benchmark data is obtained using the General Astrodynamics Library (Willmott, 2011).
    // Tests satellite at random distance from the Sun.

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedRadiationPressureForce =
            Eigen::Vector3d( -3043733.21422537, 
                             -2929936.30441141, 
                             -473166.433773283 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.4058;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 514701.9505;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( 94359740.25, 
                               90831886.1, 
                               14668782.92 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Compute radiation pressure force [N].
    const Eigen::Vector3d computedRadiationPressureForce
            = electro_magnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Compare computed and expected radiation pressure force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureForce,
                                       expectedRadiationPressureForce,
                                       1.0e-14 );
}

//! Test radiation force model on a hand at approximately 1 AU from the Sun.
BOOST_AUTO_TEST_CASE( testRadiationForceModelGiancoliData )
{
    // Benchmark data is obtained from Physics for Scientists and Engineers
    // with Modern Physics Volume 2 (Ch. 31, Ex. 7) (Giancoli, 1985).

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedRadiationPressureForce =
            Eigen::Vector3d( 6.0e-8, 0.0, 0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.0;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 0.02;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Compute radiation pressure force [N].
    const Eigen::Vector3d computedRadiationPressureForce
            = electro_magnetism::computeCannonBallRadiationPressureForce(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient );

    // Compare computed and expected radiation pressure force vectors.            
    BOOST_CHECK_SMALL( std::fabs( computedRadiationPressureForce.x( )
                                  - expectedRadiationPressureForce.x( ) ),
                       1.0e-6 );
    BOOST_CHECK_SMALL( computedRadiationPressureForce.y( ),
                       std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( computedRadiationPressureForce.z( ),
                       std::numeric_limits< double >::min( ) );
}


//! Test radiation acceleration model on Ulysses satellite at 1AU.
BOOST_AUTO_TEST_CASE( testRadiationAccelerationModelUlysses )
{
    // Benchmark data obtained from (Irizarry, 2001).

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedRadiationPressureAcceleration 
        = Eigen::Vector3d( -1.713e-7, 0.0, 0.0 );

    // Set radiation pressure coefficient (1 + emissivity).
    const double radiationPressureCoefficient = 1.0 + 0.327;

    // Set area on target that is subject to radiation pressure [m^2].
    const double areaSubjectToRadiationPressure = 10.59;

    // Set position vector [m].
    Eigen::Vector3d positionVectorToSource
            = Eigen::Vector3d( astronomicalUnitInMeters, 0.0, 0.0 );

    // Set radiation pressure at target [N/m^2].
    const double radiationPressureAtTarget = radiationPressureAtOneAU
            * astronomicalUnitInMeters * astronomicalUnitInMeters
            / positionVectorToSource.squaredNorm( );

    // Normalize position vector to get vector pointing to source in non-dimensional units.
    positionVectorToSource.normalize( );

    // Set mass of accelerated body [kg].
    const double mass = 370.0;

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration
            = electro_magnetism::computeCannonBallRadiationPressureAcceleration(
                radiationPressureAtTarget, positionVectorToSource,
                areaSubjectToRadiationPressure, radiationPressureCoefficient, mass );

    // Compare computed and expected radiation pressure acceleration vectors. 
    BOOST_CHECK_SMALL( std::fabs( computedRadiationPressureAcceleration.x( )
                                  - expectedRadiationPressureAcceleration.x( ) ),
                       1.0e-8 );
    BOOST_CHECK_SMALL( computedRadiationPressureAcceleration.y( ),
                       std::numeric_limits< double >::min( ) );
    BOOST_CHECK_SMALL( computedRadiationPressureAcceleration.z( ),
                       std::numeric_limits< double >::min( ) );
}

//! Test class implementation of radiation pressure acceleration model.

// Set position of source of radiation pressure at origin [m].
static Eigen::Vector3d sourcePosition = Eigen::Vector3d::Zero( );

// Get position of source of radiation pressure [m].
Eigen::Vector3d getSourcePosition( ) { return sourcePosition; }

// Set position of accelerated body [m].
static Eigen::Vector3d acceleratedBodyPosition 
    = Eigen::Vector3d( -astronomicalUnitInMeters, 0.0, 0.0 );

// Get position of accelerated body [m].
Eigen::Vector3d getAcceleratedBodyPosition( ) { return acceleratedBodyPosition; }

// Get vector from accelerated body to source [m].
Eigen::Vector3d getVectorToSource( ) 
{ 
    return getSourcePosition( ) - getAcceleratedBodyPosition( ); 
}

// Set radiation pressure at location of acceleration body [N/m^2].
static double radiationPressure = radiationPressureAtOneAU
    * astronomicalUnitInMeters * astronomicalUnitInMeters / getVectorToSource( ).squaredNorm( );

// Get radiation pressure at location of acceleration body [N/m^2].
double getRadiationPressure( ) { return radiationPressure; }

// Set radiation pressure coefficient.
static double radiationPressureCoefficient = 1.0 + 0.3;

// Get radiation pressure coefficient.
double getRadiationPressureCoefficient( ) { return radiationPressureCoefficient; }

// Set area subject to radiation pressure [m^2].
static double areaSubjectToRadiationPressure = 2.0;

// Get area subject to radiation pressure [m^2].
double getAreaSubjectToRadiationPressure( ) { return areaSubjectToRadiationPressure; }

// Set mass of accelerated body [kg].
static double massOfAcceleratedBody = 4.0;

// Get mass of acceleration body [kg].
double getMassOfAcceleratedBody( ) { return massOfAcceleratedBody; }

//! Test radiation pressure acceleration model constructor. 
BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationModelClassConstructor )
{
    // Declare and initialize cannon-ball radiation pressure acceleration model.
    electro_magnetism::CannonBallRadiationPressurePointer radiationPressureModel
        = boost::make_shared< electro_magnetism::CannonBallRadiationPressure >( 
            &getSourcePosition, &getAcceleratedBodyPosition, &getRadiationPressure,
            &getRadiationPressureCoefficient, &getAreaSubjectToRadiationPressure, 
            &getMassOfAcceleratedBody );

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedRadiationPressureAcceleration 
        = Eigen::Vector3d( -2.964e-06, 0.0, 0.0 );

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration 
        = radiationPressureModel->getAcceleration( );

    // Compare computed and expected radiation pressure acceleration vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration,
                                       expectedRadiationPressureAcceleration,
                                       1.0e-15 );      
}

//! Test radiation pressure acceleration model update-members function. 
BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationModelClassUpdateMembers )
{
    // Declare and initialize cannon-ball radiation pressure acceleration model.
    electro_magnetism::CannonBallRadiationPressurePointer radiationPressureModel
        = boost::make_shared< electro_magnetism::CannonBallRadiationPressure >( 
            &getSourcePosition, &getAcceleratedBodyPosition, &getRadiationPressure,
            &getRadiationPressureCoefficient, &getAreaSubjectToRadiationPressure, 
            &getMassOfAcceleratedBody );

    // Set expected radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d expectedRadiationPressureAcceleration 
        = Eigen::Vector3d( -2.05147517201883e-05, 
                           -2.05147517201883e-05, 
                           0.0 );

    // Set the distance from the Sun to Venus [m].
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;      

    // Update position of accelerated body [m].
    acceleratedBodyPosition = Eigen::Vector3d( -distanceSunToVenus / std::sqrt( 2.0 ), 
                                               -distanceSunToVenus / std::sqrt( 2.0 ),
                                               0.0 );
    
    // Update radiation pressure at location of accelerated body [N/m^2].
    radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters 
        / getVectorToSource( ).squaredNorm( );

    // Update radiation pressure coefficient.
    radiationPressureCoefficient = 1.0 + 0.5;  

    // Update area subject to radiation pressure [m^2].
    areaSubjectToRadiationPressure = 0.005;

    // Update mass of accelerated body [kg].
    massOfAcceleratedBody = 0.0022;

    // Update class members.
    radiationPressureModel->updateMembers( );

    // Compute radiation pressure acceleration [m/s^2].
    const Eigen::Vector3d computedRadiationPressureAcceleration 
        = radiationPressureModel->getAcceleration( );

    // Compare computed and expected radiation pressure acceleration vectors.
    BOOST_CHECK_SMALL( computedRadiationPressureAcceleration.z( ),
                       std::numeric_limits< double >::min( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRadiationPressureAcceleration.segment( 0, 2 ),
                                       expectedRadiationPressureAcceleration.segment( 0, 2 ),
                                       1.0e-14 );      
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
