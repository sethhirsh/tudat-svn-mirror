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
   *      YYMMDD    Author              Comment
   *      130219    D. Dirkx            File created.
   *      130227    R.C.A. Boon         Removed using directives (using namespace), added necessary
   *                                    using declarations (using ...::member), removed obsolete
   *                                    include statements
   *      130302    D. Dirkx            Updated file by putting tests in local scope;
   *                                    expanded tests.
   *
   *    References
   *
   *    Notes
   *      The unit tests use rotation matrices generated by spice, the code for the spice interface
   *      used to generate these testing values are included in this file, but commented.
   */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace unit_tests
{

using tudat::unit_conversions::convertDegreesToRadians;
using tudat::ephemerides::SimpleRotationalEphemeris;
using tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000;
using tudat::basic_astrodynamics::JULIAN_DAY_AT_0_MJD;
using tudat::reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion;

BOOST_AUTO_TEST_SUITE( test_simple_rotational_ephemeris )

// Test simple rotational ephemeris class by compariing results to Venus-fixed frame with Spice.
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemeris )
{    
    // Data from pck00010.tpc spice kernel.
    const double venusPoleRightAscension = convertDegreesToRadians( 272.76 );
    const double venusPoleDeclination = convertDegreesToRadians( 67.16 );
    const double venusPrimeMeridianAtJ2000 = convertDegreesToRadians( 160.20 );
    const double venusRotationRate = convertDegreesToRadians( -1.4813688 ) /
            basic_astrodynamics::physical_constants::JULIAN_DAY;
    
    // Define names of frames.
    const std::string baseFrame = "J2000";
    const std::string targetFrame = "IAU_VENUS";

    // Calculate initial rotation quaternion to frame.
    const Eigen::Quaterniond initialRotationToTargetFrame =
            getInertialToPlanetocentricFrameTransformationQuaternion(
                venusPoleDeclination, venusPoleRightAscension, venusPrimeMeridianAtJ2000 );

    // Create rotational ephemeris objects from (one for each type of constructor)
    SimpleRotationalEphemeris venusRotationalEphemerisFromAngles(
                venusPoleRightAscension, venusPoleDeclination, venusPrimeMeridianAtJ2000,
                venusRotationRate, 0.0, JULIAN_DAY_ON_J2000, baseFrame, targetFrame );
    SimpleRotationalEphemeris venusRotationalEphemerisFromInitialState(
                initialRotationToTargetFrame,
                venusRotationRate, 0.0, JULIAN_DAY_ON_J2000, baseFrame, targetFrame );


    // Test initial rotation matrix, both from the direct constructed object, and the one from
    // constituent angles.
    {
        // The following code block can be used to retrieve the benchmark data from Spice.
        /*
        loadSpiceKernelInTudat( getTudatRootPath( ) +
              "Astrodynamics/Ephemerides/UnitTests/pck00010.tpc" );
        Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
        Eigen::Quaterniond spiceInitialRotationToTargetFrame =
                computeRotationQuaternionBetweenFrames( baseFrame, targetFrame, 0.0 );
       */

        // Set initial rotation, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
        spiceInitialRotationToTargetFrameMatrix
                << -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
                -0.296591573568662, -0.882413772579987, -0.3652114078848295,
                0.01869081416890202, -0.3877088083617989, 0.9215923900425705;
        
        // Check calculated with spice initial state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceInitialRotationToTargetFrameMatrix,
                                           Eigen::Matrix3d( initialRotationToTargetFrame ),
                                           1.0E-15 );
        
        // Check ephemeris calculations with Spice initial state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    Eigen::Matrix3d(
                        venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                            0.0, JULIAN_DAY_ON_J2000 ) ),
                    Eigen::Matrix3d( initialRotationToTargetFrame ),
                    1.0E-15 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    Eigen::Matrix3d(
                        venusRotationalEphemerisFromInitialState.getRotationToTargetFrame(
                            0.0, JULIAN_DAY_ON_J2000 ) ),
                    Eigen::Matrix3d( initialRotationToTargetFrame ),
                    1.0E-15 );
    }

    // Set time at which rotational ephemeris it to be called for subsequent tests.
    const double secondsSinceJ2000 = 1.0E6;

    // Test rotation to target frame at specified time.
    {
        // The following code block can be used to retrieve the benchmark data from Spice.
        /*
        Eigen::Quaterniond spiceInitialRotationToTargetFrame =
               computeRotationQuaternionBetweenFrames( baseFrame, targetFrame, secondsSinceJ2000 );
        */

        // Set rotation at given time, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceRotationMatrix;
        spiceRotationMatrix << -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
                -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
                0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        // Calculate rotations to frame at certain time;
        // Check Spice result with ephemerides results.
        Eigen::Quaterniond ephemerisRotation =
                venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                    secondsSinceJ2000, JULIAN_DAY_ON_J2000 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Matrix3d( ephemerisRotation  ),
                                           spiceRotationMatrix, 5.0E-15 );

        ephemerisRotation =
                venusRotationalEphemerisFromInitialState.getRotationToTargetFrame(
                    secondsSinceJ2000, JULIAN_DAY_ON_J2000 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Matrix3d( ephemerisRotation ),
                                           spiceRotationMatrix, 5.0E-15 );
    }

    // Test rotation from target frame at specified time (is checked by checking if it is inverse
    // of opposite rotation).
    {
        // Test orthonormality of matrices from object created from initial angles.
        Eigen::Matrix3d productOfOppositeRotations =
                Eigen::Matrix3d( venusRotationalEphemerisFromAngles.getRotationToTargetFrame(
                                     secondsSinceJ2000, JULIAN_DAY_ON_J2000 ) ) *
                Eigen::Matrix3d( venusRotationalEphemerisFromAngles.getRotationToBaseFrame(
                                     secondsSinceJ2000, JULIAN_DAY_ON_J2000 ) );

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( productOfOppositeRotations( i, j ) -
                                   Eigen::Matrix3d::Identity( )( i, j ), 1.0E-15 );
            }
        }

        // Test orthonormality of matrices from object created from initial state.
        productOfOppositeRotations =
                Eigen::Matrix3d( venusRotationalEphemerisFromInitialState.getRotationToTargetFrame(
                                     secondsSinceJ2000, JULIAN_DAY_ON_J2000 ) ) *
                Eigen::Matrix3d( venusRotationalEphemerisFromInitialState.getRotationToBaseFrame(
                                     secondsSinceJ2000, JULIAN_DAY_ON_J2000 ) );

        for( int i = 0; i < 3; i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( productOfOppositeRotations( i, j ) -
                                   Eigen::Matrix3d::Identity( )( i, j ), 1.0E-15 );
            }
        }
    }

    // Test rotational ephemeris in case input epoch is not reference epoch.
    {
        // The following code block can be used to retrieve the benchmark data from Spice.
        /*
        Eigen::Quaterniond spiceInitialRotationToTargetFrame =
               computeRotationQuaternionBetweenFrames( baseFrame, targetFrame, secondsSinceJ2000 );
        */

        // Set rotation at given time, as calculated with Spice (see above commented lines)
        Eigen::Matrix3d spiceRotationMatrix;
        spiceRotationMatrix
                << -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
                -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
                0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        // Test alternative input reference epoch on object from initial angles.
        double secondsSinceMjd0 = ( JULIAN_DAY_ON_J2000 - JULIAN_DAY_AT_0_MJD ) *
                basic_astrodynamics::physical_constants::JULIAN_DAY + secondsSinceJ2000;

        Eigen::Quaterniond ephemerisRotation = venusRotationalEphemerisFromAngles.
                getRotationToTargetFrame( secondsSinceMjd0, JULIAN_DAY_AT_0_MJD );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Matrix3d( ephemerisRotation ),
                                           spiceRotationMatrix, 5.0E-15 );

        // Test alternative input reference epoch on object from initial state.
        ephemerisRotation = venusRotationalEphemerisFromInitialState.getRotationToTargetFrame(
                    secondsSinceMjd0, JULIAN_DAY_AT_0_MJD );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( Eigen::Matrix3d( ephemerisRotation ),
                                           spiceRotationMatrix, 5.0E-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
