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
 *    Dahle, M. magField_of_the_earth.py, available at 
 *       http://folk.ntnu.no/magnud/ElMag_project/python_files/, 2012.
 *    Öztürk, K. Python program to follow particle trajectories in tail-like magnetic fields,  
 *       speiserorbits.py, available at https://sites.google.com/site/mkaanozturk/programs, 2011.
 *      
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

#include "Tudat/Astrodynamics/ElectroMagnetism/localMagneticDipoleField.h"
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_local_magnetic_dipole_field_model )

//! Test implementation of radiation force model.
BOOST_AUTO_TEST_CASE( testAlignedDipoleModelEarth )
{
    // Benchmark data is obtained from MATLAB script (Ganeff, 2012).

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedLocalMagneticDipoleField =
            Eigen::Vector3d( -0.0004237415112954681, 
                             -0.0004237415112954681, 
                             0.000538084458787896 );

    // Set radius of Earth
    const double radiusOfEarth = 6378137.0;

    // Set position of body subject to acceleration
    const Eigen::Vector3d positionOfBodySubjectToAcceleration = 
            Eigen::Vector3d(7000.0, 7000.0, 3000.0);

    // Set origin of magnetic field
    const Eigen::Vector3d originOfMagneticField = Eigen::Vector3d(0.0,0.0,0.0);

    // Set axis of magnetic field
    const Eigen::Vector3d axisOfMagneticField = Eigen::Vector3d(0.0,0.0,1.0);

    // Approximate mean value of the magnetic field at the magnetic equator on the Earth's surface
    const double meanMagneticFieldAtEquator = 3.07e-5;

    // Set dipole moment of magnetic field (what's the deal with the minus sign?)
    const double dipoleMomentOfMagneticField = meanMagneticFieldAtEquator * -1.0 * 
                                               radiusOfEarth*radiusOfEarth*radiusOfEarth;


    const Eigen::Vector3d computedLocalDipoleField
                = electro_magnetism::computeLocalStaticMagneticDipoleField(
                    positionOfBodySubjectToAcceleration, originOfMagneticField, 
                    axisOfMagneticField, dipoleMomentOfMagneticField );


    // Compare computed and expected local magnetic dipole field vectors.
    TUDAT_CHECK_MATRIX_CLOSE( computedLocalDipoleField,
                              expectedLocalMagneticDipoleField, 1.0e-9 );

}

//! Test implementation of radiation force model.
BOOST_AUTO_TEST_CASE( testTiltedMagneticDipole )
{
    // Benchmark data is obtained from MATLAB script (Ganeff, 2012).

    // Set expected radiation pressure force [N].
    const Eigen::Vector3d expectedLocalMagneticDipoleField =
            Eigen::Vector3d( -1.15013283e-5, 
                             0.0, 
                             -2.95846024e-5 );

    // Set radius of Earth
    const double radiusOfEarth = 6.37e6;

    // Set position of body subject to acceleration
    const Eigen::Vector3d positionOfBodySubjectToAcceleration = 
            Eigen::Vector3d(radiusOfEarth, 0.0, 0.0);

    // Set origin of magnetic field
    const Eigen::Vector3d originOfMagneticField = Eigen::Vector3d(0.0,0.0,0.0);

    // Angle between the z-axis and Earth's magnetic dipole moment
    const double declinationOfMagneticDipoleAxis = -11.0 / 180.0 * 
                        tudat::basic_mathematics::mathematical_constants::PI; 

    // Set axis of magnetic field
    const Eigen::Vector3d axisOfMagneticField = Eigen::Vector3d( sin(declinationOfMagneticDipoleAxis) ,0.0, cos(declinationOfMagneticDipoleAxis));

    // Set dipole moment of magnetic field (what's the deal with the minus sign?) [A m^2]
    const double dipoleMomentOfMagneticField = 7.79e22;

    const Eigen::Vector3d computedLocalMagneticDipoleField
                = electro_magnetism::computeLocalStaticMagneticDipoleField(
                    positionOfBodySubjectToAcceleration, originOfMagneticField, 
                    axisOfMagneticField, dipoleMomentOfMagneticField );


    // Compare computed and expected local magnetic dipole field vectors.
    BOOST_CHECK_SMALL(computedLocalMagneticDipoleField.x( ) 
                    - expectedLocalMagneticDipoleField.x( ), 1.0e-9);

    BOOST_CHECK_SMALL( computedLocalMagneticDipoleField.y( ),
                       std::numeric_limits< double >::min( ) );

    BOOST_CHECK_SMALL(computedLocalMagneticDipoleField.z( ) 
                    - expectedLocalMagneticDipoleField.z( ), 1.0e-9);

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
