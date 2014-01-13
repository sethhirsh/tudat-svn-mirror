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
 *      140110    S. Hirsh          File created.
 *
 *    References
 *    Hameduddin, I. Rotate vector(s) about axis, rodrigues_rot.m, available at
 *        http://www.mathworks.com/matlabcentral/fileexchange/34426-rotate-vectors-about-axis,
 *        2012, last accessed: 11th January, 2014.
 *    Murray, G. Rotation About an Arbitrary Axis in 3 Dimensions, 
 *        http://twist-and-shout.appspot.com/, 2013. last accessed: 11th January, 2014. 
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h"

#include "TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for rotations about about arbitrary axes.
BOOST_AUTO_TEST_SUITE( test_RotationAboutArbitraryAxis )

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_PointRotationWithCommonOrigin )
{

  //Benchmark data is obtained using Matlab Script (Hameduddin, 2012)

  //Set expected rotated vector
  const Eigen::Vector3d expectedRotationPosition = Eigen::Vector3d( 0.0, 1.414213562373095, 1.0 );

  //Set origin of rotation
  const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 0.0, 0.0, 0.0 );

  //Set angle of rotation [rad]
  const double angleOfRotation = tudat::basic_mathematics::mathematical_constants::PI / 4.0;

  //Set axis of rotation
  const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( 0.0, 0.0, 1.0 );

  //Set initial position of point
  const Eigen::Vector3d initialPositionOfPoint = Eigen::Vector3d( 1.0, 1.0, 1.0 );

  //Compute rotated position.
  const Eigen::Vector3d computedRotatedPosition 
              = basic_mathematics::computeRotationOfPointAboutArbitraryAxis(
                  originOfRotation,  angleOfRotation, axisOfRotation, initialPositionOfPoint );

  // Compare computed and expected radiation pressure acceleration vectors
  BOOST_CHECK_SMALL( computedRotatedPosition.x( ),
                     std::numeric_limits<double>::epsilon( ) );

  TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedPosition.segment( 1, 2 ),
                                     expectedRotationPosition.segment( 1, 2),
                                     std::numeric_limits<double>::epsilon( ) );

}

BOOST_AUTO_TEST_CASE( test_RotationAboutArbitraryAxis_PointRotationWithDifferentOrigins )
{

  //Benchmark data is obtained using http://twist-and-shout.appspot.com/ (Murray, 2013).

  //Set expected rotated vector
  const Eigen::Vector3d expectedRotationPosition = Eigen::Vector3d( 3.1566, 
                                                                    -5.9715, 
                                                                    -4.4187 );
  //Set origin of rotation
  const Eigen::Vector3d originOfRotation = Eigen::Vector3d( 4.0, 1.0, -1.0 );

  //Set angle of rotation [rad]
  const double angleOfRotation = 7.0;

  //Set axis of rotation
  const Eigen::Vector3d axisOfRotation = Eigen::Vector3d( 2.0, -2.0, 3.0 );

  //Set initial position of point
  const Eigen::Vector3d initialPositionOfPoint = Eigen::Vector3d( -1.0, -5.0, -1.0 );

  //Compute rotated position.
  const Eigen::Vector3d computedRotatedPosition 
              = basic_mathematics::computeRotationOfPointAboutArbitraryAxis(
                  originOfRotation,  angleOfRotation, axisOfRotation, initialPositionOfPoint );

  // Compare computed and expected radiation pressure acceleration vectors
  TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedRotatedPosition,
                                     expectedRotationPosition,
                                     1.0e-4 );

}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
