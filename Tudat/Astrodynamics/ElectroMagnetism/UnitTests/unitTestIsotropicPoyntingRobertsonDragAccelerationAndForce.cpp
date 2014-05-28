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
 *      131031    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/ElectroMagnetism/isotropicPoyntingRobertsonDragForce.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/isotropicPoyntingRobertsonDragAcceleration.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_poynting_robertson_drag_acceleration_and_force_models )

//! Test implementation of Poynting-Robertson drag force model.
BOOST_AUTO_TEST_CASE( testPoyntingRobertsonDragForce )
{
	// Set expected magnitude of Poynting-Robertson drag force
	const Eigen::Vector3d expectedPoyntingRobertsonDragForce = 
				Eigen::Vector3d( 0.0, -1.0e-11, 0.0 );

	// Set radiation pressure force
	const double radiationPressureForce = 1.0e-6;

	// Gravitational parameter of the Sun
	const double solarGravitationalParameter = 1.327124400180e20;

	// Set astronomical unit in meters
	const double astronomicalUnitInMeters = 149597870.0e3;

	// Set velocity of accelerated body
	const Eigen::Vector3d velocityOfAcceleratedBody = Eigen::Vector3d(
                     0.0, sqrt(solarGravitationalParameter / astronomicalUnitInMeters), 0.0 );

	const Eigen::Vector3d vectorFromSource = Eigen::Vector3d( 1.0, 0.0, 0.0 );

	
	const Eigen::Vector3d computedPoyntingRobertsonDragForce = 
		electro_magnetism::computeIsotropicPoyntingRobertsonDragForce(
				velocityOfAcceleratedBody, vectorFromSource, radiationPressureForce);

	// Compare computed and expected radiation pressure force vectors.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedPoyntingRobertsonDragForce,
                                       expectedPoyntingRobertsonDragForce,
                                       10.0 );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat