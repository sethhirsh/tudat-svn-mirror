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
 *      131109	  S. Hirsh	    File Created. 
 *
 *    References
 *		
 *
 *    Notes
 *
 */

#ifndef TUDAT_ROTATING_MAGNETIC_DIPOLE_ACCELERATION_H
#define TUDAT_ROTATING_MAGNETIC_DIPOLE_ACCELERATION_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute acceleration of particle due to rotating magnetic dipole field.
/*!
 * Computes acceleration of particle due to static magnetic dipole field and the 
 * corresponding corotational electric field. Assumes particle to be point-like and 
 * the magnetic field rotates at a constant angular velocity.
 * \param positionOfBodySubjectToAcceleration position of body accelerated 
 *              by the magnetic field                                                            [m]
 * \param velocityOfBodySubjectToAcceleration velocity of body accelerated 
 *              by the magnetic field                                                            [m/s]
 * \param originOfMagneticField origin of magnetic field of body that produces
 *              magnetic dipole moment                                                           [m]
 * \param originOfRotation origin of axis about which the field the is rotating                  [m]
 * \param axisOfMagneticField unitvector that points in the direction of the symmetry axis of
 *              magnetic dipole field.                                                           [-]
 * \param angularVelocityOfRotatingBody angular velocity of body which is producing the 
 *               corotation electric field. Note: This value is assumed to remain constant.      [rad/s]
 * \param permeabilityOfSpaceBetweenBodies permeability of space between bodies.                 [V·s/(A·m)]                                                                          [rad/s]                                                                                            [-]
 * \param dipoleMomentOfMagneticField magnetic dipole moment from body exerting acceleration     [T·m]
 * \param massOfBodySubjectToAcceleration mass of body accelerated by the magnetic field         [kg]  	
 * \param timeElapsed time elapsed since static magnetic field was in initial orientation        [s]				
 * \return acceleration of particle due to rotating magnetic field                               [N]
 */
//! Compute acceleration of particle due to rotating magnetic dipole field
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
        const double timeElapsed );

//! Rotating magnetic dipole field acceleration model class.
/*!
 * Class that can be used to compute the rotating magnetic dipole field acceleration. 
 * The particle subject to the magnetic field is assumed to be point-like and the angular
 * velocity of the rotating body is assumed to be constant.
 */
class RotatingMagneticDipoleField : public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef boost::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param acceleratedBodyPositionFunction returning position of accelerated body.
     * \param acceleratedBodyVelocityFunction Function returning velocity of accelerated body.
     * \param originOfMagneticFieldFunction Function returning origin of magnetic field of body producing the
     *              magnetic dipole moment
     * \param originOfRotationFunction Function returning origin of axis about which the field the is rotating
     * \param axisOfMagneticFieldFunction Function returning unitvector that points in the direction of the 
     *              symmetry axis of magnetic dipole field.
     * \param angularVelocityOfRotatingBodyFunction Function returning angular velocity of body which producing the 
     *              corotation electric field.
     * \param chargeOfAcceleratedBodyFunction Function returning charge of body accelerated by magnetic field.
     * \param permeabilityOfSpaceBetweenBodiesFunction Function returning permeability of space between bodies.
     * \param dipoleMomentOfMagneticFieldFunction Function returning magnetic dipole moment from body exerting acceleration.
     * \param massFunction Function returning current mass of body undergoing acceleration.
     * \param timeElapsed Function returning current time elapsed since static magnetic field was in initial orientation 
     */
    RotatingMagneticDipoleField(
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction,
            Vector3dReturningFunction originOfMagneticFieldFunction,
            Vector3dReturningFunction originOfRotationFunction,
            Vector3dReturningFunction axisOfMagneticFieldFunction,
            Vector3dReturningFunction angularVelocityOfRotatingBodyFunction,            
            DoubleReturningFunction chargeOfAcceleratedBodyFunction,
            DoubleReturningFunction permeabilityOfSpaceBetweenBodiesFunction,
            DoubleReturningFunction dipoleMomentOfMagneticFieldFunction,
            DoubleReturningFunction massFunction,
            DoubleReturningFunction timeElapsedFunction )
        : acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          originOfMagneticFieldFunction_( originOfMagneticFieldFunction ),
          originOfRotationFunction_( originOfRotationFunction ),
          axisOfMagneticFieldFunction_( axisOfMagneticFieldFunction ),
          angularVelocityOfRotatingBodyFunction_( angularVelocityOfRotatingBodyFunction ),
          chargeOfAcceleratedBodyFunction_( chargeOfAcceleratedBodyFunction ),
          permeabilityOfSpaceBetweenBodiesFunction_( permeabilityOfSpaceBetweenBodiesFunction ),
          dipoleMomentOfMagneticFieldFunction_( dipoleMomentOfMagneticFieldFunction ),         
          massFunction_( massFunction ),
          timeElapsedFunction_( timeElapsedFunction )
    {
        this->updateMembers( );
    }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for position and velocity vectors and elapsed time
     * and constant values for other parameters. 
     * \param acceleratedBodyPosition Function returning position of accelerated body.
     * \param acceleratedBodyVelocityFunction Function returning velocity of accelerated body.
     * \param originOfMagneticField Constant origin of magnetic field of body producing the
     *              magnetic dipole moment.                                                         
     * \param originOfRotation Constant origin of axis about which the magnetic field is rotating    
     * \param axisOfMagneticField Constant unitvector that points in the direction of the 
     *              symmetry axis of magnetic dipole field.                                       
     * \param angularVelocityOfRotatingBody Constant angular velocity of body which producing the 
     *              corotation electric field.                                                    
     * \param chargeOfAcceleratedBody Constant charge of body accelerated by magnetic field.         
     * \param permeabilityOfSpaceBetweenBodies Constant permeability of space between bodies.        
     * \param dipoleMomentOfMagneticField Constant magnetic dipole moment from body exerting acceleration.
     * \param mass Constant mass of body undergoing acceleration.                                
     * \param timeElapsedFunction Function returning current time elapsed since static magnetic field 
     *              was in initial orientation 
     */
    RotatingMagneticDipoleField(
            Vector3dReturningFunction acceleratedBodyPositionFunction,
            Vector3dReturningFunction acceleratedBodyVelocityFunction,
            const Eigen::Vector3d& originOfMagneticField,
            const Eigen::Vector3d& originOfRotation,
            const Eigen::Vector3d& axisOfMagneticField,
            const Eigen::Vector3d& angularVelocityOfRotatingBody,            
            const double chargeOfAcceleratedBody,
            const double permeabilityOfSpaceBetweenBodies,
            const double dipoleMomentOfMagneticField,
            const double mass,
            DoubleReturningFunction timeElapsedFunction )
        : acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          acceleratedBodyVelocityFunction_( acceleratedBodyVelocityFunction ),
          originOfMagneticFieldFunction_( 
                boost::lambda::constant( originOfMagneticField ) ),
          originOfRotationFunction_( 
                boost::lambda::constant( originOfRotation ) ),
          axisOfMagneticFieldFunction_( 
                boost::lambda::constant( axisOfMagneticField ) ),
          angularVelocityOfRotatingBodyFunction_( 
                boost::lambda::constant( angularVelocityOfRotatingBody ) ),
          chargeOfAcceleratedBodyFunction_( 
                boost::lambda::constant( chargeOfAcceleratedBody ) ),
          permeabilityOfSpaceBetweenBodiesFunction_(
                boost::lambda::constant( permeabilityOfSpaceBetweenBodies ) ),
          dipoleMomentOfMagneticFieldFunction_(
                boost::lambda::constant( dipoleMomentOfMagneticField ) ),
          massFunction_( boost::lambda::constant( mass ) ),
          timeElapsedFunction_( timeElapsedFunction )
    {
        this->updateMembers( );
    }

    //! Get rotating magnetic dipole field acceleration.
    /*!
     * Returns the rotating magnetic dipole field acceleration. No arguments are passed to this 
     * function. Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc., which are to be set in a derived class and evaluated by the
     * updateMembers() function below. This function is essentially a wrapper for the free
     * function that computes the rotating magnetic dipole field acceleration.
     * \return rotating magnetic dipole field acceleration.
     * \sa computeRotatingMagneticDipoleFieldAcceleration().
     */
    Eigen::Vector3d getAcceleration( );

    //! Update member variables used by the rotating magnetic dipole field acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration() function.
     */
    void updateMembers( );

private:

    //! Function pointer returning position of accelerated body.
    /*!
     * Function pointer returning position of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    //! Function pointer returning velocity of accelerated body.
    /*!
     * Function pointer returning velocity of accelerated body (3D vector).
     */
    const Vector3dReturningFunction acceleratedBodyVelocityFunction_;    

    //! Function pointer returning origin of magnetic field.
    /*!
     * Function pointer returning origin of magnetic field (3D vector).
     */
    const Vector3dReturningFunction originOfMagneticFieldFunction_;

    //! Function pointer returning origin of rotation.
    /*!
     * Function pointer returning origin of rotation (3D vector).
     */
    const Vector3dReturningFunction originOfRotationFunction_;

    //! Function pointer returning axis of magnetic field.
    /*!
     * Function pointer returning axis of magnetic field (3D vector).
     */
    const Vector3dReturningFunction axisOfMagneticFieldFunction_;

    //! Function pointer returning angular velocity of rotating body.
    /*!
     * Function pointer returning angular velocity of rotating body (3D vector).
     */
    const Vector3dReturningFunction angularVelocityOfRotatingBodyFunction_;

    //! Function pointer returning charge of accelerated body.
    /*!
     * Function pointer returning charge of accelerated body [C].
     */
    const DoubleReturningFunction chargeOfAcceleratedBodyFunction_;

    //! Function pointer returning permeability of space between bodies.
    /*!
     * Function pointer returning permeability of space between bodies [V·s/(A·m)].
     */
    const DoubleReturningFunction permeabilityOfSpaceBetweenBodiesFunction_;

    //! Function pointer returning dipole moment of magnetic field.
    /*!
     * Function pointer returning dipole moment of magnetic field [T·m].
     */
    const DoubleReturningFunction dipoleMomentOfMagneticFieldFunction_;

    //! Function pointer returning mass of accelerated body.
    /*!
     * Function pointer returning mass of accelerated body [kg.
     */
    const DoubleReturningFunction massFunction_;

    //! Function pointer returning time elapsed since initial orientation of magnetic field.
    /*!
     * Function pointer returning time elapsed since initial orientation of magnetic field [s].
     */
    const DoubleReturningFunction timeElapsedFunction_;

    //! Current vector of position of accelerated body.
    /*!
     * Current vector of position of accelerated body (3D vector).
     */
    Eigen::Vector3d currentAcceleratedBodyPosition_;

    //! Current vector of velocity of accelerated body.
    /*!
     * Current vector of velocity of accelerated body (3D vector).
     */
    Eigen::Vector3d currentAcceleratedBodyVelocity_;

    //! Current origin of magnetic field.
    /*!
     * Current origin of magnetic field (3D vector).
     */
    Eigen::Vector3d currentOriginOfMagneticField_;

    //! Current origin of rotation.
    /*!
     * Current origin of rotation (3D vector).
     */
    Eigen::Vector3d currentOriginOfRotation_;

    //! Current axis of magnetic field.
    /*!
     * Current axis of magnetic field (3D vector).
     */
    Eigen::Vector3d currentAxisOfMagneticField_;

    //! Current vector of angular velocity of rotating body.
    /*!
     * Current vector of angular velocity of rotating body (3D vector).
     */
    Eigen::Vector3d currentAngularVelocityOfRotatingBody_;

    //! Current charge of accelerated body.
    /*!
     * Current charge of accelerated body [C].
     */
    double currentChargeOfAcceleratedBody_;

    //! Current permeability of space between bodies.
    /*!
     * Current permeability of space between bodies [V·s/(A·m)].
     */
    double currentPermeabilityOfSpaceBetweenBodies_;


    //! Current dipole moment of magnetic field.
    /*!
     * Current dipole moment of magnetic field [T·m].
     */
    double currentDipoleMomentOfMagneticField_;


    //! Current mass of accelerated body.
    /*!
     * Current mass of accelerated body [kg].
     */
    double currentMass_;

    //! Current time elapsed since initial orientation of magnetic field.
    /*!
     * Current time elapsed since initial orientation of magnetic field [s].
     */
    double currentTimeElapsed_;

};


} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_ROTATING_MAGNETIC_DIPOLE_ACCELERATION_H
