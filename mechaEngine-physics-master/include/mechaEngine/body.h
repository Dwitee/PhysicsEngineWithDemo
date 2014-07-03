/*
 * Interface file for the rigid body class.
 *
 * Part of the MechaEngine physics system.
 *
 * Copyright (c) Dwitee Krishna Panda 2012. All Rights Reserved.
 */

/**
 * @file
 *
 * This file contains the definitions for the rigid body class, the
 * basic building block of all the physics system.
 */
#ifndef MECHAENGINE_BODY_H
#define MECHAENGINE_BODY_H

#include "core.h"


namespace mechaEngine {

    
    class RigidBody
    {

    protected:

        real inverseMass;
        real linearDamping;
        real angularDamping;
        Matrix3 inverseInertiaTensor;
        Matrix3 inverseInertiaTensorWorld;
        Vector3 position;
        Quaternion orientation;
        Vector3 velocity;
        Vector3 rotation;
        real motion; //holds the ammount of energy with the body this can decide if the body can be put to sleep
        bool isAwake;
        bool canSleep;
        Matrix4 transformMatrix; // to conver the body from object space to world space
        Vector3 forceAccum; // resultant of all the forces acting on the body
        Vector3 torqueAccum;// resultant of all the torques acting on the body
        Vector3 acceleration; // sum of all the acceleration acting on the body

        Vector3 lastFrameAcceleration;// Holds the linear acceleration of the rigid body, for the body

    public:

        void calculateDerivedData();

        // This function integrates the  body forward in time by the given amount.
        // It uses a Newton-Euler integration method, which is a
        // linear approximation to the correct integral.
        void integrate(real duration);

        void setMass(const real mass);

        real getMass() const;

        void setInverseMass(const real inverseMass);

        real getInverseMass() const;

        bool hasFiniteMass() const;

        void setInertiaTensor(const Matrix3 &inertiaTensor);


        void getInertiaTensor(Matrix3 *inertiaTensor) const;

        Matrix3 getInertiaTensor() const;


        void getInertiaTensorWorld(Matrix3 *inertiaTensor) const;


        Matrix3 getInertiaTensorWorld() const;


        void setInverseInertiaTensor(const Matrix3 &inverseInertiaTensor);

        void getInverseInertiaTensor(Matrix3 *inverseInertiaTensor) const;


        Matrix3 getInverseInertiaTensor() const;

  
        void getInverseInertiaTensorWorld(Matrix3 *inverseInertiaTensor) const;

      
        Matrix3 getInverseInertiaTensorWorld() const;

      
        void setDamping(const real linearDamping, const real angularDamping);

      
        void setLinearDamping(const real linearDamping);

      
        real getLinearDamping() const;

      
        void setAngularDamping(const real angularDamping);

      
        real getAngularDamping() const;

      
        void setPosition(const Vector3 &position);

       
        void setPosition(const real x, const real y, const real z);

       
        void getPosition(Vector3 *position) const;

       
        Vector3 getPosition() const;

       
        void setOrientation(const Quaternion &orientation);

       
        void setOrientation(const real r, const real i,
            const real j, const real k);

       
        void getOrientation(Quaternion *orientation) const;

       
        Quaternion getOrientation() const;

       
        void getOrientation(Matrix3 *matrix) const;

       
        void getOrientation(real matrix[9]) const;

       
        void getTransform(Matrix4 *transform) const;

       
        void getTransform(real matrix[16]) const;

       
        void getGLTransform(float matrix[16]) const;

       
        Matrix4 getTransform() const;

       
        Vector3 getPointInLocalSpace(const Vector3 &point) const;

       
        Vector3 getPointInWorldSpace(const Vector3 &point) const;

       
        Vector3 getDirectionInLocalSpace(const Vector3 &direction) const;

       
        Vector3 getDirectionInWorldSpace(const Vector3 &direction) const;

       
        void setVelocity(const Vector3 &velocity);

       
        void setVelocity(const real x, const real y, const real z);

       
        void getVelocity(Vector3 *velocity) const;

       
        Vector3 getVelocity() const;

        void addVelocity(const Vector3 &deltaVelocity);

        void setRotation(const Vector3 &rotation);

        void setRotation(const real x, const real y, const real z);

        
        void getRotation(Vector3 *rotation) const;

        
        Vector3 getRotation() const;

        
        void addRotation(const Vector3 &deltaRotation);

        
        bool getAwake() const
        {
            return isAwake;
        }

        
        void setAwake(const bool awake=true);

        
        bool getCanSleep() const
        {
            return canSleep;
        }

        
        void setCanSleep(const bool canSleep=true);

        
        void getLastFrameAcceleration(Vector3 *linearAcceleration) const;

        Vector3 getLastFrameAcceleration() const;

        void clearAccumulators();

        void addForce(const Vector3 &force);

        
        void addForceAtPoint(const Vector3 &force, const Vector3 &point);

        
        void addForceAtBodyPoint(const Vector3 &force, const Vector3 &point);

        
        void addTorque(const Vector3 &torque);

        
        void setAcceleration(const Vector3 &acceleration);

        
        void setAcceleration(const real x, const real y, const real z);

        
        void getAcceleration(Vector3 *acceleration) const;

        
        Vector3 getAcceleration() const;


    };

} // namespace mechaEngine

#endif // MECHAENGINE_BODY_H