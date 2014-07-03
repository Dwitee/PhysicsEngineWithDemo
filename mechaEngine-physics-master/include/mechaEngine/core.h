/*
 * Interface file for math helperutils
 *
 * Copyright (c) Dwitee Krishna Panda. All Rights Reserved.
 *
*/

#include <math.h>

/**
 * basic mathematical helpers and utils.
 */
#ifndef MECHAENGINE_CORE_H
#define MECHAENGINE_CORE_H

#include "precision.h"


namespace mechaEngine {

    
    //value of kinetic energy under which the body will be put to sleep.
    //default value is 0.1
    
    extern real minValForSleep;

    
    //the kinetic energy under
    //which a body may be put to sleep.
    //The value is global; all bodies will use it.
    
    real getminValForSleep();
    void setminValForSleep(real value);

    
    //Holds a vector in 3D. Four data members are allocated
    //to ensure alignment in an array.
    class Vector3
    {
    public:
        //value along the x axis.
        real x;

        /** value along the y axis. */
        real y;

        /** value along the z axis. */
        real z;

    private:
        
        /** Padding to ensure 4 word alignment. */
        real pad;

    public:
        /** The default constructor . */
        Vector3() : x(0), y(0), z(0) {}

        /** The explicit constructor creates a vector with given values */
        Vector3(const real x, const real y, const real z)
            : x(x), y(y), z(z) {}

        const static Vector3 UP;
        const static Vector3 RIGHT;
        const static Vector3 GRAVITY;
        const static Vector3 HIGH_GRAVITY;
		
		const static Vector3 X;
        const static Vector3 Y;
        const static Vector3 Z;

        /** Adds thw passed vector to this vector */
        void operator+=(const Vector3& vec)
        {
            this->x += vec.x;
            this->y += vec.y;
            this->z += vec.z;
        }
        
        /** for v -= v1*/
        void operator-=(const Vector3& v)
        {
            this->x -= v.x;
            this->y -= v.y;
            this->z -= v.z;
        }
        
        /** overloading + operator. */
        Vector3 operator+(const Vector3& vec) const
        {
            return Vector3(x+vec.x, y+vec.y, z+vec.z);
        }

        /** get vector components by index */
        real& operator[](unsigned i)
        {
            if (i == 0)
            {
                return x;
            }
            if (i == 1)
            {
                return y;
            }
            return z;
        }
        
        /** get vector components by index */
        real operator[](unsigned i) const
        {
            if (i == 0)
            {
                return x;
            }
            if (i == 1)
            {
                return y;
            }
            return z;
        }
        
        /**  for vR = v -v1 */
        Vector3 operator-(const Vector3& v) const
        {
            return Vector3(x-v.x, y-v.y, z-v.z);
        }

        Vector3 operator*(const real scaleFactor) const
        {
            return Vector3( x*scaleFactor, y*scaleFactor, z*scaleFactor);
        }

        /** Multiplies this vector by the passed scalar. */
        void operator*=(const real scaleFactor)
        {
            this->x *= scaleFactor;
            this->y *= scaleFactor;
            this->z *= scaleFactor;
        }

        void componentProductUpdate(const Vector3 &vector)
        {
            this->x *= vector.x;
            this->y *= vector.y;
            this->z *= vector.z;
        }
		
        Vector3 componentProduct(const Vector3 &vector) const
        {
            return Vector3(x * vector.x, y * vector.y, z * vector.z);
        }



        Vector3 vectorProduct(const Vector3 &vector)  const
        {
	        return Vector3(y*vector.z-z*vector.y,
	        z*vector.x-x*vector.z,
	        x*vector.y-y*vector.x);
        }

        Vector3 operator%(const Vector3 &vector) const
        {
            return Vector3(y*vector.z-z*vector.y,
                           z*vector.x-x*vector.z,
                           x*vector.y-y*vector.x);
        }
		
        void operator %=(const Vector3 &vector)
        {
            *this = vectorProduct(vector);
        }



        real operator *(const Vector3 &vector) const
        {
            return x*vector.x + y*vector.y + z*vector.z;
        }
		
        real scalarProduct(const Vector3 &vector) const
        {
            return x*vector.x + y*vector.y + z*vector.z;
        }



        void addScaledVector(const Vector3& vector, real scale)
        {
            x += vector.x * scale;
            y += vector.y * scale;
            z += vector.z * scale;
        }

        /** size is max limit */
        void trim(real size)
        {
            if (squareMagnitude() > size*size)
            {
                normalise();
                x *= size;
                y *= size;
                z *= size;
            }
        }
        real magnitude() const
        {
            return real_sqrt(x*x+y*y+z*z);
        }
        
        real squareMagnitude() const
        {
            return x*x + y*y + z*z;
        }



        /** converts a non zero vector into a unit vector. */
        void normalise()
        {
            real l = magnitude();
            if (l > 0)
            {
                (*this) *= ((real)1)/l;
            }
        }

        /** normalise a vector. */
        Vector3 unit()  const
        {
            Vector3 resultVec = *this;
            resultVec.normalise();
            return resultVec;
        }

        /** Checks if the two vectors are identical . */
        bool operator==(const Vector3& otherVec) const
        {
            return  x == otherVec.x &&
                    y == otherVec.y &&
                    z == otherVec.z;
        }

        /** Checks if the two vectors are non-identical . */
        bool operator!=(const Vector3& otherVec) const
        {
             return !(*this == otherVec);
        }
        
        /**
         * Checks if this vector is less than or equal to the other.
         */
        bool operator<=(const Vector3& otherVec) const
        {
            return x <= otherVec.x && y <= otherVec.y && z <= otherVec.z;
        }
        
        /**
         * Checks if this vector is greater or equal to the other.
         */
        bool operator>=(const Vector3& otherVec) const
        {
            return x >= otherVec.x && y >= otherVec.y && z >= otherVec.z;
        }
        /**
         * Checks if this vector is greater than the other.
         */
        bool operator>(const Vector3& otherVec) const
        {
            return x > otherVec.x && y > otherVec.y && z > otherVec.z;
        }
        
        /** Checks if this vector is  less than the other.*/
        bool operator<(const Vector3& otherVec) const
        {
            return x < otherVec.x && y < otherVec.y && z < otherVec.z;
        }

        /** Flips the vector. */
        void invert()
        {
            x =  -x;
            y =  -y;
            z =  -z;
        }

        /** Zeros the vector. */
        void clear()
        {
            x = y = z = 0;
        }



    };

    /**
     * Holds a three degree of freedom orientation.
     * Quaternions are used to representation orientation, it require four  data to
     * hold  three degrees of freedom. These four items of data are coefficients of a complex number with three
     * imaginary parts. Useful for 3D rotations.
     *A quaternion is valid if its  length is 1.
     */
    class Quaternion
    {
    public:
        union
        {
            struct
            {
                /**first complex component of the quaternion.*/
                real i;
                
                /**second complex component of the quaternion.*/
                real j;
                
                /** third complex component of the quaternion.*/
                real k;
                
                /**real component of the quaternion.*/
                real r;
            };
            
            /**
             * Holds the quaternion data in array form.
             */
            real data[4];
        };
        
        
        
        /**
         * The default constructor creates a quaternion with zero rotation.
         */
        Quaternion() : r(1), i(0), j(0), k(0) {}
        
        Quaternion(const real r, const real i, const real j, const real k)
        : r(r), i(i), j(j), k(k)
        {
        }
        
        void rotateByVector(const Vector3& vector)
        {
            Quaternion q(0, vector.x, vector.y, vector.z);
            (*this) *= q;
        }
        
        void normalise()
        {
            real d = r*r + i*i + j*j + k*k;
            
            if (d < real_epsilon) {
                r = 1;
                return;
            }
            
            d = ((real)1.0)/real_sqrt(d);
            r *= d;
            i *= d;
            j *= d;
            k *= d;
        }
        
        
        void operator *=(const Quaternion &multFactor)
        {
            Quaternion quat = *this;
            
            r = quat.r*multFactor.r - quat.i*multFactor.i -
            quat.j*multFactor.j - quat.k*multFactor.k;
            
            i = quat.r*multFactor.i + quat.i*multFactor.r +
            quat.j*multFactor.k - quat.k*multFactor.j;
            
            j = quat.r*multFactor.j + quat.j*multFactor.r +
            quat.k*multFactor.i - quat.i*multFactor.k;
            
            k = quat.r*multFactor.k + quat.k*multFactor.r +
            quat.i*multFactor.j - quat.j*multFactor.i;
        }
        
        /**
         * This is used to update the orientation quaternion by a rotation
         * and time.
         */
        void addScaledVector(const Vector3& vector, real scale)
        {
            Quaternion quat(0,
                            vector.x * scale,
                            vector.y * scale,
                            vector.z * scale);
            quat *= *this;
            r += quat.r * ((real)0.5);
            i += quat.i * ((real)0.5);
            j += quat.j * ((real)0.5);
            k += quat.k * ((real)0.5);
        }
        
    };
    
    /**
     homogenous 4x4 matrix. used for transformation
     */
    class Matrix3
    {
    public:
        /**
         * Holds the tensor matrix data in array form.
         */
        real data[9];
        
        // ... Other Matrix3 code as before ...
        
        /**
         * Creates a new matrix.
         */
        Matrix3()
        {
            data[0] = data[1] = data[2] = data[3] = data[4] = data[5] =
            data[6] = data[7] = data[8] = 0;
        }
        
        /**
         * Creates a new matrix with the given three vectors making
         * up its columns.
         */
        Matrix3(const Vector3 &compOne, const Vector3 &compTwo,
                const Vector3 &compThree)
        {
            setComponents(compOne, compTwo, compThree);
        }
        
        /**
         * Creates a new matrix with explicit coefficients.
         */
        Matrix3(real c0, real c1, real c2, real c3, real c4, real c5,
                real c6, real c7, real c8)
        {
            data[0] = c0; data[1] = c1; data[2] = c2;
            data[3] = c3; data[4] = c4; data[5] = c5;
            data[6] = c6; data[7] = c7; data[8] = c8;
        }
        
        /**
         * Sets the matrix to be a diagonal matrix with the given
         * values along the leading diagonal.
         */
        void setDiagonal(real a, real b, real c)
        {
            setInertiaTensorCoeffs(a, b, c);
        }
        
        /**
         * Sets the value of the matrix from inertia tensor values.
         */
        void setInertiaTensorCoeffs(real ix, real iy, real iz,
                                    real ixy=0, real ixz=0, real iyz=0)
        {
            data[0] = ix;
            data[1] = data[3] = -ixy;
            data[2] = data[6] = -ixz;
            data[4] = iy;
            data[5] = data[7] = -iyz;
            data[8] = iz;
        }
        
        /**
         * Sets the value of the matrix as an inertia tensor of
         * a rectangular block aligned with the body's coordinate
         * system with the given axis half-sizes and mass.
         */
        void setBlockInertiaTensor(const Vector3 &halfSizes, real mass)
        {
            Vector3 squares = halfSizes.componentProduct(halfSizes);
            setInertiaTensorCoeffs(0.3f*mass*(squares.y + squares.z),
                                   0.3f*mass*(squares.x + squares.z),
                                   0.3f*mass*(squares.x + squares.y));
        }
        
        /**
         * Sets the matrix to be a skew symmetric matrix based on
         * the given vector. The skew symmetric matrix is the equivalent
         * of the vector product. So if a,b are vectors. a x b = A_s b
         * where A_s is the skew symmetric form of a.
         */
        void setSkewSymmetric(const Vector3 vector)
        {
            data[0] = data[4] = data[8] = 0;
            data[1] = -vector.z;
            data[2] = vector.y;
            data[3] = vector.z;
            data[5] = -vector.x;
            data[6] = -vector.y;
            data[7] = vector.x;
        }
        
        /**
         * Sets the matrix values from the given three vector components.
         * These are arranged as the three columns of the vector.
         */
        void setComponents(const Vector3 &compOne, const Vector3 &compTwo,
                           const Vector3 &compThree)
        {
            data[0] = compOne.x;
            data[1] = compTwo.x;
            data[2] = compThree.x;
            data[3] = compOne.y;
            data[4] = compTwo.y;
            data[5] = compThree.y;
            data[6] = compOne.z;
            data[7] = compTwo.z;
            data[8] = compThree.z;
            
        }
        
        /**
         * Transform the given vector by this matrix.
         *
         * @param vector The vector to transform.
         */
        Vector3 operator*(const Vector3 &vector) const
        {
            return Vector3(
                           vector.x * data[0] + vector.y * data[1] + vector.z * data[2],
                           vector.x * data[3] + vector.y * data[4] + vector.z * data[5],
                           vector.x * data[6] + vector.y * data[7] + vector.z * data[8]
                           );
        }
        
        /**
         * Transform the given vector by this matrix.
         *
         * @param vector The vector to transform.
         */
        Vector3 transform(const Vector3 &vector) const
        {
            return (*this) * vector;
        }
        
        /**
         * Transform the given vector by the transpose of this matrix.
         *
         * @param vector The vector to transform.
         */
        Vector3 transformTranspose(const Vector3 &vector) const
        {
            return Vector3(
                           vector.x * data[0] + vector.y * data[3] + vector.z * data[6],
                           vector.x * data[1] + vector.y * data[4] + vector.z * data[7],
                           vector.x * data[2] + vector.y * data[5] + vector.z * data[8]
                           );
        }
        
        /**
         * Gets a vector representing one row in the matrix.
         *
         * @param i The row to return.
         */
        Vector3 getRowVector(int i) const
        {
            return Vector3(data[i*3], data[i*3+1], data[i*3+2]);
        }
        
        /**
         * Gets a vector representing one axis (i.e. one column) in the matrix.
         *
         * @param i The row to return.
         *
         * @return The vector.
         */
        Vector3 getAxisVector(int i) const
        {
            return Vector3(data[i], data[i+3], data[i+6]);
        }
        
        /**
         * Sets the matrix to be the inverse of the given matrix.
         *
         * @param m The matrix to invert and use to set this.
         */
        void setInverse(const Matrix3 &m)
        {
            real t4 = m.data[0]*m.data[4];
            real t6 = m.data[0]*m.data[5];
            real t8 = m.data[1]*m.data[3];
            real t10 = m.data[2]*m.data[3];
            real t12 = m.data[1]*m.data[6];
            real t14 = m.data[2]*m.data[6];
            
            // Calculate the determinant
            real t16 = (t4*m.data[8] - t6*m.data[7] - t8*m.data[8]+
                        t10*m.data[7] + t12*m.data[5] - t14*m.data[4]);
            
            // Make sure the determinant is non-zero.
            if (t16 == (real)0.0f) return;
            real t17 = 1/t16;
            
            data[0] = (m.data[4]*m.data[8]-m.data[5]*m.data[7])*t17;
            data[1] = -(m.data[1]*m.data[8]-m.data[2]*m.data[7])*t17;
            data[2] = (m.data[1]*m.data[5]-m.data[2]*m.data[4])*t17;
            data[3] = -(m.data[3]*m.data[8]-m.data[5]*m.data[6])*t17;
            data[4] = (m.data[0]*m.data[8]-t14)*t17;
            data[5] = -(t6-t10)*t17;
            data[6] = (m.data[3]*m.data[7]-m.data[4]*m.data[6])*t17;
            data[7] = -(m.data[0]*m.data[7]-t12)*t17;
            data[8] = (t4-t8)*t17;
        }
        
        /** Returns a new matrix containing the inverse of this matrix. */
        Matrix3 inverse() const
        {
            Matrix3 result;
            result.setInverse(*this);
            return result;
        }
        
        /**
         * Inverts the matrix.
         */
        void invert()
        {
            setInverse(*this);
        }
        
        /**
         * Sets the matrix to be the transpose of the given matrix.
         *
         * @param m The matrix to transpose and use to set this.
         */
        void setTranspose(const Matrix3 &m)
        {
            data[0] = m.data[0];
            data[1] = m.data[3];
            data[2] = m.data[6];
            data[3] = m.data[1];
            data[4] = m.data[4];
            data[5] = m.data[7];
            data[6] = m.data[2];
            data[7] = m.data[5];
            data[8] = m.data[8];
        }
        
        /** Returns a new matrix containing the transpose of this matrix. */
        Matrix3 transpose() const
        {
            Matrix3 result;
            result.setTranspose(*this);
            return result;
        }
        
        /**
         * Returns a matrix which is this matrix multiplied by the given
         * other matrix.
         */
        Matrix3 operator*(const Matrix3 &otherMatrix) const
        {
            return Matrix3(
                           data[0]*otherMatrix.data[0] + data[1]*otherMatrix.data[3] + data[2]*otherMatrix.data[6],
                           data[0]*otherMatrix.data[1] + data[1]*otherMatrix.data[4] + data[2]*otherMatrix.data[7],
                           data[0]*otherMatrix.data[2] + data[1]*otherMatrix.data[5] + data[2]*otherMatrix.data[8],
                           
                           data[3]*otherMatrix.data[0] + data[4]*otherMatrix.data[3] + data[5]*otherMatrix.data[6],
                           data[3]*otherMatrix.data[1] + data[4]*otherMatrix.data[4] + data[5]*otherMatrix.data[7],
                           data[3]*otherMatrix.data[2] + data[4]*otherMatrix.data[5] + data[5]*otherMatrix.data[8],
                           
                           data[6]*otherMatrix.data[0] + data[7]*otherMatrix.data[3] + data[8]*otherMatrix.data[6],
                           data[6]*otherMatrix.data[1] + data[7]*otherMatrix.data[4] + data[8]*otherMatrix.data[7],
                           data[6]*otherMatrix.data[2] + data[7]*otherMatrix.data[5] + data[8]*otherMatrix.data[8]
                           );
        }
        
        /**
         * Multiplies this matrix in place by the given other matrix.
         */
        void operator*=(const Matrix3 &otherMatrix)
        {
            real t1;
            real t2;
            real t3;
            
            t1 = data[0]*otherMatrix.data[0] + data[1]*otherMatrix.data[3] + data[2]*otherMatrix.data[6];
            t2 = data[0]*otherMatrix.data[1] + data[1]*otherMatrix.data[4] + data[2]*otherMatrix.data[7];
            t3 = data[0]*otherMatrix.data[2] + data[1]*otherMatrix.data[5] + data[2]*otherMatrix.data[8];
            data[0] = t1;
            data[1] = t2;
            data[2] = t3;
            
            t1 = data[3]*otherMatrix.data[0] + data[4]*otherMatrix.data[3] + data[5]*otherMatrix.data[6];
            t2 = data[3]*otherMatrix.data[1] + data[4]*otherMatrix.data[4] + data[5]*otherMatrix.data[7];
            t3 = data[3]*otherMatrix.data[2] + data[4]*otherMatrix.data[5] + data[5]*otherMatrix.data[8];
            data[3] = t1;
            data[4] = t2;
            data[5] = t3;
            
            t1 = data[6]*otherMatrix.data[0] + data[7]*otherMatrix.data[3] + data[8]*otherMatrix.data[6];
            t2 = data[6]*otherMatrix.data[1] + data[7]*otherMatrix.data[4] + data[8]*otherMatrix.data[7];
            t3 = data[6]*otherMatrix.data[2] + data[7]*otherMatrix.data[5] + data[8]*otherMatrix.data[8];
            data[6] = t1;
            data[7] = t2;
            data[8] = t3;
        }
        
        /**
         * Multiplies this matrix in place by the given scalar.
         */
        void operator*=(const real scalar)
        {
            data[0] *= scalar; data[1] *= scalar; data[2] *= scalar;
            data[3] *= scalar; data[4] *= scalar; data[5] *= scalar;
            data[6] *= scalar; data[7] *= scalar; data[8] *= scalar;
        }
        
        /**
         * Does a component-wise addition of this matrix and the given
         * matrix.
         */
        void operator+=(const Matrix3 &otherMatrix)
        {
            data[0] += otherMatrix.data[0]; data[1] += otherMatrix.data[1]; data[2] += otherMatrix.data[2];
            data[3] += otherMatrix.data[3]; data[4] += otherMatrix.data[4]; data[5] += otherMatrix.data[5];
            data[6] += otherMatrix.data[6]; data[7] += otherMatrix.data[7]; data[8] += otherMatrix.data[8];
        }
        
        /**
         * Sets this matrix to be the rotation matrix corresponding to
         * the given quaternion.
         */
        void setOrientation(const Quaternion &quat)
        {
            data[0] = 1 - (2*quat.j*quat.j + 2*quat.k*quat.k);
            data[1] = 2*quat.i*quat.j + 2*quat.k*quat.r;
            data[2] = 2*quat.i*quat.k - 2*quat.j*quat.r;
            data[3] = 2*quat.i*quat.j - 2*quat.k*quat.r;
            data[4] = 1 - (2*quat.i*quat.i  + 2*quat.k*quat.k);
            data[5] = 2*quat.j*quat.k + 2*quat.i*quat.r;
            data[6] = 2*quat.i*quat.k + 2*quat.j*quat.r;
            data[7] = 2*quat.j*quat.k - 2*quat.i*quat.r;
            data[8] = 1 - (2*quat.i*quat.i  + 2*quat.j*quat.j);
        }
        
        /**
         * Interpolates a couple of matrices.
         */
        static Matrix3 linearInterpolate(const Matrix3& a, const Matrix3& b, real prop);
    };

    /**
    * homogenous 4x4 matrix. used for transformation
    */
    class Matrix4
    {
    public:
        /** Holds the transform matrix data in array form. */
        real data[12];

        /** sets diagonal matrix with the given coefficients. */
        void setDiagonal(real a, real b, real c)
        {
            data[0] = a;
            data[5] = b;
            data[10] = c;
        }

        /** makes identity matrix. */
        Matrix4()
        {
            data[1] = data[2] = data[3] = data[4] = data[6] = data[7] = data[8] = data[9] = data[11] = 0;
            data[0] = data[5] = data[10] = 1;
        }

        /** Returns a matrix which is this matrix multiplied by the given other matrix.*/
        Matrix4 operator*(const Matrix4 &otherMatrix) const
        {
            Matrix4 result;
            result.data[0] = (otherMatrix.data[0]*data[0]) + (otherMatrix.data[4]*data[1]) + (otherMatrix.data[8]*data[2]);
            result.data[4] = (otherMatrix.data[0]*data[4]) + (otherMatrix.data[4]*data[5]) + (otherMatrix.data[8]*data[6]);
            result.data[8] = (otherMatrix.data[0]*data[8]) + (otherMatrix.data[4]*data[9]) + (otherMatrix.data[8]*data[10]);

            result.data[1] = (otherMatrix.data[1]*data[0]) + (otherMatrix.data[5]*data[1]) + (otherMatrix.data[9]*data[2]);
            result.data[5] = (otherMatrix.data[1]*data[4]) + (otherMatrix.data[5]*data[5]) + (otherMatrix.data[9]*data[6]);
            result.data[9] = (otherMatrix.data[1]*data[8]) + (otherMatrix.data[5]*data[9]) + (otherMatrix.data[9]*data[10]);

            result.data[2] = (otherMatrix.data[2]*data[0]) + (otherMatrix.data[6]*data[1]) + (otherMatrix.data[10]*data[2]);
            result.data[6] = (otherMatrix.data[2]*data[4]) + (otherMatrix.data[6]*data[5]) + (otherMatrix.data[10]*data[6]);
            result.data[10] = (otherMatrix.data[2]*data[8]) + (otherMatrix.data[6]*data[9]) + (otherMatrix.data[10]*data[1]);

            result.data[3] = (otherMatrix.data[3]*data[0]) + (otherMatrix.data[7]*data[1]) + (otherMatrix.data[11]*data[2])
            + data[3];
            
            result.data[7] = (otherMatrix.data[3]*data[4]) + (otherMatrix.data[7]*data[5]) + (otherMatrix.data[11]*data[6])
            + data[7];
            
            result.data[11] = (otherMatrix.data[3]*data[8]) + (otherMatrix.data[7]*data[9]) + (otherMatrix.data[11]*data[10])
            + data[11];

            return result;
        }

        /** Transform the given vector by this matrix */
        Vector3 operator*(const Vector3 &vector) const
        {
            return Vector3(
            vector.x * data[0] +
            vector.y * data[1] +
            vector.z * data[2] + data[3],

            vector.x * data[4] +
            vector.y * data[5] +
            vector.z * data[6] + data[7],

            vector.x * data[8] +
            vector.y * data[9] +
            vector.z * data[10] + data[11]
            );
        }

        /** Transform the given vector by this matrix.*/
        Vector3 transform(const Vector3 &vector) const
        {
             return (*this) * vector;
        }
        
        /**Sets the matrix to be the inverse of the given matrix.*/
        void setInverse(const Matrix4 &m);
        
        /** Returns a new matrix containing the inverse of this matrix. */
        Matrix4 inverse() const
        {
            Matrix4 result;
            result.setInverse(*this);
            return result;
        }
        
        /**
         * Inverts the matrix.
         */
        void invert()
        {
            setInverse(*this);
        }
        
        /** Returns the determinant of the matrix. */
        real getDeterminant() const;

        /**
         * Transform the given direction vector by this matrix.
         */
        Vector3 transformDirection(const Vector3 &vector) const
        {
            return Vector3(
                vector.x * data[0] +
                vector.y * data[1] +
                vector.z * data[2],

                vector.x * data[4] +
                vector.y * data[5] +
                vector.z * data[6],

                vector.x * data[8] +
                vector.y * data[9] +
                vector.z * data[10]
            );
        }

        /**
         * Transforms the given direction vector by the
         * transformational inverse of this matrix.
         */
        Vector3 transformInverseDirection(const Vector3 &vector) const
        {
            return Vector3(
                vector.x * data[0] +
                vector.y * data[4] +
                vector.z * data[8],

                vector.x * data[1] +
                vector.y * data[5] +
                vector.z * data[9],

                vector.x * data[2] +
                vector.y * data[6] +
                vector.z * data[10]
            );
        }

        /**
         * Transform the given vector by the transformational inverse
         * of this matrix.
         */
        Vector3 transformInverse(const Vector3 &vector) const
        {
            Vector3 tmpVec = vector;
            tmpVec.x -= data[3];
            tmpVec.y -= data[7];
            tmpVec.z -= data[11];
            
            return Vector3(
                tmpVec.x * data[0] +
                tmpVec.y * data[4] +
                tmpVec.z * data[8],

                tmpVec.x * data[1] +
                tmpVec.y * data[5] +
                tmpVec.z * data[9],

                tmpVec.x * data[2] +
                tmpVec.y * data[6] +
                tmpVec.z * data[10]
            );
        }

        /**Gets a vector representing one axis in the matrix.*/
        Vector3 getAxisVector(int i) const
        {
            return Vector3(data[i], data[i+4], data[i+8]);
        }

        /**Sets this matrix to be the rotation matrix corresponding to the given quaternion.*/
        void setOrientationAndPos(const Quaternion &quat, const Vector3 &posVec)
        {
            data[0] = 1 - (2*quat.j*quat.j + 2*quat.k*quat.k);
            data[1] = 2*quat.i*quat.j + 2*quat.k*quat.r;
            data[2] = 2*quat.i*quat.k - 2*quat.j*quat.r;
            data[3] = posVec.x;

            data[4] = 2*quat.i*quat.j - 2*quat.k*quat.r;
            data[5] = 1 - (2*quat.i*quat.i  + 2*quat.k*quat.k);
            data[6] = 2*quat.j*quat.k + 2*quat.i*quat.r;
            data[7] = posVec.y;

            data[8] = 2*quat.i*quat.k + 2*quat.j*quat.r;
            data[9] = 2*quat.j*quat.k - 2*quat.i*quat.r;
            data[10] = 1 - (2*quat.i*quat.i  + 2*quat.j*quat.j);
            data[11] = posVec.z;
        }

        /**
         * Fills the given array with this transform matrix, so it is
         * usable as an open-gl transform matrix. OpenGL uses a column
         * major format, so that the values are transposed as they are
         * written.
         */
        void fillGLArray(float array[16]) const
        {
            array[0] = (float)data[0];
            array[1] = (float)data[4];
            array[2] = (float)data[8];
            array[3] = (float)0;

            array[4] = (float)data[1];
            array[5] = (float)data[5];
            array[6] = (float)data[9];
            array[7] = (float)0;

            array[8] = (float)data[2];
            array[9] = (float)data[6];
            array[10] = (float)data[10];
            array[11] = (float)0;

            array[12] = (float)data[3];
            array[13] = (float)data[7];
            array[14] = (float)data[11];
            array[15] = (float)1;
        }
    };


    


}

#endif // MECHAENGINE_CORE_H