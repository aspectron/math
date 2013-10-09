#pragma once
#ifndef __MATH_QUATERNION_HPP__
#define __MATH_QUATERNION_HPP__

#include <string.h>

#include "math.hpp"

namespace aspect
{

namespace math
{

class matrix;
class euler_angles;

class MATH_API quat
{
	public:

		union
		{
			struct
			{
				double x, y, z, w;
			};

			double v[4];
		};

	public:

		enum ORD		{ X, Y, Z, W };

		// constructors...

		quat()					{ unity(); }
		quat(const quat& q)	{ *this = q; }
		quat(const matrix& m)	{ from_matrix(m); }
		quat(const euler_angles &ang);
		quat(double x, double y, double z, double w)
			: x(x), y(y), z(z), w(w)
		{
		}

		// accessors...
		
		double	get(ORD ord) const			{ return v[ord]; }
		void	set(ORD ord, double f)	{ v[ord] = f; }

		// operators...

		inline bool operator == (const quat& q) const;
		operator double*() { return(v); }

		// misc functions...

		// Multiply by scalar f. 
		void scale(double f);

		// Caculate conjugate of quaternion.
		void conjugate()	{ v[X]=-v[X]; v[Y]=-v[Y]; v[Z]=-v[Z]; }

		// Place on the unit sphere
		void unity()		{ ::memset(this, 0, sizeof(quat)); v[W]=1.0f; }
		void normalize();

		// Quaternion product this * qR.  Note: order is important!
		// To combine rotations, use the product qFirst.Multiply(qSecond),
		// which gives the effect of rotating by qFirst then qSecond.
		void multiply(const quat &qR);
		inline quat & operator *= (const quat &Q);
		inline quat   operator *  (const quat &Q) const;

		// Matrix conversion
		void to_matrix(matrix &m) const;
		matrix to_matrix() const;
		void from_matrix(const matrix &m);

		// Spherical linear interpolation of unit quaternions with spins.  A & B -> this
		void slerp(double alpha, const quat &a, const quat &b, int spin = 0);
		
};

inline bool quat::operator == (const quat& q) const 
{ 
	return is_zero(q.v[X] - v[X]) 
		&& is_zero(q.v[Y] - v[Y]) 
		&& is_zero(q.v[Z] - v[Z]) 
		&& is_zero(q.v[W] - v[W]); 
}

inline quat & quat::operator *= (const quat &Q) 
{ 
	multiply(Q); 
	return *this; 
}

inline quat   quat::operator *  (const quat &Q) const 
{ 
	quat qr(*this); 
	qr.multiply(Q); 
	return qr; 
}

} // math

} // aspect

#endif // __MATH_QUATERNION_HPP__
