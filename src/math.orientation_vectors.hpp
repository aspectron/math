#pragma once
#ifndef __MATH_ORIENTATION_VECTORS_HPP__
#define __MATH_ORIENTATION_VECTORS_HPP__

#include "math.vector3.hpp"

namespace aspect
{

namespace math
{

class vec3;
class quat;
class matrix;

//============================================================================
// Orienation vectors
//
class MATH_API orientation_vectors
{
	public:

		vec3 m_locO;
		vec3 m_locX;
		vec3 m_locY;

	public:

		orientation_vectors() {}
		orientation_vectors(const quat& q) { from_quaternion(q); }

		void from_quaternion(const quat& q);
		void to_quaternion(quat& q) const;
		void apply_matrix(matrix &m);
};

} // math

} // aspect

#endif // MATH_ORIENTATION_VECTORS_HPP__