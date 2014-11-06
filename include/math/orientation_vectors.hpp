#ifndef MATH_ORIENTATION_VECTORS_HPP_INCLUDED
#define MATH_ORIENTATION_VECTORS_HPP_INCLUDED

#include "math/vector3.hpp"

namespace aspect { namespace math {

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

}} // aspect::math

#endif // MATH_ORIENTATION_VECTORS_HPP_INCLUDED
