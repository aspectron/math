//#include "aspect.hpp"


#include "math.hpp"
// 
// /*
// #include "math.point2d.hpp"
// #include "math.vec3.hpp"
// #include "math.point4d.hpp"
// */
// 
// #include "math.orientation_vectors.hpp"
// #include "math.quaternion.hpp"
// #include "math.matrix.hpp"

namespace aspect
{

namespace math
{


//============================================================================
// Orientation Vectors
//
void orientation_vectors::from_quaternion(const quat& q)
{
	matrix m;
	q.to_matrix(m);

	m_locO = m * vec3(0.f, 0.f, 0.f);
	m_locX = m * vec3(1.f, 0.f, 0.f);
	m_locY = m * vec3(0.f, 1.f, 0.f);
}
void orientation_vectors::to_quaternion(quat& q) const
{
	// Displace to the origin
	vec3 locX = m_locX - m_locO;
	vec3 locY = m_locY - m_locO;

	matrix m;

	// Rotate the Y vector to the Y-axis
	vec3 locAxis;
	locY.cross(locAxis, vec3(0.f, 1.f, 0.f));
	if (!locAxis.is_zero_length())
	{
		matrix m1, m2;
		locAxis.normalize();	// Make unit vector
		double fAngle2 = locY.get_angle(vec3(0.f, 1.f, 0.f));
		m.set_identity();
		m.set_rotation(locAxis, fAngle2);
		
		// Rotate the X vector to the X-axis
		vec3 loc;
		loc = m * locX;
		double fAngle = -loc.get_angle(vec3(1.f, 0.f, 0.f));
		m1.set_identity();
		m1.set_rot_y(loc.z > 0.f ? fAngle : -fAngle);

		m2.set_identity();
		m2.set_rotation(locAxis, -fAngle2);
		m = m1 * m2;
	}
	else
	{
		// Cross-product is 0: already on the Y-axis
		if (locY.y > 0.f)
		{
			double fAngle = -locX.get_angle(vec3(1.f, 0.f, 0.f));
			m.set_rot_y(locX.z > 0.f ? fAngle : -fAngle);
		}
		else
		{
			static double fPi = (double)::acosf(-1);
			m.set_rot_x(fPi);
			vec3 loc;
			loc = m * locX;
			double fAngle = -loc.get_angle(vec3(1.f, 0.f, 0.f));
			matrix m1, m2;
			m1.set_rot_y(loc.z > 0.f ? fAngle : -fAngle);
			m2.set_rot_x(fPi);
			m = m1 * m2;
		}
	}
	q.from_matrix(m);
	q.normalize();

  
}

void orientation_vectors::apply_matrix(matrix &m)
{
	m_locO = m * m_locO;
	m_locX = m * m_locX;
	m_locY = m * m_locY;

}

} // math

} // aspect