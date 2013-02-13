#include "math.hpp"
// #include "math.hpp"
// 
// 
// //#include "math.point2d.hpp"
// //#include "math.vec3.hpp"
// //#include "math.point4d.hpp"
// #include "math.matrix.hpp"
// #include "math.quaternion.hpp"
// #include "math.euler_angles.hpp"


namespace aspect
{

namespace math
{


double matrix::ms_identity[16] = 
{
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f,
};



matrix::matrix(const quat &Q)
{
	set_identity();
	set_orientation(Q);
}

matrix::matrix(const euler_angles &ang)
{
	set_identity();
	set_orientation(ang);
}



void matrix::set_orientation(const quat &quat)
{
	quat.to_matrix(*this);
}

void matrix::set_orientation(const euler_angles &ang)
{
	ang.to_matrix(*this);
}

void matrix::apply_orientation(const quat &q)
{
	matrix m;
	q.to_matrix(m);
	*this = m * (*this);
}

// These optimizations come from Graphics Gems Volume 1 
// ("Efficient Post-Concatenation of Transformation matrices", p. 476)
// and, by adding the effect in-place to an existing matrix, represent 
// a significant performance improvement over the technique of creating the 
// matrix and multiplying it with an existing one.

void matrix::apply_translation(const vec3 &t)
{
	for (int i = 0; i < 4; ++i) 
	{
		m[i][0] = m[i][0] + m[i][3]*t.x;
		m[i][1] = m[i][1] + m[i][3]*t.y;
		m[i][2] = m[i][2] + m[i][3]*t.z;
	}
}


void matrix::pre_translate( const math::vec3 &src )
{
	matrix m; m.set_identity();
	m.apply_translation(src);
	*this = m * (*this);
}


void matrix::apply_rot_x(double a)
{
	double c = ::cos(a);
	double s = ::sin(a);
	for (int i = 0; i < 4; ++i) 
	{
		double f = m[i][1];
		m[i][1] = f*c - m[i][2]*s;
		m[i][2] = f*s - m[i][2]*c;
	}
}

void matrix::apply_rot_y(double a)
{
	double c = ::cos(a);
	double s = ::sin(a);
	for (int i = 0; i < 4; ++i) 
	{
		double f = m[i][0];
		m[i][0] = f*c - m[i][2]*s;
		m[i][2] = f*s - m[i][2]*c;
	}
}

void matrix::apply_rot_z(double a)
{
	// THIS FUNCTION IS INCORRECT

/*	double c = ::cosf(a);
	double s = ::sinf(a);
	for (int i = 0; i < 4; ++i) 
	{
		double f = m[i][0];
		m[i][0] = f*c - m[i][1]*s;
		m[i][1] = f*s - m[i][1]*c;
	}
*/
}

void matrix::apply_scale(const vec3& scale)
{
	for (int i = 0; i < 4; ++i) 
	{
		m[i][0] = m[i][0]*scale.x;
		m[i][1] = m[i][1]*scale.y;
		m[i][2] = m[i][2]*scale.z;
	}
}

//============================================================================
// GAPI_InvertAffineMatrix
//
// Computes the inverse of a 3D affine matrix; i.e. a matrix with a dimen-
// sionality of 4 where the right column has the entries (0, 0, 0, 1).
//
// This procedure treats the 4 by 4 matrix as a block matrix and
// calculates the inverse of one submatrix for a significant perform-
// ance improvement over a general procedure that can invert any non-
// singular matrix:
//          --        --          --          --
//          |          | -1       |    -1      |
//          | A      0 |          |   A      0 |
//    -1    |          |          |            |
//   M   =  |          |     =    |     -1     |
//          | C      1 |          | -C A     1 |
//          |          |          |            |
//          --        --          --          --
//
//  where     M is a 4 by 4 matrix,
//            A is the 3 by 3 upper left submatrix of M,
//            C is the 1 by 3 lower left submatrix of M.
//
// Input:
//   in   - 3D affine matrix
//
// Output:
//   out  - inverse of 3D affine matrix
//
// Returned value:
//   TRUE   if input matrix is nonsingular
//   FALSE  otherwise
//
//
bool matrix::invert_affine(const matrix &in)
{
	matrix &out = *this;

    double    det_1;
    double    pos, neg, temp;

    //
    // Calculate the determinant of submatrix A and determine if the
    // the matrix is singular as limited by the double precision
    // doubleing-point data representation.
    //
    pos = neg = 0.0;
    temp =  in[0][0] * in[1][1] * in[2][2];
    if (temp >= 0.0) pos += temp; else neg += temp;
    temp =  in[0][1] * in[1][2] * in[2][0];
    if (temp >= 0.0) pos += temp; else neg += temp;
    temp =  in[0][2] * in[1][0] * in[2][1];
    if (temp >= 0.0) pos += temp; else neg += temp;
    temp = -in[0][2] * in[1][1] * in[2][0];
    if (temp >= 0.0) pos += temp; else neg += temp;
    temp = -in[0][1] * in[1][0] * in[2][2];
    if (temp >= 0.0) pos += temp; else neg += temp;
    temp = -in[0][0] * in[1][2] * in[2][1];
    if (temp >= 0.0) pos += temp; else neg += temp;
    det_1 = pos + neg;

    // Is the submatrix A singular? 
	double f = det_1 / (pos - neg);
    if (is_zero(det_1) || is_zero(f)) {
        // Matrix M has no inverse 
        return false;
    }

    // Calculate inverse(A) = adj(A) / det(A) 
    det_1 = (double)1.0 / det_1;
    out[0][0] =   ( in[1][1] * in[2][2] - in[1][2] * in[2][1] ) * det_1;
    out[1][0] = - ( in[1][0] * in[2][2] - in[1][2] * in[2][0] ) * det_1;
    out[2][0] =   ( in[1][0] * in[2][1] - in[1][1] * in[2][0] ) * det_1;
    out[0][1] = - ( in[0][1] * in[2][2] - in[0][2] * in[2][1] ) * det_1;
    out[1][1] =   ( in[0][0] * in[2][2] - in[0][2] * in[2][0] ) * det_1;
    out[2][1] = - ( in[0][0] * in[2][1] - in[0][1] * in[2][0] ) * det_1;
    out[0][2] =   ( in[0][1] * in[1][2] - in[0][2] * in[1][1] ) * det_1;
    out[1][2] = - ( in[0][0] * in[1][2] - in[0][2] * in[1][0] ) * det_1;
    out[2][2] =   ( in[0][0] * in[1][1] - in[0][1] * in[1][0] ) * det_1;

    // Calculate -C * inverse(A) 
    out[3][0] = - ( in[3][0] * out[0][0] + in[3][1] * out[1][0] + in[3][2] * out[2][0] );
    out[3][1] = - ( in[3][0] * out[0][1] + in[3][1] * out[1][1] + in[3][2] * out[2][1] );
    out[3][2] = - ( in[3][0] * out[0][2] + in[3][1] * out[1][2] + in[3][2] * out[2][2] );

    // Fill in last column 
    out[0][3] = out[1][3] = out[2][3] = 0.0;
    out[3][3] = 1.0;

//	in = out;

    return true;
}

void matrix::get_orientation( math::quat &q )
{
	q.from_matrix(*this);
}

void matrix::get_translation( math::vec3 &loc )
{
	loc = math::vec3(m_41,m_42,m_43);
}

void matrix::look_at( const vec3& eye, const vec3& target, const vec3& up )
{

	// Please note that look_at matrix is already inverted, intended to be combined with viewport.

	//m_eye = eye;
	set_identity();

	vec3 z = eye - target;
	z.normalize();

//	z = -z;

	// vec3 view_dir = -z;

	vec3 x; x.cross(up, z);
	x.normalize();

	vec3 y; y.cross(z, x);
	y.normalize();
//	x.normalize();

	m[0][0] = x.x;
	m[1][0] = x.y;
	m[2][0] = x.z;
	//	m[3][0] = -x.dot(eye); //vec::dot(m_xAxis, eye);
	m[3][0] = x.dot(-eye); //vec::dot(m_xAxis, eye);
//	m[3][0] = -eye.x; //vec::dot(m_xAxis, eye);

	m[0][1] = y.x;
	m[1][1] = y.y;
	m[2][1] = y.z;
	//m[3][1] = -eye.y; //Vector3::dot(m_yAxis, eye);
	//m[3][1] = y.dot(eye); //Vector3::dot(m_yAxis, eye);
	m[3][1] = y.dot(-eye); //Vector3::dot(m_yAxis, eye);

	m[0][2] = z.x;
	m[1][2] = z.y;
	m[2][2] = z.z;    
	//m[3][2] = -eye.z; //Vector3::dot(m_zAxis, eye);
	m[3][2] = z.dot(-eye); //Vector3::dot(m_zAxis, eye);
	//m[3][2] = z.dot(eye); //Vector3::dot(m_zAxis, eye);

	/*
	vec3 x,y,z;

	z = (eye - target).normalize();
	// z.subVectors( eye, target ).normalize();

	if ( z.length() == 0.0 ) {
	z.z = 1;
	}

	// x.crossVectors( up, z ).normalize();
	x.cross(up,z).normalize();

	if ( x.length() == 0.0 ) {

	z.x += 0.0001;
	//x.crossVectors( up, z ).normalize();
	x.cross(up,z).normalize();
	}

	//	y.crossVectors( z, x );
	y.cross(z,x);

	set_identity();
	v[0] = x.x; v[4] = y.x; v[8] = z.x;
	v[1] = x.y; v[5] = y.y; v[9] = z.y;
	v[2] = x.z; v[6] = y.z; v[10] = z.z;
	*/
}


} // math

} // aspect
