#ifndef MATH_MATRIX_HPP_INCLUDED
#define MATH_MATRIX_HPP_INCLUDED

#include "math/math.hpp"
#include "math/vector2.hpp"
#include "math/vector3.hpp"
#include "math/vector4.hpp"
#include "math/quaternion.hpp"

namespace aspect { namespace math {

class quat;
class euler_angles;

struct MATH_API matrix_row
{
	double v[4];
	double& operator[] (int i) { return v[i]; }
	double operator[] (int i) const { return v[i]; }
};

class MATH_API matrix
{
public:
	union
	{
		struct
		{
			double m_11, m_12, m_13, m_14;
			double m_21, m_22, m_23, m_24;
			double m_31, m_32, m_33, m_34;
			double m_41, m_42, m_43, m_44;
		};

		double m[4][4];

		matrix_row m_row[4];

		double v[16];
	};

public:

	static double ms_identity[16];

	// constructors...

	matrix() {}

	matrix(const double f11, const double f12, const double f13, const double f14,
				const double f21, const double f22, const double f23, const double f24,
				const double f31, const double f32, const double f33, const double f34,
				const double f41, const double f42, const double f43, const double f44)
		: m_11(f11), m_12(f12), m_13(f13), m_14(f14),
			m_21(f21), m_22(f22), m_23(f23), m_24(f24),
			m_31(f31), m_32(f32), m_33(f33), m_34(f34),
			m_41(f41), m_42(f42), m_43(f43), m_44(f44)
	{
	}

/*		matrix( const __m128 &ms1, const __m128 &ms2, const __m128 &ms3, const __m128 &ms4)
		: m1(ms1), m2(ms2), m3(ms3), m4(ms4)
	{
	}
*/
	explicit matrix(const quat &Q)
	{
		set_identity();
		set_orientation(Q);
	}

	explicit matrix(const euler_angles &ang)
	{
		set_identity();
		set_orientation(ang);
	}

	// initializers...
	void set_identity()
	{
		memcpy(&(m[0][0]),ms_identity,sizeof(ms_identity));
//			m_12 = m_13 = m_14 = m_21 = m_23 = m_24 = 
//			m_31 = m_32 = m_34 = m_41 = m_42 = m_43 = 0.0f;
//			m_11 = m_22 = m_33 = m_44 = 1.0f;
	}

	void zero() { ::memset(this,0,sizeof(matrix)); }
 
	// accessors and casts...
	double & operator () (int x, int y) { return m[x][y]; }
	double * operator [] (int i) { return m[i]; }
	const vec4& operator [] (int i) const { return((vec4&)(*m[i])); }

	// operators...
	matrix & operator *= (const matrix &M);//	{ *this = (*this) * (M); return *this; }

	double& at(int i) { return v[i]; }

	// ---
	void set_translation(const vec3 &Dist);
	void set_scale(const vec3 &Scale);
	void set_rot_x(double fAngle);
	void set_rot_y(double fAngle);
	void set_rot_z(double fAngle);
	void set_rotation(const vec3& axis, double angle);
	void set_orientation(const quat &quat);
	void set_orientation(const euler_angles &ang);
		
	double det();
	void adjoint(const matrix& in);

	bool invert(const matrix &src);
	bool invert_affine(const matrix &src);

	void look_at(const vec3& eye, const vec3& target, const vec3& up);

	void apply_translation(const vec3 &t);
	void apply_rot_x(double a);
	void apply_rot_y(double a);
	void apply_rot_z(double a);
	void apply_scale(const vec3& scale);
	void apply_orientation(const quat &quat);

	void pre_translate(const math::vec3 &src);

//	void get_location(math::vec3 &loc) { loc = math::vec3(0.0f,0.0f,0.0f); return *this * loc; }

	void get_scale(math::vec3 &scale)
	{
		scale.x = math::vec3(m_11,m_21,m_31).length();
		scale.y = math::vec3(m_12,m_22,m_32).length();
		scale.z = math::vec3(m_13,m_23,m_33).length();
	}

	void get_translation(math::vec3 &loc);
	void get_orientation(math::quat &q);

	math::vec3 scale() const
	{
		return math::vec3(
			math::vec3(m_11, m_21, m_31).length(),
			math::vec3(m_12, m_22, m_32).length(),
			math::vec3(m_13, m_23, m_33).length());
	}

	math::vec3 translation() const
	{
		return math::vec3(m_41, m_42, m_43);
	}

	math::quat orientation() const
	{
		return math::quat(*this);
	}
};

inline
//__forceinline
matrix operator * (const matrix &A, const matrix &B)
{
/*
#ifndef _DEBUG
	if(! (((size_t)&A | (size_t)&B) & 15) )
	{
		return matrix(
			_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(A.m1, A.m1, _MM_SHUFFLE(0,0,0,0)), B.m1), _mm_mul_ps(_mm_shuffle_ps(A.m1, A.m1, _MM_SHUFFLE(1,1,1,1)), B.m2)), _mm_mul_ps(_mm_shuffle_ps(A.m1, A.m1, _MM_SHUFFLE(2,2,2,2)), B.m3)), _mm_mul_ps(_mm_shuffle_ps(A.m1, A.m1, _MM_SHUFFLE(3,3,3,3)), B.m4)),
			_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(A.m2, A.m2, _MM_SHUFFLE(0,0,0,0)), B.m1), _mm_mul_ps(_mm_shuffle_ps(A.m2, A.m2, _MM_SHUFFLE(1,1,1,1)), B.m2)), _mm_mul_ps(_mm_shuffle_ps(A.m2, A.m2, _MM_SHUFFLE(2,2,2,2)), B.m3)), _mm_mul_ps(_mm_shuffle_ps(A.m2, A.m2, _MM_SHUFFLE(3,3,3,3)), B.m4)),
			_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(A.m3, A.m3, _MM_SHUFFLE(0,0,0,0)), B.m1), _mm_mul_ps(_mm_shuffle_ps(A.m3, A.m3, _MM_SHUFFLE(1,1,1,1)), B.m2)), _mm_mul_ps(_mm_shuffle_ps(A.m3, A.m3, _MM_SHUFFLE(2,2,2,2)), B.m3)), _mm_mul_ps(_mm_shuffle_ps(A.m3, A.m3, _MM_SHUFFLE(3,3,3,3)), B.m4)),
			_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(A.m4, A.m4, _MM_SHUFFLE(0,0,0,0)), B.m1), _mm_mul_ps(_mm_shuffle_ps(A.m4, A.m4, _MM_SHUFFLE(1,1,1,1)), B.m2)), _mm_mul_ps(_mm_shuffle_ps(A.m4, A.m4, _MM_SHUFFLE(2,2,2,2)), B.m3)), _mm_mul_ps(_mm_shuffle_ps(A.m4, A.m4, _MM_SHUFFLE(3,3,3,3)), B.m4))
		);
	}
#endif
*/
	return matrix(
		A.m[0][0]*B.m[0][0] + A.m[0][1]*B.m[1][0] + A.m[0][2]*B.m[2][0] + A.m[0][3]*B.m[3][0],
		A.m[0][0]*B.m[0][1] + A.m[0][1]*B.m[1][1] + A.m[0][2]*B.m[2][1] + A.m[0][3]*B.m[3][1],
		A.m[0][0]*B.m[0][2] + A.m[0][1]*B.m[1][2] + A.m[0][2]*B.m[2][2] + A.m[0][3]*B.m[3][2],
		A.m[0][0]*B.m[0][3] + A.m[0][1]*B.m[1][3] + A.m[0][2]*B.m[2][3] + A.m[0][3]*B.m[3][3],

		A.m[1][0]*B.m[0][0] + A.m[1][1]*B.m[1][0] + A.m[1][2]*B.m[2][0] + A.m[1][3]*B.m[3][0],
		A.m[1][0]*B.m[0][1] + A.m[1][1]*B.m[1][1] + A.m[1][2]*B.m[2][1] + A.m[1][3]*B.m[3][1],
		A.m[1][0]*B.m[0][2] + A.m[1][1]*B.m[1][2] + A.m[1][2]*B.m[2][2] + A.m[1][3]*B.m[3][2],
		A.m[1][0]*B.m[0][3] + A.m[1][1]*B.m[1][3] + A.m[1][2]*B.m[2][3] + A.m[1][3]*B.m[3][3],

		A.m[2][0]*B.m[0][0] + A.m[2][1]*B.m[1][0] + A.m[2][2]*B.m[2][0] + A.m[2][3]*B.m[3][0],
		A.m[2][0]*B.m[0][1] + A.m[2][1]*B.m[1][1] + A.m[2][2]*B.m[2][1] + A.m[2][3]*B.m[3][1],
		A.m[2][0]*B.m[0][2] + A.m[2][1]*B.m[1][2] + A.m[2][2]*B.m[2][2] + A.m[2][3]*B.m[3][2],
		A.m[2][0]*B.m[0][3] + A.m[2][1]*B.m[1][3] + A.m[2][2]*B.m[2][3] + A.m[2][3]*B.m[3][3],

		A.m[3][0]*B.m[0][0] + A.m[3][1]*B.m[1][0] + A.m[3][2]*B.m[2][0] + A.m[3][3]*B.m[3][0],
		A.m[3][0]*B.m[0][1] + A.m[3][1]*B.m[1][1] + A.m[3][2]*B.m[2][1] + A.m[3][3]*B.m[3][1],
		A.m[3][0]*B.m[0][2] + A.m[3][1]*B.m[1][2] + A.m[3][2]*B.m[2][2] + A.m[3][3]*B.m[3][2],
		A.m[3][0]*B.m[0][3] + A.m[3][1]*B.m[1][3] + A.m[3][2]*B.m[2][3] + A.m[3][3]*B.m[3][3]
	);
}

inline
//__forceinline 
matrix &matrix::operator *= (const matrix &Src)
{
/*
#ifndef _DEBUG
	if(! (((size_t)this | (size_t)&Src) & 15) )
	{
		m1 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(0,0,0,0)), Src.m1), _mm_mul_ps(_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(1,1,1,1)), Src.m2)), _mm_mul_ps(_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(2,2,2,2)), Src.m3)), _mm_mul_ps(_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(3,3,3,3)), Src.m4));
		m2 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(m2, m2, _MM_SHUFFLE(0,0,0,0)), Src.m1), _mm_mul_ps(_mm_shuffle_ps(m2, m2, _MM_SHUFFLE(1,1,1,1)), Src.m2)), _mm_mul_ps(_mm_shuffle_ps(m2, m2, _MM_SHUFFLE(2,2,2,2)), Src.m3)), _mm_mul_ps(_mm_shuffle_ps(m2, m2, _MM_SHUFFLE(3,3,3,3)), Src.m4));
		m3 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(m3, m3, _MM_SHUFFLE(0,0,0,0)), Src.m1), _mm_mul_ps(_mm_shuffle_ps(m3, m3, _MM_SHUFFLE(1,1,1,1)), Src.m2)), _mm_mul_ps(_mm_shuffle_ps(m3, m3, _MM_SHUFFLE(2,2,2,2)), Src.m3)), _mm_mul_ps(_mm_shuffle_ps(m3, m3, _MM_SHUFFLE(3,3,3,3)), Src.m4));
		m4 = _mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_shuffle_ps(m4, m4, _MM_SHUFFLE(0,0,0,0)), Src.m1), _mm_mul_ps(_mm_shuffle_ps(m4, m4, _MM_SHUFFLE(1,1,1,1)), Src.m2)), _mm_mul_ps(_mm_shuffle_ps(m4, m4, _MM_SHUFFLE(2,2,2,2)), Src.m3)), _mm_mul_ps(_mm_shuffle_ps(m4, m4, _MM_SHUFFLE(3,3,3,3)), Src.m4));
		return *this;
	}
#endif
*/
	*this = (*this) * (Src); 
	return *this;
}

inline void matrix::set_translation(const vec3 &Dist)
{
	// WARNING:  This function assumes identity matrix.
	m[3][0] = Dist.x;
	m[3][1] = Dist.y;
	m[3][2] = Dist.z;
}

inline void matrix::set_scale(const vec3 &Scale)
{
	// WARNING:  This function assumes identity matrix.
	m[0][0] = Scale.x;
	m[1][1] = Scale.y;
	m[2][2] = Scale.z;
}

inline void matrix::set_rot_x(double fAngle)
{
	// WARNING:  This function assumes identity matrix.
	m[2][1] = -(m[1][2] = sin(fAngle));
	m[2][2] = m[1][1] = cos(fAngle);
}

inline void matrix::set_rot_y(double fAngle)
{
	// WARNING:  This function assumes identity matrix.
	m[0][0] = m[2][2] = cos(fAngle);
	m[0][2] = -(m[2][0] = sin(fAngle));
}

inline void matrix::set_rot_z(double fAngle)
{
	// WARNING:  This function assumes identity matrix.

	m[1][0] = -(m[0][1] = sin(fAngle));
	m[0][0] = m[1][1] = cos(fAngle);
}

inline vec2& vec2::operator *= (const matrix &m)
{
	vec3 p(x, y, 0.0f);
	p *= m;
	x = p.x;
	y = p.y;

	return *this;
}

inline vec3&	vec3::operator *= (const matrix &A)//			{   *this = m * (*this); return *this; }
{
	vec3 l;
	vec3 &B = *this;

	double fX, fY, fZ, fW;

	fX = A.m[0][0] * B.x + A.m[1][0] * B.y + A.m[2][0] * B.z + A.m[3][0];
	fY = A.m[0][1] * B.x + A.m[1][1] * B.y + A.m[2][1] * B.z + A.m[3][1];
	fZ = A.m[0][2] * B.x + A.m[1][2] * B.y + A.m[2][2] * B.z + A.m[3][2];
	fW = A.m[0][3] * B.x + A.m[1][3] * B.y + A.m[2][3] * B.z + A.m[3][3];

	if(fW == 0.0f)
		l = vec3(0,0,0);
	else if (fW == 1.0f)
		l = vec3(fX, fY, fZ);
	else
		l = vec3(fX/fW, fY/fW, fZ/fW);

	*this = l;
	return *this;
}

inline vec3 operator * (const matrix &A, const vec3 &B)
{
	vec3 l;
	double fX, fY, fZ, fW;

	fX = A.m[0][0] * B.x + A.m[1][0] * B.y + A.m[2][0] * B.z + A.m[3][0];
	fY = A.m[0][1] * B.x + A.m[1][1] * B.y + A.m[2][1] * B.z + A.m[3][1];
	fZ = A.m[0][2] * B.x + A.m[1][2] * B.y + A.m[2][2] * B.z + A.m[3][2];
	fW = A.m[0][3] * B.x + A.m[1][3] * B.y + A.m[2][3] * B.z + A.m[3][3];

	if(fW == 0.0f)
		l = vec3(0,0,0);
	else if (fW == 1.0f)
		l = vec3(fX, fY, fZ);
	else
		l = vec3(fX/fW, fY/fW, fZ/fW);

	return l;
}

inline vec3 operator * (const vec3 &B, const matrix &A)
{
	vec3 l;
	double fX, fY, fZ, fW;

	fX = A.m[0][0] * B.x + A.m[1][0] * B.y + A.m[2][0] * B.z + A.m[3][0];
	fY = A.m[0][1] * B.x + A.m[1][1] * B.y + A.m[2][1] * B.z + A.m[3][1];
	fZ = A.m[0][2] * B.x + A.m[1][2] * B.y + A.m[2][2] * B.z + A.m[3][2];
	fW = A.m[0][3] * B.x + A.m[1][3] * B.y + A.m[2][3] * B.z + A.m[3][3];

	if(fW == 0.0f)
		l = vec3(0,0,0);
	else if (fW == 1.0f)
		l = vec3(fX, fY, fZ);
	else
		l = vec3(fX/fW, fY/fW, fZ/fW);

	return l;
}

// Rotation at origin by angle about unit-vector axis.
inline void matrix::set_rotation(const vec3& axis, double angle)
{
	double c = (double)::cos(angle);
	double s = (double)::sin(angle);
	double t = 1.f - c;
	double x = axis.x;
	double y = axis.y;
	double z = axis.z;

	m[0][0] = t*x*x + c;	m[0][1] = t*x*y + s*z;	m[0][2] = t*x*z - s*y;	m[0][3] = 0.f;
	m[1][0] = t*x*y - s*z;	m[1][1] = t*y*y + c;	m[1][2] = t*y*z + s*x;	m[1][3] = 0.f;
	m[2][0] = t*x*z + s*y;	m[2][1] = t*y*z - s*x;	m[2][2] = t*z*z + c;	m[2][3] = 0.f;
	m[3][0] = 0.f;			m[3][1] = 0.f;			m[3][2] = 0.f;			m[3][3] = 1.f;
}



//============================================================================
// Matrix Inversion
// by Richard Carling
// from "Graphics Gems", Academic Press, 1990
//

// Calculate the determinent of a 2x2 matrix.
static inline double det2x2( double a, double b, double c, double d)
{
	return a * d - b * c;
}

// calculate the determinent of a 3x3 matrix in the form
//
//     | a1,  b1,  c1 |
//     | a2,  b2,  c2 |
//     | a3,  b3,  c3 |
//
static inline double det3x3( 
	double a1, double a2, double a3, 
	double b1, double b2, double b3, 
	double c1, double c2, double c3 
) {
	return (a1 * det2x2( b2, b3, c2, c3 )
			- b1 * det2x2( a2, a3, c2, c3 )
			+ c1 * det2x2( a2, a3, b2, b3 ));
}

// Calculate the determinent of a 4x4 matrix.
inline double matrix::det(void)
{
	double ans;
	double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

	// assign to individual variable names to aid selecting 
	// correct elements 

	a1 = m[0][0]; b1 = m[0][1]; 
	c1 = m[0][2]; d1 = m[0][3];

	a2 = m[1][0]; b2 = m[1][1]; 
	c2 = m[1][2]; d2 = m[1][3];

	a3 = m[2][0]; b3 = m[2][1]; 
	c3 = m[2][2]; d3 = m[2][3];

	a4 = m[3][0]; b4 = m[3][1]; 
	c4 = m[3][2]; d4 = m[3][3];

	ans = a1 * det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4)
		- b1 * det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4)
		+ c1 * det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4)
		- d1 * det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

	return ans;
}

// 
//     calculate the adjoint of a 4x4 matrix
//
//      Let  a   denote the minor determinant of matrix A obtained by
//           ij
//
//      deleting the ith row and jth column from A.
//
//                    i+j
//     Let  b   = (-1)    a
//          ij            ji
//
//    The matrix B = (b  ) is the adjoint of A
//                     ij
//
inline void matrix::adjoint(const matrix &in)
{
//	matrix &out = *this;

	double a1, a2, a3, a4, b1, b2, b3, b4;
	double c1, c2, c3, c4, d1, d2, d3, d4;

	// assign to individual variable names to aid  
	// selecting correct values  

	a1 = in[0][0]; b1 = in[0][1]; 
	c1 = in[0][2]; d1 = in[0][3];

	a2 = in[1][0]; b2 = in[1][1]; 
	c2 = in[1][2]; d2 = in[1][3];

	a3 = in[2][0]; b3 = in[2][1];
	c3 = in[2][2]; d3 = in[2][3];

	a4 = in[3][0]; b4 = in[3][1]; 
	c4 = in[3][2]; d4 = in[3][3];

	// row column labeling reversed since we transpose rows & columns 

	m[0][0]  =   det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
	m[1][0]  = - det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
	m[2][0]  =   det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
	m[3][0]  = - det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

	m[0][1]  = - det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
	m[1][1]  =   det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
	m[2][1]  = - det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
	m[3][1]  =   det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);

	m[0][2]  =   det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
	m[1][2]  = - det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
	m[2][2]  =   det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
	m[3][2]  = - det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);

	m[0][3]  = - det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
	m[1][3]  =   det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
	m[2][3]  = - det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
	m[3][3]  =   det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);
}


// 
//    calculate the inverse of a 4x4 matrix
//
//     -1     
//     A  = ___1__ adjoint A
//         det A
//

inline
//__forceinline 
bool matrix::invert( const matrix &in )
{
//	matrix &out = *this;

//    int i, j;

	// calculate the adjoint matrix 

	//out.
	adjoint(in);

//    GAPI_AdjointMatrix( out, in );

	//  calculate the 4x4 determinent
	//  if the determinent is zero, 
	//  then the inverse matrix is not unique.

	double fDet = det();  //GAPI_DeterminantMatrix( out );
	//    double fDet = out.fDeterminant();  //GAPI_DeterminantMatrix( out );

	if (is_zero(fDet)) {
		return false;
	}

	// scale the adjoint matrix to get the inverse 

//    for (i=0; i<4; i++)
//        for(j=0; j<4; j++)
//		    out[i][j] /= fDet;

	m[0][0] /= fDet;
	m[0][1] /= fDet;
	m[0][2] /= fDet;
	m[0][3] /= fDet;

	m[1][0] /= fDet;
	m[1][1] /= fDet;
	m[1][2] /= fDet;
	m[1][3] /= fDet;

	m[2][0] /= fDet;
	m[2][1] /= fDet;
	m[2][2] /= fDet;
	m[2][3] /= fDet;

	m[3][0] /= fDet;
	m[3][1] /= fDet;
	m[3][2] /= fDet;
	m[3][3] /= fDet;

	return true;
}


}} // aspect::math

#endif // MATH_MATRIX_HPP_INCLUDED
