#include "aspect.hpp"

// #include "math.hpp"
// #include "math.quaternion.hpp"
// #include "math.euler_angles.hpp"
// #include "math.matrix.hpp"

namespace aspect
{

namespace math
{


quat::quat(const euler_angles &ang) 
{ 
	ang.to_quaternion(*this); 
}


// Quaternion product this * qR.  Note: order is important!
// To combine rotations, use the product qFirst.Multiply(qSecond),
// which gives the effect of rotating by qFirst then qSecond.
void quat::multiply(const quat &qR)
{
	quat qL(*this);
    v[W] = qL.v[W]*qR.v[W] - qL.v[X]*qR.v[X] - qL.v[Y]*qR.v[Y] - qL.v[Z]*qR.v[Z];
    v[X] = qL.v[W]*qR.v[X] + qL.v[X]*qR.v[W] + qL.v[Y]*qR.v[Z] - qL.v[Z]*qR.v[Y];
    v[Y] = qL.v[W]*qR.v[Y] + qL.v[Y]*qR.v[W] + qL.v[Z]*qR.v[X] - qL.v[X]*qR.v[Z];
    v[Z] = qL.v[W]*qR.v[Z] + qL.v[Z]*qR.v[W] + qL.v[X]*qR.v[Y] - qL.v[Y]*qR.v[X];
}

// Multiply by scalar f. 
void quat::scale(axScalar f)
{
    v[W] *= f;
    v[X] *= f;
    v[Y] *= f;
    v[Z] *= f;
}

void quat::normalize()
{
	axScalar fLength = v[W]*v[W] + v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z];
	if (is_zero(fLength)) 
	{
		// Set to unit quaternion
		v[W] = 0.0;
		v[X] = 0.0;
		v[Y] = 0.0;
		v[Z] = 1.0;
	}
 	else 
	{
		axScalar f = 1.0 - fLength;
		if (is_zero(f)) scale(1.0 / (axScalar)::sqrt(fLength));
	}
}


void quat::to_matrix(matrix &m) const
{
    axScalar Nq = v[X]*v[X]+v[Y]*v[Y]+v[Z]*v[Z]+v[W]*v[W];
    axScalar s = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
    axScalar xs = v[X]*s,	  ys = v[Y]*s,	 zs = v[Z]*s;
    axScalar wx = v[W]*xs,	  wy = v[W]*ys,	 wz = v[W]*zs;
    axScalar xx = v[X]*xs,	  xy = v[X]*ys,	 xz = v[X]*zs;
    axScalar yy = v[Y]*ys,	  yz = v[Y]*zs,	 zz = v[Z]*zs;
    m[X][X] = (axScalar)(1.0 - (yy + zz)); 
	m[X][Y] = (axScalar)(xy - wz); 
	m[X][Z] = (axScalar)(xz + wy);
    m[Y][X] = (axScalar)(xy + wz); 
	m[Y][Y] = (axScalar)(1.0 - (xx + zz)); 
	m[Y][Z] = (axScalar)(yz - wx);
    m[Z][X] = (axScalar)(xz - wy); 
	m[Z][Y] = (axScalar)(yz + wx); 
	m[Z][Z] = (axScalar)(1.0 - (xx + yy));
    m[W][X]=m[W][Y]=m[W][Z]=m[X][W]=m[Y][W]=m[Z][W]=0.0; m[W][W]=1.0;
}

// Construct a unit quaternion from rotation matrix.  Assumes matrix is
// used to multiply column vector on the left: vnew = mat vold.	 Works
// correctly for right-handed coordinate system and right-handed rotations.
// Translation and perspective components ignored. 
void quat::from_matrix(const matrix &m)
{
    // This algorithm avoids near-zero divides by looking for a large component
    // - first w, then x, y, or z.  When the trace is greater than zero,
    // |w| is greater than 1/2, which is as small as a largest component can be.
    // Otherwise, the largest diagonal entry corresponds to the largest of |x|,
    // |y|, or |z|, one of which must be larger than |w|, and at least 1/2. 
    axScalar tr, s;

    tr = m[X][X] + m[Y][Y]+ m[Z][Z];
    if (tr >= 0.0) {
	    s = sqrt(tr + m[W][W]);
	    v[W] = (axScalar)(s*0.5);
	    s = axScalar(0.5) / s;
	    v[X] = (axScalar)((m[Z][Y] - m[Y][Z]) * s);
	    v[Y] = (axScalar)((m[X][Z] - m[Z][X]) * s);
	    v[Z] = (axScalar)((m[Y][X] - m[X][Y]) * s);
	} else {
	    int h = X;
	    if (m[Y][Y] > m[X][X]) h = Y;
	    if (m[Z][Z] > m[h][h]) h = Z;
	    switch (h) {
#define caseMacro(i,j,k,I,J,K) \
		    case I:\
				s = sqrt( (m[I][I] - (m[J][J]+m[K][K])) + m[W][W] );\
				v[i] = (axScalar)(s*axScalar(0.5));\
				s = axScalar(0.5) / s;\
				v[j] = (axScalar)(s * (m[I][J] + m[J][I]));\
				v[k] = (axScalar)(s * (m[K][I] + m[I][K]));\
				v[W] = (axScalar)(s * (m[K][J] - m[J][K]));\
				break
		    caseMacro(X,Y,Z,X,Y,Z);
		    caseMacro(Y,Z,X,Y,Z,X);
			caseMacro(Z,X,Y,Z,X,Y);
#undef caseMacro
	    }
	}
    if (m[W][W] != 1.0) scale((axScalar)(1.0/sqrt(m[W][W])));
}

// Spherical linear interpolation of unit quaternions with spins.  A & B -> this
void quat::slerp(axScalar _alpha, const quat &a, const quat &b, int spin)
{
	axScalar alpha = _alpha; // keep calculations in axScalar type to avoid data loss ?

	axScalar beta;			// complementary interp parameter 
	axScalar theta;			// angle between A and B 
	axScalar sin_t, cos_t;	// sine, cosine of theta 
	axScalar phi;				// theta plus spins 
	bool bflip;				// use negation of B? 

	// cosine theta = dot product of A and B 
	cos_t = a.v[X]*b.v[X] + a.v[Y]*b.v[Y] + a.v[Z]*b.v[Z] + a.v[W]*b.v[W];

	// if B is on opposite hemisphere from A, use -B instead 
 	if (cos_t < 0.0) 
	{
		cos_t = -cos_t;
		bflip = true;
	} 
	else
		bflip = false;

	// if B is (within precision limits) the same as A,
	// just linear interpolate between A and B.
	// Can't do spins, since we don't know what direction to spin.

	if (1.0 - cos_t < 1.0E-6)  // EPSILON
	{
		beta = 1.0 - alpha;
 	} else 
	{				// normal case 
 		theta = acos(cos_t);
 		phi = theta + spin * 3.14159265358979323846; // M_PI
 		sin_t = sin(theta);
 		beta = sin(theta - alpha*phi) / sin_t;
 		alpha = sin(alpha*phi) / sin_t;
 	}

	if (bflip)
		alpha = -alpha;

	// interpolate 
 	v[X] = (axScalar) (beta*a.v[X] + alpha*b.v[X]);
 	v[Y] = (axScalar) (beta*a.v[Y] + alpha*b.v[Y]);
 	v[Z] = (axScalar) (beta*a.v[Z] + alpha*b.v[Z]);
 	v[W] = (axScalar) (beta*a.v[W] + alpha*b.v[W]);
}


} // math

} // aspect