#ifndef _MATH_HPP_
#define _MATH_HPP_

#include <math.h>

//#define bASPECT_IS_ZERO(f)   (fabs(f)<GAPI_PRECISION_LIMIT)

// ---
/*
#define ASPECT_PI				(3.1415926535f)
#define ASPECT_TWOPI			(6.283185307f)
#define ASPECT_HALFPI			(1.570796326794895f)
#define ASPECT_DEG_TO_RAD		(ASPECT_PI/180.0f)
#define ASPECT_RAD_TO_DEG		(180.0f/ASPECT_PI)
#define ASPECT_DegToRad(fDeg)	(((float)fDeg)*ASPECT_DEG_TO_RAD)
#define ASPECT_RadToDeg(fRad)	(((float)fRad)*ASPECT_RAD_TO_DEG)
*/

typedef float axScalar;

namespace aspect
{
	namespace math
	{
		const float precision_limit = 0.000001f;
		inline bool is_zero(axScalar f)   { return ( fabs(f) < precision_limit ); }

		// ---

		const double _e        = 2.71828182845904523536;
		const double _log2e    = 1.44269504088896340736;
		const double _log10e   = 0.434294481903251827651;
		const double _ln2      = 0.693147180559945309417;
		const double _ln10     = 2.30258509299404568402;
		const double _twopi	   = 6.283185307;
		const double _pi       = 3.14159265358979323846;
		const double _pi_2     = 1.57079632679489661923;
		const double _halfpi   = 1.57079632679489661923;
		const double _pi_4     = 0.785398163397448309616;
		const double _1_pi     = 0.318309886183790671538;
		const double _2_pi     = 0.636619772367581343076;
		const double _2_sqrtpi = 1.12837916709551257390;
		const double _sqrt2    = 1.41421356237309504880;
		const double _sqrt1_2  = 0.707106781186547524401;

		const double _deg_to_rad = _pi/180.0;
		const double _rad_to_deg = 180.0/_pi;
		inline double deg_to_rad(double deg)	{ return deg * _deg_to_rad; }
		inline double rad_to_deg(double rad) { return rad * _rad_to_deg; }

		template <class T> inline void swap(T& a, T& b) { T c = a; a = b; b = c; }

		inline bool feq(double f1, double f2) { return (fabs(f1-f2) < precision_limit); }

		template<class T> inline T clamp(T n, T _min, T _max) { if(n < _min) return _min; else if(n > _max) return _max; return n; }

		inline float random_neg_1_to_1(void) { return (float)rand() / (float)(RAND_MAX / 2) - 1.0f; }
		inline float random_0_to_1(void) { return (float)rand() / (float)RAND_MAX; }

	} // math
} // aspect

// #include "math.vector2.hpp"
// #include "math.vector3.hpp"
// #include "math.vector4.hpp"
// #include "math.matrix.hpp"
// #include "math.quaternion.hpp"
// #include "math.euler_angles.hpp"
// #include "math.orientation_vectors.hpp"

#endif // _MATH_HPP_