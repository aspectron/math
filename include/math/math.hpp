#ifndef MATH_HPP_INCLUDED
#define MATH_HPP_INCLUDED

#include <math.h>
#include <stdlib.h>

#include "jsx/core.hpp"
#include "jsx/v8_core.hpp"

#if OS(WINDOWS)
//	#pragma warning ( disable : 4251 )
#if defined(MATH_EXPORTS)
#define MATH_API __declspec(dllexport)
#else
#define MATH_API __declspec(dllimport)
#endif
#elif __GNUC__ >= 4
# define MATH_API __attribute__((visibility("default")))
#else
#define MATH_API // nothing, symbols in a shared library are exported by default
#endif

namespace aspect { namespace math {
const float precision_limit = 0.000001f;
inline bool is_zero(double f)   { return ( fabs(f) < precision_limit ); }

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

}} // aspect::math

#include "math/vector2.hpp"
#include "math/vector3.hpp"
#include "math/vector4.hpp"
#include "math/matrix.hpp"
#include "math/quaternion.hpp"
#include "math/euler_angles.hpp"
#include "math/orientation_vectors.hpp"

#endif // MATH_HPP_INCLUDED
