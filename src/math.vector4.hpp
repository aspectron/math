#pragma once
#ifndef __MATH_vec4_HPP__
#define __MATH_vec4_HPP__

#include "math.hpp"

namespace aspect
{

	namespace math
	{

		class MATH_API vec4 // bracket accessor to the matrix row data
		{
			public:

				union
				{
					struct
					{
						double x, y, z, w;
					};

					struct
					{
						double r, g, b, a;
					};

					struct
					{
						double left, top, width, height;
					};

					double m[4];

				};

			public:

				vec4()
				{

				}

				vec4(double _x, double _y, double _z, double _w)
					: x(_x), y(_y), z(_z), w(_w)
				{

				}

				// accessorts and casts...

				double& operator [] (int i) { return (&x)[i]; }     
				const double& operator [] (int i) const { return (&x)[i]; }  
				operator double*() { return(&x); }

				void set(double _x, double _y, double _z, double _w)
				{
					x = _x;  y = _y;  z = _z;  w = _w;
				}

				void set(const vec3 &src, double _w)
				{
					x = src.x; y = src.y; z = src.z; w = _w;
				}

				inline vec4	operator *  (double f)	const	{ return vec4(x*f, y*f, z*f, w*f); }
				inline vec4	operator /  (double f)	const	{ return vec4(x/f, y/f, z/f, w/f); }
				inline vec4 operator - (const vec4 &src) const { return vec4(x-src.x,y-src.y,z-src.z,w-src.w); }
				inline vec4 operator + (const vec4 &src) const { return vec4(x+src.x,y+src.y,z+src.z,w+src.w); }
				inline vec4& operator += (const vec4 &src) { x += src.x; y += src.y; z += src.z; w += src.w; return *this; }

				double right(void) const { return left+width; }
				double bottom(void) const { return top+height; }
		};

		inline bool operator == (vec4 &a, vec4 &b)
		{
			return (a.x-b.x) < math::precision_limit && (a.y-b.y) < math::precision_limit && (a.z-b.z) < math::precision_limit && (a.w-b.w) < math::precision_limit;
		}

	} // math

} // aspect

#endif // __MATH_vec4_HPP__