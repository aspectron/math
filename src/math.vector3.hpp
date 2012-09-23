#pragma once
#ifndef __MATH_vec3_HPP__
#define __MATH_vec3_HPP__

#include "math.hpp"
#include "math.vector2.hpp"

namespace aspect
{

	namespace math
	{

		class MATH_API vec3
		{

			public:				

				union
				{
					struct { double x,  y,  z; };
					struct { double r,  g,  b; };
					struct { double Y, Cb, Cr; };
					double v[3];
				};

				// constructors...

				vec3	() { x = y = z = 0.0; }
				vec3	(double _x, double _y, double _z = 0.0) {x = _x; y = _y; z = _z;}
				vec3	(const vec2& src, double _z) { x = src.x; y = src.y; z = _z; }
//				vec3	(int _x, int _y, int _z = 0) {x = (double)_x; y = (double)_y; z = (double)_z;}

				// initializers...

				void zero(void) { x = y = z = 0.0; }

				// accessors and casts...

				double &operator [] (int i) { return v[i]; }
				operator double*() { return(v); }
				operator const double*() { return(v); }

				// operators...

				vec3&	operator =  (const vec2 &pt) { x = pt.x; y = pt.y; z = 0.0f; return *this; }

				vec3&	operator *= (const vec3 &w);			
				vec3&	operator *= (const matrix &w);			
				vec3&	operator += (const vec3 &w);			
				vec3&	operator -= (const vec3 &w);			
				vec3	operator -  (const vec3 &w)	const;	
				vec3	operator +  (const vec3 &w)	const;	
				vec3	operator /  (const vec3 &w)	const;	
				vec3	operator *  (double f)				const;	
				vec3	operator /  (double f)				const;	

				vec3&	operator *= (const double f);
				vec3&	operator += (const double f);
				vec3&	operator -= (const double f);

				vec3 operator -() const; 
				vec3 operator +() const;

				bool       operator == (const vec3 &w) const	{return (w-*this).is_zero_length();}
				bool       operator != (const vec3 &w) const	{return !(w-*this).is_zero_length();}

				// misc functions...

				bool		is_zero_length	() const	{ return is_zero(x) && is_zero(y) && is_zero(z); }
				double		length		() const;							// Compute the length of a vector
				double		distance		(const vec3& loc) const;		// Compute the distance between two points
				double		distance		(const vec2& loc) const;		// Compute the distance between two points
				double		sqr_distance	(const vec3& loc) const;		// Compute SQR distance between two points
				double		get_angle		(const vec3& loc) const;		// Get the angle between 2 vectors
				double		max_value		() const	{ return x > y ? (x > z ? x : z) : (y > z ? y : z); }

				// Get a point along line AB, at a distance from A expressed in relation to
				// the distance from A to B.  Use this to project(extrapolate) or interpolate a line.
				void		get_point_on_line	(vec3& locA, vec3& locB, const double f);

				void		get_nearest_point_on_line(const vec3 &pt, const vec3 &ptLine1, const vec3 &ptLine2);


				double		dot				(const vec3& loc) const;
				void		cross			(vec3& locOut, const vec3& loc) const;
				void		normalize		();
				void		clamp			(double _min, double _max);

				void		from_string(const char *sz);
				void		to_string(char *sz, int iMaxLen);

				// color operations

				/*
		mat3           rgb_to_ycbcr_mat = mat3( 0.257, 0.504, 0.098, -0.148, -0.291, 0.439, 0.439, -0.368, -0.071 );
		vec3           rgb_to_ycbcr_vec = vec3( 16.0/255.0, 128.0/255.0, 128.0/255.0 );

				*/
/*
				const vec3 RGB_to_YCbCr(void)
				{
					static const vec3 rgb_to_ycbcr_vec(16.0/255.0, 128.0/255.0, 128.0/255.0 );
					vec3 rgb = (*this);// * 255.0;
					vec3 pt;
					pt.Y = (0.257 * rgb.r + 0.504 * rgb.g + 0.098 * rgb.b);
					pt.Cb = ((-0.148) * rgb.r - 0.291 * rgb.g + 0.439 * rgb.b);
					pt.Cr = (0.439 * rgb.r - 0.368 * rgb.g - 0.071 * rgb.b);

		//			pt.Y = (0.299 * rgb.r + 0.587 * rgb.g + 0.114 * rgb.b) / 255.0;
		//			pt.Cb = (((-0.1687) * rgb.r - 0.3313 * rgb.g + 0.5 * rgb.b)+128.0) / 255.0;
		//			pt.Cr = ((0.5 * rgb.r - 0.4187 * rgb.g - 0.0813 * rgb.b)+128.0) / 255.0;
					return pt + rgb_to_ycbcr_vec;
				}

				const vec3 YCbCr_to_RGB(void)
				{
					static const vec3 ycbcr_to_rgb_vec(-0.87416470588235, +0.53132549019608, -1.08599215686274);
					vec3 n = (*this);
					vec3 pt;
					pt.r = (1.164 * n.Y + 1.596 * n.Cr);
					pt.g = (1.164 * n.Y - 0.391 * n.Cb - 0.813 * n.Cr);
					pt.b = (1.164 * n.Y + 2.018 * n.Cb);

					return pt + ycbcr_to_rgb_vec;
				}
*/
		};

		inline vec3&	vec3::operator *= (const vec3 &w)			{ x *= w.x; y *= w.y; z *= w.z; return *this; }
		inline vec3&	vec3::operator += (const vec3 &w)			{ x += w.x; y += w.y; z += w.z; return *this; }
		inline vec3&	vec3::operator -= (const vec3 &w)			{ x -= w.x; y -= w.y; z -= w.z; return *this; }
		inline vec3	vec3::operator -  (const vec3 &w)	const	{ return vec3(x-w.x, y-w.y, z-w.z); }
		inline vec3	vec3::operator +  (const vec3 &w)	const	{ return vec3(x+w.x, y+w.y, z+w.z); }
		inline vec3	vec3::operator /  (const vec3 &w)	const	{ return vec3(x/w.x, y/w.y, z/w.z); }
		inline vec3	vec3::operator *  (double f)				const	{ return vec3(x*f, y*f, z*f); }
		inline vec3	vec3::operator /  (double f)				const	{ return vec3(x/f, y/f, z/f); }

		inline vec3&	vec3::operator *= (const double f)			{ x *= f; y *= f; z *= f; return *this; }
		inline vec3&	vec3::operator += (const double f)			{ x += f; y += f; z += f; return *this; }
		inline vec3&	vec3::operator -= (const double f)			{ x -= f; y -= f; z -= f; return *this; }

		inline vec3	vec3::operator -() const { return (vec3(-x, -y, -z)); }
		inline vec3	vec3::operator +() const { return (*this); }

		inline vec3	operator * (const vec3 &p1, const vec3 &p2) { return vec3(p1.x*p2.x,p1.y*p2.y,p1.z*p2.z); } 
		inline vec3	operator / (const double f, const vec3 &p) { return vec3(f/p.x,f/p.y,f/p.z); } 

		inline vec3 operator ^(const vec3& A, const vec3& B)
		{
			return vec3(	A.y * B.z - A.z * B.y,
								A.z * B.x - A.x * B.z,
								A.x * B.y - A.y * B.x);
		}

		// vec3::fDistance()
		//   Compute the distance between two locs
		inline double vec3::distance(const vec3& loc) const 
		{
			return sqrt( (loc.x-x) * (loc.x-x)
						 +(loc.y-y) * (loc.y-y)
						 +(loc.z-z) * (loc.z-z));
		}

		inline double vec3::distance(const vec2& loc) const 
		{
			return sqrt( (loc.x-x) * (loc.x-x)
						 +(loc.y-y) * (loc.y-y));
		}

		inline double vec3::sqr_distance(const vec3& loc) const 
		{
			return ( (loc.x-x) * (loc.x-x)
					+(loc.y-y) * (loc.y-y)
					+(loc.z-z) * (loc.z-z));
		}

		// vec3::fLength()
		//   Compute the length of the vector formed by the loc with the origin
		inline double vec3::length() const 
		{
			return sqrt(x*x + y*y + z*z);
		}

		inline double vec3::get_angle(const vec3& loc) const
		{
			double f = length() * loc.length();
			if (is_zero(f)) return 0.f;
			f = dot(loc) / f;
			if (f < -1.0) f = -1.0;
			if (f >  1.0) f =  1.0;
			return ::acos(f);
		}

		// vec3::GetPointOnLine()
		//   Get a point along line AB, at a distance from A expressed in relation to
		//   the distance from A to B.  Use this to project or interpolate a line.
		inline void vec3::get_point_on_line(vec3& locA, vec3& locB, const double f)
		{
			x = locA.x + (f * (locB.x - locA.x));
			y = locA.y + (f * (locB.y - locA.y));
			z = locA.z + (f * (locB.z - locA.z));
		}

		inline double vec3::dot(const vec3& loc) const
		{
			return x * loc.x + y * loc.y + z * loc.z;
		}

		inline void vec3::cross(vec3& locOut, const vec3& loc) const
		{
			locOut.x = y * loc.z - z * loc.y;
			locOut.y = z * loc.x - x * loc.z;
			locOut.z = x * loc.y - y * loc.x;
		}

		inline void vec3::normalize()	// Convert to unit-vector
		{
			double f = distance(vec3(0.f, 0.f, 0.f));
			if (f > 0.0f) 
			{
				x /= f;
				y /= f;
				z /= f;
			}
			else
			{
				x = 0.0f; 
				y = 0.0f; 
				z = -1.0f; 
			}
		}

		inline void vec3::get_nearest_point_on_line(const vec3 &pt, const vec3 &ptA, const vec3 &ptB)
		{
			vec3 ptDiff = pt-ptA;
			vec3 ptDir = ptB-ptA;
			double t = ptDiff.dot(ptDir) / ptDir.dot(ptDir);
			*this = ptA + ptDir * t;
		}

		inline void vec3::clamp(double fMin, double fMax)
		{
			x = x < fMin ? fMin : x > fMax ? fMax : x;
			y = y < fMin ? fMin : y > fMax ? fMax : y;
			z = z < fMin ? fMin : z > fMax ? fMax : z;
		}


		//GAPI_MATH GAPI_Matrix operator * (const GAPI_Matrix &A, const GAPI_Matrix &B);
		//GAPI_MATH vec3 operator * (const GAPI_Matrix &A, const vec3 &B);
		//inline vec3 operator * (const vec3 &B, const GAPI_Matrix &A) { return operator *(A,B); }

		inline void vec3::from_string(const char *sz)
		{
		//	sscanf_s(sz,"%f,%f,%f",&x,&y,&z);
		}

		inline void vec3::to_string(char *szDest, int iMaxLen)
		{
		//	sprintf_s(szDest,iMaxLen,"%f,%f,%f",x,y,z);
		}

		/*
		inline vec2& operator = (vec2 &Dst, const vec3 &Src)
		{
			Dst.x = Src.x;
			Dst.y = Src.y;
			return Dst;
		}
		*/

	} // math

} // aspect

#endif // __MATH_vec3_HPP__