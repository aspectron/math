#pragma once
#ifndef __MATH_vec2_HPP__
#define __MATH_vec2_HPP__

#include "math.hpp"

namespace aspect
{

namespace math
{

class matrix;


class vec2
{

	public:				

		union
		{
			struct
			{
				axScalar x,y;
			};

			axScalar v[2];
		};

		// constructors...

		vec2	() { x = axScalar(0.0); y = axScalar(0.0); } 
		vec2	(axScalar _x, axScalar _y) {x = _x; y = _y; }

		// initializers...

		void zero(void) { x = y = 0.0; }

		// accessors and casts...

		axScalar &operator [] (int i) { return v[i]; }
		operator axScalar*() { return(v); }

		void	set_x(const axScalar _x)    { x = _x; }
		void	set_y(const axScalar _y)    { y = _y; }

		axScalar   get_x() const			{ return x; }
		axScalar   get_y() const			{ return y; }

		// operators...

		inline vec2&	operator *= (const vec2 &w);			
		inline vec2&	operator *= (const matrix &m);			
		inline vec2&	operator += (const vec2 &w);			
		inline vec2&	operator -= (const vec2 &w);			
		inline vec2	operator -  (const vec2 &w)	const;	
		inline vec2	operator +  (const vec2 &w)	const;	
		inline vec2	operator *  (const vec2 &w)	const;	
		inline vec2	operator /  (const vec2 &w)	const;	
		inline vec2	operator -  (axScalar f)	const;	
		inline vec2	operator +  (axScalar f)	const;	
		inline vec2	operator *  (axScalar f)	const;	

		inline vec2 operator -() const;
		inline vec2 operator +() const;

		bool       operator == (const vec2 &w) const	{return (w-*this).is_zero_length();}
		bool       operator != (const vec2 &w) const	{return !(w-*this).is_zero_length();}

		// misc functions...

		bool		is_zero_length	() const	{ return is_zero(x) && is_zero(y); }
		axScalar		get_length		() const;							// Compute the length of a vector
		axScalar		distance		(const vec2& loc) const;		// Compute the distance between two points
		axScalar		sqr_distance	(const vec2& loc) const;		// Compute SQR distance between two points
//		axScalar		get_angle		(const vec2& loc) const;		// Get the angle between 2 vectors

		// Get a point along line AB, at a distance from A expressed in relation to
		// the distance from A to B.  Use this to project(extrapolate) or interpolate a line.
//		void		get_point_on_line	(vec2& locResult, const axScalar f, const vec2& locB) const;
		void get_point_on_line(vec2& locA, vec2& locB, const axScalar f);

//		axScalar		dot				(const vec2& loc) const;
		void		cross			(vec2& locOut, const vec2& loc) const;
		void		normalize		();

		void		rotate(axScalar f)	
		{	
			axScalar fX = x * cos(f) - y * sin(f); 
			axScalar fY = x * sin(f) + y * cos(f); 
			x = fX;
			y = fY;
		}

		inline axScalar dot(const vec2& loc) const;
		inline axScalar get_angle(const vec2& loc) const;


};

inline vec2&	vec2::operator -= (const vec2 &w)			{ x -= w.x; y -= w.y; return *this; }
inline vec2&	vec2::operator += (const vec2 &w)			{ x += w.x; y += w.y; return *this; }
inline vec2&	vec2::operator *= (const vec2 &w)			{ x *= w.x; y *= w.y; return *this; }
inline vec2	vec2::operator -  (const vec2 &w)	const	{ return vec2(x-w.x, y-w.y); }
inline vec2	vec2::operator +  (const vec2 &w)	const	{ return vec2(x+w.x, y+w.y); }
inline vec2	vec2::operator *  (const vec2 &w)	const	{ return vec2(x*w.x, y*w.y); }
inline vec2	vec2::operator /  (const vec2 &w)	const	{ return vec2(x/w.x, y/w.y); }
inline vec2	vec2::operator -  (axScalar f)				const	{ return vec2(x-f, y-f); }
inline vec2	vec2::operator +  (axScalar f)				const	{ return vec2(x+f, y+f); }
inline vec2	vec2::operator *  (axScalar f)				const	{ return vec2(x*f, y*f); }

inline vec2	vec2::operator -() const { return (vec2(-x, -y)); }
inline vec2	vec2::operator +() const { return (*this); }

inline vec2	operator *(const axScalar f, const vec2 &p) { return p * f; } //vec2(p.x*f, p.y*f); }


// vec2::fDistance()
//   Compute the distance between two locs
inline axScalar vec2::distance(const vec2& loc) const 
{
	return sqrt( (loc.x-x) * (loc.x-x)
				 +(loc.y-y) * (loc.y-y));
}


inline void vec2::get_point_on_line(vec2& locA, vec2& locB, const axScalar f)
{
	x = locA.x + (f * (locB.x - locA.x));
	y = locA.y + (f * (locB.y - locA.y));
}

inline axScalar vec2::get_length() const 
{
	return sqrt(x*x + y*y);
}


inline axScalar vec2::dot(const vec2& loc) const
{
	return x * loc.x + y * loc.y;
}

inline axScalar vec2::get_angle(const vec2& loc) const
{
	axScalar f = get_length() * loc.get_length();
	if (is_zero(f)) return 0.f;
	f = dot(loc) / f;
	if (f < -1.0) f = -1.0;
	if (f >  1.0) f =  1.0;
	return ::acos(f);
}



} // math

} // aspect

#endif // __MATH_vec2_HPP__