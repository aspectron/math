#ifndef MATH_VEC2_HPP_INCLUDED
#define MATH_VEC2_HPP_INCLUDED

#include "math/math.hpp"

#include <v8pp/convert.hpp>
#include <v8pp/object.hpp>

namespace aspect { namespace math {

class matrix;

class MATH_API vec2
{
public:

	union
	{
		struct
		{
			double x,y;
		};

		double v[2];
	};

	// constructors...

	vec2	() { x = double(0.0); y = double(0.0); } 
	vec2	(double _x, double _y) {x = _x; y = _y; }

	// initializers...
	void zero() { x = y = 0.0; }

	// accessors and casts...
	double &operator [] (int i) { return v[i]; }
	operator double*() { return(v); }

	void set_x(const double _x) { x = _x; }
	void set_y(const double _y) { y = _y; }

	double get_x() const { return x; }
	double get_y() const { return y; }

	// operators...

	inline vec2& operator *= (const vec2 &w);
	inline vec2& operator *= (const matrix &m);
	inline vec2& operator += (const vec2 &w);
	inline vec2& operator -= (const vec2 &w);
	inline vec2  operator -  (const vec2 &w) const;
	inline vec2  operator +  (const vec2 &w) const;
	inline vec2  operator *  (const vec2 &w) const;
	inline vec2  operator /  (const vec2 &w) const;
	inline vec2  operator -  (double f) const;
	inline vec2  operator +  (double f) const;
	inline vec2  operator *  (double f) const;

	inline vec2 operator -() const;
	inline vec2 operator +() const;

	bool operator == (const vec2 &w) const {return (w-*this).is_zero_length();}
	bool operator != (const vec2 &w) const {return !(w-*this).is_zero_length();}

	// misc functions...

	bool is_zero_length() const { return is_zero(x) && is_zero(y); }
	double get_length() const; // Compute the length of a vector
	double distance(const vec2& loc) const; // Compute the distance between two points
	double sqr_distance(const vec2& loc) const; // Compute SQR distance between two points
//	double get_angle(const vec2& loc) const; // Get the angle between 2 vectors

	// Get a point along line AB, at a distance from A expressed in relation to
	// the distance from A to B.  Use this to project(extrapolate) or interpolate a line.
//	void get_point_on_line(vec2& locResult, const double f, const vec2& locB) const;
	void get_point_on_line(vec2& locA, vec2& locB, const double f);

//	double dot(const vec2& loc) const;
	void cross(vec2& locOut, const vec2& loc) const;
	void normalize();

	void rotate(double f)
	{
		double fX = x * cos(f) - y * sin(f);
		double fY = x * sin(f) + y * cos(f);
		x = fX;
		y = fY;
	}

	inline double dot(const vec2& loc) const;
	inline double get_angle(const vec2& loc) const;
};

inline vec2& vec2::operator -= (const vec2 &w) { x -= w.x; y -= w.y; return *this; }
inline vec2& vec2::operator += (const vec2 &w) { x += w.x; y += w.y; return *this; }
inline vec2& vec2::operator *= (const vec2 &w) { x *= w.x; y *= w.y; return *this; }
inline vec2  vec2::operator -  (const vec2 &w) const { return vec2(x-w.x, y-w.y); }
inline vec2  vec2::operator +  (const vec2 &w) const { return vec2(x+w.x, y+w.y); }
inline vec2  vec2::operator *  (const vec2 &w) const { return vec2(x*w.x, y*w.y); }
inline vec2  vec2::operator /  (const vec2 &w) const { return vec2(x/w.x, y/w.y); }
inline vec2  vec2::operator -  (double f) const { return vec2(x-f, y-f); }
inline vec2  vec2::operator +  (double f) const { return vec2(x+f, y+f); }
inline vec2  vec2::operator *  (double f) const { return vec2(x*f, y*f); }

inline vec2  vec2::operator -() const { return (vec2(-x, -y)); }
inline vec2  vec2::operator +() const { return (*this); }

inline vec2  operator *(const double f, const vec2 &p) { return p * f; } //vec2(p.x*f, p.y*f); }


// vec2::fDistance()
//   Compute the distance between two locs
inline double vec2::distance(const vec2& loc) const
{
	return sqrt( (loc.x-x) * (loc.x-x)
				 +(loc.y-y) * (loc.y-y));
}


inline void vec2::get_point_on_line(vec2& locA, vec2& locB, const double f)
{
	x = locA.x + (f * (locB.x - locA.x));
	y = locA.y + (f * (locB.y - locA.y));
}

inline double vec2::get_length() const
{
	return sqrt(x*x + y*y);
}


inline double vec2::dot(const vec2& loc) const
{
	return x * loc.x + y * loc.y;
}

inline double vec2::get_angle(const vec2& loc) const
{
	double f = get_length() * loc.get_length();
	if (is_zero(f)) return 0.f;
	f = dot(loc) / f;
	if (f < -1.0) f = -1.0;
	if (f >  1.0) f =  1.0;
	return ::acos(f);
}


}} // aspect::math

template<>
struct v8pp::convert<aspect::math::vec2>
{
	using from_type = aspect::math::vec2;
	using to_type = v8::Handle<v8::Object>;

	static bool is_valid(v8::Isolate*, v8::Handle<v8::Value> value)
	{
		return value->IsObject();
	}

	static from_type from_v8(v8::Isolate* isolate, v8::Handle<v8::Value> value)
	{
		v8::HandleScope scope(isolate);

		from_type result;
		if (value->IsArray())
		{
			v8::Local<v8::Array> arr = value.As<v8::Array>();
			v8::Local<v8::Value> x = arr->Get(0), y = arr->Get(1);
			if (arr->Length() != 2 || !x->IsNumber() || !y->IsNumber())
			{
				throw std::invalid_argument("expected [x, y] array");
			}
			result.x = v8pp::from_v8<double>(isolate, x);
			result.y = v8pp::from_v8<double>(isolate, y);
		}
		else if (value->IsObject())
		{
			v8::Local<v8::Object> obj = value->ToObject();
			if (!get_option(isolate, obj, "x", result.x) || !get_option(isolate, obj, "y", result.y))
			{
				throw std::invalid_argument("expected {x, y} object");
			}
		}
		else
		{
			throw std::invalid_argument("expected [x, y] array or {x, y} object");
		}
		return result;
	}

	static to_type to_v8(v8::Isolate* isolate, aspect::math::vec2 const& value)
	{
		v8::EscapableHandleScope scope(isolate);

		v8::Local<v8::Object> obj = v8::Object::New(isolate);
		set_option(isolate, obj, "x", value.x);
		set_option(isolate, obj, "y", value.y);

		return scope.Escape(obj);
	}
};

#endif // MATH_VEC2_HPP_INCLUDED
