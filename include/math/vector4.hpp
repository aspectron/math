#ifndef MATH_VEC4_HPP_INCLUDED
#define MATH_VEC4_HPP_INCLUDED

#include "math/math.hpp"
#include "math/vector3.hpp"

namespace aspect { namespace math {

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

	inline vec4  operator *  (double f) const { return vec4(x*f, y*f, z*f, w*f); }
	inline vec4  operator /  (double f) const { return vec4(x/f, y/f, z/f, w/f); }
	inline vec4  operator - (const vec4 &src) const { return vec4(x-src.x,y-src.y,z-src.z,w-src.w); }
	inline vec4  operator + (const vec4 &src) const { return vec4(x+src.x,y+src.y,z+src.z,w+src.w); }
	inline vec4& operator += (const vec4 &src) { x += src.x; y += src.y; z += src.z; w += src.w; return *this; }

	double right(void) const { return left+width; }
	double bottom(void) const { return top+height; }
};

inline bool operator == (vec4 const& a, vec4 const& b)
{
	return (a.x-b.x) < math::precision_limit
		&& (a.y-b.y) < math::precision_limit
		&& (a.z-b.z) < math::precision_limit
		&& (a.w-b.w) < math::precision_limit;
}

}} // aspect::math

namespace v8pp {

using aspect::get_option;
using aspect::set_option;

template<>
struct convert<aspect::math::vec4>
{
	typedef aspect::math::vec4 result_type;

	static bool is_valid(v8::Isolate*, v8::Handle<v8::Value> value)
	{
		return value->IsObject();
	}

	static result_type from_v8(v8::Isolate* isolate, v8::Handle<v8::Value> value)
	{
		v8::HandleScope scope(isolate);

		result_type result;
		if (value->IsArray())
		{
			v8::Local<v8::Array> arr = value.As<v8::Array>();
			v8::Local<v8::Value> x = arr->Get(0), y = arr->Get(1);
			v8::Local<v8::Value> z = arr->Get(2), w = arr->Get(3);
			if (arr->Length() != 4 || !x->IsNumber() || !y->IsNumber()
				|| !z->IsNumber() || !w->IsNumber())
			{
				throw std::invalid_argument("expected [x, y, z, w] array");
			}
			result.x = v8pp::from_v8<double>(isolate, x);
			result.y = v8pp::from_v8<double>(isolate, y);
			result.z = v8pp::from_v8<double>(isolate, z);
			result.w = v8pp::from_v8<double>(isolate, w);
		}
		else if (value->IsObject())
		{
			v8::Local<v8::Object> obj = value->ToObject();
			if (!get_option(isolate, obj, "x", result.x)
				|| !get_option(isolate, obj, "y", result.y)
				|| !get_option(isolate, obj, "z", result.z)
				|| !get_option(isolate, obj, "w", result.w))
			{
				throw std::invalid_argument("expected {x, y, z, w} object");
			}
		}
		else
		{
			throw std::invalid_argument("expected [x, y, z, w] array or {x, y, z, w} object");
		}
		return result;
	}

	static v8::Handle<v8::Value> to_v8(v8::Isolate* isolate, aspect::math::vec4 const& value)
	{
		v8::EscapableHandleScope scope(isolate);

		v8::Local<v8::Object> obj = v8::Object::New(isolate);

		set_option(isolate, obj, "x", value.x);
		set_option(isolate, obj, "y", value.y);
		set_option(isolate, obj, "z", value.z);
		set_option(isolate, obj, "w", value.w);

		return scope.Escape(obj);
	}
};

} // v8pp

#endif // MATH_VEC4_HPP_INCLUDED
