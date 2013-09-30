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

namespace v8pp {

namespace detail {

template<>
struct from_v8<::aspect::math::vec4>
{
	typedef ::aspect::math::vec4 result_type;

	static result_type exec(v8::Handle<v8::Value> value)
	{
		v8::HandleScope scope;

		result_type result;
		if (value->IsArray())
		{
			v8::Handle<v8::Array> arr = value.As<v8::Array>();
		
			if (arr->Length() != 4)
			{
				throw std::invalid_argument("array must contain 4 coordinates");
			}
			result.x = v8pp::from_v8<double>(arr->Get(0));
			result.y = v8pp::from_v8<double>(arr->Get(1));
			result.z = v8pp::from_v8<double>(arr->Get(2));
			result.w = v8pp::from_v8<double>(arr->Get(3));
		}
		else if (value->IsObject())
		{
			v8::Handle<v8::Object> obj = value->ToObject();
			result.x = v8pp::from_v8<double>(obj->Get(to_v8("x")));
			result.y = v8pp::from_v8<double>(obj->Get(to_v8("y")));
			result.z = v8pp::from_v8<double>(obj->Get(to_v8("z")));
			result.w = v8pp::from_v8<double>(obj->Get(to_v8("w")));
		}
		else
		{
			throw std::invalid_argument("expecting object(x,y,z,w) or array[4]");
		}
		return result;
	}
};

template<typename U>
struct from_v8_ref<::aspect::math::vec4, U> : from_v8<::aspect::math::vec4> {};

} // detail

inline v8::Handle<v8::Value> to_v8(::aspect::math::vec4 const& value)
{
	v8::HandleScope scope;

	v8::Handle<v8::Object> obj = v8::Object::New();

	obj->Set(to_v8("x"), to_v8(value.x));
	obj->Set(to_v8("y"), to_v8(value.y));
	obj->Set(to_v8("z"), to_v8(value.z));
	obj->Set(to_v8("w"), to_v8(value.w));

	return scope.Close(obj);
}

} // v8pp

#endif // __MATH_vec4_HPP__
