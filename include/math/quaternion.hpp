#ifndef MATH_QUATERNION_HPP_INCLUDED
#define MATH_QUATERNION_HPP_INCLUDED

#include "math/math.hpp"

namespace aspect { namespace math {

class matrix;
class euler_angles;

class MATH_API quat
{
public:

	union
	{
		struct
		{
			double x, y, z, w;
		};

		double v[4];
	};

public:

	enum ORD { X, Y, Z, W };

	// constructors...

	quat() { unity(); }

	explicit quat(const matrix& m) { from_matrix(m); }
	explicit quat(const euler_angles &ang);

	quat(double x, double y, double z, double w)
		: x(x), y(y), z(z), w(w)
	{
	}

	// accessors...
	double get(ORD ord) const { return v[ord]; }
	void set(ORD ord, double f) { v[ord] = f; }

	// operators...
	inline bool operator == (const quat& q) const;
	operator double*() { return(v); }

	// misc functions...

	// Multiply by scalar f. 
	void scale(double f);

	// Caculate conjugate of quaternion.
	void conjugate() { v[X]=-v[X]; v[Y]=-v[Y]; v[Z]=-v[Z]; }

	// Place on the unit sphere
	void unity() { ::memset(this, 0, sizeof(quat)); v[W]=1.0f; }
	void normalize();

	// Quaternion product this * qR.  Note: order is important!
	// To combine rotations, use the product qFirst.Multiply(qSecond),
	// which gives the effect of rotating by qFirst then qSecond.
	void multiply(const quat &qR);
	inline quat & operator *= (const quat &Q);
	inline quat   operator *  (const quat &Q) const;

	// Matrix conversion
	void to_matrix(matrix &m) const;
	matrix to_matrix() const;
	void from_matrix(const matrix &m);

	// Spherical linear interpolation of unit quaternions with spins.  A & B -> this
	void slerp(double alpha, const quat &a, const quat &b, int spin = 0);
};

inline bool quat::operator == (const quat& q) const
{
	return is_zero(q.v[X] - v[X])
		&& is_zero(q.v[Y] - v[Y])
		&& is_zero(q.v[Z] - v[Z])
		&& is_zero(q.v[W] - v[W]);
}

inline quat & quat::operator *= (const quat &Q)
{
	multiply(Q);
	return *this;
}

inline quat   quat::operator *  (const quat &Q) const
{
	quat qr(*this);
	qr.multiply(Q);
	return qr;
}

}} // aspect::math

namespace v8pp {

using aspect::get_option;
using aspect::set_option;

template<>
struct convert<aspect::math::quat>
{
	typedef aspect::math::quat result_type;

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

	static v8::Handle<v8::Value> to_v8(v8::Isolate* isolate, aspect::math::quat const& value)
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

#endif // MATH_QUATERNION_HPP_INCLUDED
