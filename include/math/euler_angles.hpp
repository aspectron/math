#ifndef MATH_EULER_ANGLES_HPP_INCLUDED
#define MATH_EULER_ANGLES_HPP_INCLUDED

namespace aspect { namespace math {

class quat;
class matrix;

//============================================================================
// Euler Angles
//
// Support for 24 angle schemes - Adapted from Ken Shoemake, 1993 
// Order type constants, constructors, extractors:
//    There are 24 possible conventions, designated by: 
//         o EulAxI = axis used initially			   
//         o EulPar = parity of axis permutation		   
//         o EulRep = repetition of initial axis as last
//         o EulFrm = frame from which axes are taken   
//    Axes I,J,K will be a permutation of X,Y,Z.		   
//    Axis H will be either I or K, depending on EulRep.
//    Frame S takes axes from initial static frame.	   
//	Example:
//    If ord = (AxI=X, Par=Even, Rep=No, Frm=S), then   
//    {a,b,c,ord} means Rz(c)Ry(b)Rx(a), where Rz(c)v   
//    rotates v around Z by c radians.				   


class MATH_API euler_angles    // (x,y,z)=ang 1,2,3, w=order code 
{
protected:
	double v[4];

public:
	// No virtual methods allowed!

	enum AXIS		{ X, Y, Z, W };
	enum FRAME		{ EulFrmS = 0, EulFrmR = 1 };
	enum REPETITION { EulRepNo = 0, EulRepYes = 1};
	enum PARITY		{ EulParEven = 0, EulParOdd = 1};

// MAKE_ORDER creates an order value between 0 and 23 from 4-tuple choices.
#define MAKE_ORDER(i,p,r,f)	   (((((((i)<<1)+(p))<<1)+(r))<<1)+(f))

	enum ORDER {
		// Static axes 
		orderXYZs =   MAKE_ORDER(X,EulParEven,EulRepNo,EulFrmS),
		orderXYXs =   MAKE_ORDER(X,EulParEven,EulRepYes,EulFrmS),
		orderXZYs =   MAKE_ORDER(X,EulParOdd,EulRepNo,EulFrmS),
		orderXZXs =   MAKE_ORDER(X,EulParOdd,EulRepYes,EulFrmS),
		orderYZXs =   MAKE_ORDER(Y,EulParEven,EulRepNo,EulFrmS),
		orderYZYs =   MAKE_ORDER(Y,EulParEven,EulRepYes,EulFrmS),
		orderYXZs =   MAKE_ORDER(Y,EulParOdd,EulRepNo,EulFrmS),
		orderYXYs =   MAKE_ORDER(Y,EulParOdd,EulRepYes,EulFrmS),
		orderZXYs =   MAKE_ORDER(Z,EulParEven,EulRepNo,EulFrmS),
		orderZXZs =   MAKE_ORDER(Z,EulParEven,EulRepYes,EulFrmS),
		orderZYXs =   MAKE_ORDER(Z,EulParOdd,EulRepNo,EulFrmS),
		orderZYZs =   MAKE_ORDER(Z,EulParOdd,EulRepYes,EulFrmS),
		// Rotating axes
		orderZYXr =   MAKE_ORDER(X,EulParEven,EulRepNo,EulFrmR),
		orderXYXr =   MAKE_ORDER(X,EulParEven,EulRepYes,EulFrmR),
		orderYZXr =   MAKE_ORDER(X,EulParOdd,EulRepNo,EulFrmR),
		orderXZXr =   MAKE_ORDER(X,EulParOdd,EulRepYes,EulFrmR),
		orderXZYr =   MAKE_ORDER(Y,EulParEven,EulRepNo,EulFrmR),
		orderYZYr =   MAKE_ORDER(Y,EulParEven,EulRepYes,EulFrmR),
		orderZXYr =   MAKE_ORDER(Y,EulParOdd,EulRepNo,EulFrmR),
		orderYXYr =   MAKE_ORDER(Y,EulParOdd,EulRepYes,EulFrmR),
		orderYXZr =   MAKE_ORDER(Z,EulParEven,EulRepNo,EulFrmR),
		orderZXZr =   MAKE_ORDER(Z,EulParEven,EulRepYes,EulFrmR),
		orderXYZr =   MAKE_ORDER(Z,EulParOdd,EulRepNo,EulFrmR),
		orderZYZr =   MAKE_ORDER(Z,EulParOdd,EulRepYes,EulFrmR)
	};
#undef MAKE_ORDER

	euler_angles() {}
	euler_angles(double fX, double fY, double fZ, ORDER order) { from_angles(fX, fY, fZ, order); }
	euler_angles(const quat& q, ORDER order) { from_quaternion(q, order); }
	euler_angles(const matrix &m, ORDER order) { from_matrix(m, order); }
		
	// Angle conversion
	void to_angles(double& fX, double& fY, double& fZ, int& order) const	{ fX=v[X]; fY=v[Y]; fZ=v[Z]; order=(int)v[W]; }
	void from_angles(double fX, double fY, double fZ, ORDER order)			{ v[X]=fX; v[Y]=fY; v[Z]=fZ; v[W]=(double)order; }

	// Quaternion conversion
	void to_quaternion(quat& q) const;
	void from_quaternion(const quat& q, ORDER order);

	// Matrix conversion
	void to_matrix(matrix &m) const;
	void from_matrix(const matrix &m, ORDER order);

	double get(int axis) const { return v[axis]; }
	void set(int axis, double f) { v[axis] = f; }
};

}} // aspect::math

#endif // MATH_EULER_ANGLES_HPP_INCLUDED
