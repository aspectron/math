#include "math/math.hpp"
#include "math/euler_angles.hpp"

namespace aspect { namespace math {

//============================================================================
// Euler Angles
// Support for 24 angle schemes - Adapted from Ken Shoemake, 1993 
//
#define EulFrm(ord)  ((unsigned)(ord)&1)
#define EulRep(ord)  (((unsigned)(ord)>>1)&1)
#define EulPar(ord)  (((unsigned)(ord)>>2)&1)

#define EulSafe	     "\000\001\002\000"
#define EulNext	     "\001\002\000\001"
#define EulAxI(ord)  ((int)(EulSafe[(((unsigned)(ord)>>3)&3)]))
#define EulAxJ(ord)  ((int)(EulNext[EulAxI(ord)+(EulPar(ord)==euler_angles::EulParOdd)]))
#define EulAxK(ord)  ((int)(EulNext[EulAxI(ord)+(EulPar(ord)!=euler_angles::EulParOdd)]))
#define EulAxH(ord)  ((EulRep(ord)==EulRepNo)?EulAxK(ord):EulAxI(ord))

// EulGetOrd unpacks all useful information about order simultaneously.
#define EulGetOrd(ord,i,j,k,h,n,s,f) {unsigned o=ord;f=o&1;o>>=1;s=o&1;o>>=1;\
	n=o&1;o>>=1;i=EulSafe[o&3];j=EulNext[i+n];k=EulNext[i+1-n];h=s?k:i;}

// Quaternion conversion:
// Construct quaternion from Euler angles (in radians).
void euler_angles::to_quaternion(quat& q) const
{
	double a[3], ti, tj, th, ci, cj, ch, si, sj, sh, cc, cs, sc, ss;
	int i,j,k,h,n,s,f;
	int ord = (int)v[W];
	EulGetOrd(ord,i,j,k,h,n,s,f);

	euler_angles ea(*this);
	if (f==euler_angles::EulFrmR) { ea.v[X]=v[Z]; ea.v[Z]=v[X]; }
	if (n==euler_angles::EulParOdd) ea.v[Y] = -v[Y];
	ti = ea.v[X]*double(0.5); tj = ea.v[Y]*double(0.5); th = ea.v[Z]*double(0.5);
	ci = cos(ti);  cj = cos(tj);  ch = cos(th);
	si = sin(ti);  sj = sin(tj);  sh = sin(th);
	cc = ci*ch; cs = ci*sh; sc = si*ch; ss = si*sh;

	if (s==euler_angles::EulRepYes) {
		a[i] = cj*(cs + sc);	// Could speed up with 
		a[j] = sj*(cc + ss);	// trig identities. 
		a[k] = sj*(cs - sc);
		q.set(quat::W, (double)(cj*(cc - ss)));
	} else {
		a[i] = cj*sc - sj*cs;
		a[j] = cj*ss + sj*cc;
		a[k] = cj*cs - sj*sc;
		q.set(quat::W, (double)(cj*cc + sj*ss));
	}
	if (n==euler_angles::EulParOdd) a[j] = -a[j];
	q.set(quat::X, (double)a[X]); 
	q.set(quat::Y, (double)a[Y]); 
	q.set(quat::Z, (double)a[Z]);
}

// Convert quaternion to Euler angles (in radians).
void euler_angles::from_quaternion(const quat& q, ORDER order)
{
	matrix m;
	q.to_matrix(m);
	from_matrix(m, order);
}

// Matrix conversion
void euler_angles::to_matrix(matrix &m) const
{
	double ti, tj, th, ci, cj, ch, si, sj, sh, cc, cs, sc, ss;
	int i,j,k,h,n,s,f;
	int ord = (int)v[W];
	EulGetOrd(ord,i,j,k,h,n,s,f);

	euler_angles ea(*this);
	if (f==euler_angles::EulFrmR) { ea.v[X] = v[Z]; ea.v[Z] = v[X];}
	if (n==euler_angles::EulParOdd) {ea.v[X] = -v[X]; ea.v[Y] = -v[Y]; ea.v[Z] = -v[Z];}
	ti = ea.v[X];	  tj = ea.v[Y];	th = ea.v[Z];
	ci = cos(ti); cj = cos(tj); ch = cos(th);
	si = sin(ti); sj = sin(tj); sh = sin(th);
	cc = ci*ch; cs = ci*sh; sc = si*ch; ss = si*sh;
	if (s==euler_angles::EulRepYes)
	{
		m[i][i] = (double)(cj);	  m[i][j] = (double)( sj*si);    m[i][k] = (double)( sj*ci);
		m[j][i] = (double)(sj*sh);  m[j][j] = (double)(-cj*ss+cc); m[j][k] = (double)(-cj*cs-sc);
		m[k][i] = (double)(-sj*ch); m[k][j] = (double)( cj*sc+cs); m[k][k] = (double)( cj*cc-ss);
	}
	else
	{
		m[i][i] = (double)(cj*ch);  m[i][j] = (double)(sj*sc-cs);  m[i][k] = (double)(sj*cc+ss);
		m[j][i] = (double)(cj*sh);  m[j][j] = (double)(sj*ss+cc);  m[j][k] = (double)(sj*cs-sc);
		m[k][i] = (double)(-sj);	  m[k][j] = (double)(cj*si);     m[k][k] = (double)(cj*ci);
	}
	m[W][X]=m[W][Y]=m[W][Z]=m[X][W]=m[Y][W]=m[Z][W]=0.0; m[W][W]=1.0;
}

void euler_angles::from_matrix(const matrix &m, ORDER order)
{
	int i,j,k,h,n,s,f;
	EulGetOrd(order,i,j,k,h,n,s,f);
	if (s==euler_angles::EulRepYes)
	{
		double sy = sqrt(m[i][j]*m[i][j] + m[i][k]*m[i][k]);
		if (!is_zero(sy)) 
		{
			v[X] = (atan2(m[i][j], m[i][k]));
			v[Y] = (atan2(sy, m[i][i]));
			v[Z] = (atan2(m[j][i], -m[k][i]));
		}
		else
		{
			v[X] = (atan2(-m[j][k], m[j][j]));
			v[Y] = (atan2(sy, m[i][i]));
			v[Z] = 0.f;
		}
	}
	else
	{
		double cy = sqrt(m[i][i]*m[i][i] + m[j][i]*m[j][i]);
		if (!is_zero(cy))
		{
			v[X] = (atan2(m[k][j], m[k][k]));
			v[Y] = (atan2(-m[k][i], cy));
			v[Z] = (atan2(m[j][i], m[i][i]));
		}
		else
		{
			v[X] = (atan2(-m[j][k], m[j][j]));
			v[Y] = (atan2(-m[k][i], cy));
			v[Z] = 0.f;
		}
	}

	if (n==euler_angles::EulParOdd) {v[X] = -v[X]; v[Y] = - v[Y]; v[Z] = -v[Z];}
	if (f==euler_angles::EulFrmR) {double t = v[X]; v[X] = v[Z]; v[Z] = t;}
	v[W] = (double)order;
}

}} // aspect::math
