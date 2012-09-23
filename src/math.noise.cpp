#include "math.hpp"

namespace aspect
{

	namespace math
	{

		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////

		float noise3(float* vec);


		float bias(float a, float b)
		{
			return (float)(pow(a, (float) (log(b) *-0.301029995664f)));
		}

		float gain(float a, float b)
		{
			float p =(float)( log(1. - b) *-0.301029995664f);

			if (a < .001)
				return 0.0f;
			else if (a > 0.999f)
				return 1.0f;
			if (a < 0.5f)
				return (float)(pow(2.0f * a, p) *0.5f);
			else
				return (float)(1.0f - pow(2.0f * (1.0f - a), p) *0.5f);
		}


		static void noise_init(void);
// 		float noise3(float* vec);
// 		float turbulence(float *v, float freq);



		float noise(float *v, int levels)
		{
			float sum = 0.0f;
			float l,f = 1.0f;

			float tmp[3];

			for(l=(float)levels;l>=1.0f;l-=1.0f) 
			{	
				tmp[0]=v[0]*f;
				tmp[1]=v[1]*f;
				tmp[2]=v[1]*f;

				sum += noise3(tmp)/(float)f;
				f*=2.0f;
			}

			if (l>0.0f)
			{
				v[0] *= f;
				v[1] *= f;
				v[2] *= f;
				sum+=l*noise3(v)/f;
			}

			return sum;
		}

		float turbulence(float *v, float freq)
		{
			float t, vec[3]; 

			for (t = 0.0f ; freq >= 1.0f ; freq *= 0.5f) {
				vec[0] = freq * v[0];
				vec[1] = freq * v[1];
				vec[2] = freq * v[2];
				t +=(float)( fabs(noise3(vec)) / freq);
			}
			return t;
		}

		/* noise functions over 1, 2, and 3 dimensions */

		#define B 0x100
		#define BM 0xff

		#define N 0x1000
		#define NP 12   /* 2^N */
		#define NM 0xfff

		static int p[B + B + 2];
		static float g3[B + B + 2][3];
		static float g2[B + B + 2][2];
		static float g1[B + B + 2];
		static int start = 1;


		#define s_curve(t) ( t * t * (3.0f - 2.0f * t) )

		#define lerp(t, a, b) ( a + t * (b - a) )

		#define setup(i,b0,b1,r0,r1)\
			t = vec[i] + N;\
			b0 = ((int)t) & BM;\
			b1 = (b0+1) & BM;\
			r0 = t - (int)t;\
			r1 = r0 - 1.0f;



		float noise3(float* vec)
		{
			int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
			float rx0, rx1, ry0, ry1, rz0, rz1, *q, sy, sz, a, b, c, d, t, u, v;
			int i, j;

			if (start) {
				start = 0;
				noise_init();
			}

			setup(0, bx0,bx1, rx0,rx1);
			setup(1, by0,by1, ry0,ry1);
			setup(2, bz0,bz1, rz0,rz1);

			i = p[ bx0 ];
			j = p[ bx1 ];

			b00 = p[ i + by0 ];
			b10 = p[ j + by0 ];
			b01 = p[ i + by1 ];
			b11 = p[ j + by1 ];

			t  = s_curve(rx0);
			sy = s_curve(ry0);
			sz = s_curve(rz0);

		#define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )

			q = g3[ b00 + bz0 ] ; u = at3(rx0,ry0,rz0);
			q = g3[ b10 + bz0 ] ; v = at3(rx1,ry0,rz0);
			a = lerp(t, u, v);

			q = g3[ b01 + bz0 ] ; u = at3(rx0,ry1,rz0);
			q = g3[ b11 + bz0 ] ; v = at3(rx1,ry1,rz0);
			b = lerp(t, u, v);

			c = lerp(sy, a, b);

			q = g3[ b00 + bz1 ] ; u = at3(rx0,ry0,rz1);
			q = g3[ b10 + bz1 ] ; v = at3(rx1,ry0,rz1);
			a = lerp(t, u, v);

			q = g3[ b01 + bz1 ] ; u = at3(rx0,ry1,rz1);
			q = g3[ b11 + bz1 ] ; v = at3(rx1,ry1,rz1);
			b = lerp(t, u, v);

			d = lerp(sy, a, b);

			return lerp(sz, c, d);
		}

		static void normalize2(float v[2])
		{
			float s;

			s = (float)1.0f/(float)sqrt(v[0] * v[0] + v[1] * v[1]);
			v[0] = v[0] * s;
			v[1] = v[1] * s;
		}

		static void normalize3(float v[3])
		{
			float s;

			s =(float)1.0f/(float) sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
			v[0] = v[0] * s;
			v[1] = v[1] * s;
			v[2] = v[2] * s;
		}

		#define random() rand()

		static void noise_init(void)
		{
			int i, j, k;
			srand(0);

			for (i = 0 ; i < B ; i++) {
				p[i] = i;

				for (j = 0 ; j < 3 ; j++)
					g3[i][j] = (float)((random() % (B + B)) - B) / B;

				normalize3(g3[i]);
			}

			while (--i) {
				k = p[i];
				p[i] = p[j = random() % B];
				p[j] = k;
			}

			for (i = 0 ; i < B + 2 ; i++) {
				p[B + i] = p[i];

				for (j = 0 ; j < 3 ; j++)
					g3[B + i][j] = g3[i][j];
			}
		}



//////////////////////////////////////////////////////////////////////////




		/* This macro is a *lot* faster than using (long)floor() on an x86 CPU.
		   It actually speeds up the entire Worley() call with almost 10%.
		   Added by Stefan Gustavson, October 2003. */
		#define LFLOOR(x) ((x)<0 ? ((long)x-1) : ((long)x) )

		/* A hardwired lookup table to quickly determine how many feature
		   points should be in each spatial cube. We use a table so we don't
		   need to make multiple slower tests.  A random number indexed into
		   this array will give an approximate Poisson distribution of mean
		   density 2.5. Read the book for the longwinded explanation. */

		static int Poisson_count[256]=
		{4,3,1,1,1,2,4,2,2,2,5,1,0,2,1,2,2,0,4,3,2,1,2,1,3,2,2,4,2,2,5,1,2,3,2,2,2,2,2,3,
		 2,4,2,5,3,2,2,2,5,3,3,5,2,1,3,3,4,4,2,3,0,4,2,2,2,1,3,2,2,2,3,3,3,1,2,0,2,1,1,2,
		 2,2,2,5,3,2,3,2,3,2,2,1,0,2,1,1,2,1,2,2,1,3,4,2,2,2,5,4,2,4,2,2,5,4,3,2,2,5,4,3,
		 3,3,5,2,2,2,2,2,3,1,1,4,2,1,3,3,4,3,2,4,3,3,3,4,5,1,4,2,4,3,1,2,3,5,3,2,1,3,1,3,
		 3,3,2,3,1,5,5,4,2,2,4,1,3,4,1,5,3,3,5,3,4,3,2,2,1,1,1,1,1,2,4,5,4,5,4,2,1,5,1,1,
		 2,3,3,3,2,5,2,3,3,2,0,2,1,1,4,2,1,3,2,1,2,2,3,2,5,5,3,4,5,5,2,4,4,5,3,2,2,2,1,4,
		 2,3,3,4,2,5,4,2,4,2,2,2,4,5,3,2};

		/* This constant is manipulated to make sure that the mean value of F[0]
		   is 1.0. This makes an easy natural "scale" size of the cellular features. */
		#define DENSITY_ADJUSTMENT  0.398150

		/* the function to merge-sort a "cube" of samples into the current best-found
		   list of values. */
		static void AddSamples(long xi, long yi, long zi, long max_order,
					double at[3], double *F,
					double (*delta)[3], unsigned long *ID);


		/* The main function! */
		void cellular(double at[3], long max_order,
				double *F, double (*delta)[3], unsigned long *ID)
		{
		  double x2,y2,z2, mx2, my2, mz2;
		  double new_at[3];
		  long int_at[3], i;
		  
		  /* Initialize the F values to "huge" so they will be replaced by the
			 first real sample tests. Note we'll be storing and comparing the
			 SQUARED distance from the feature points to avoid lots of slow
			 sqrt() calls. We'll use sqrt() only on the final answer. */
		  for (i=0; i<max_order; i++) F[i]=999999.9;
		  
		  /* Make our own local copy, multiplying to make mean(F[0])==1.0  */
		  new_at[0]=DENSITY_ADJUSTMENT*at[0];
		  new_at[1]=DENSITY_ADJUSTMENT*at[1];
		  new_at[2]=DENSITY_ADJUSTMENT*at[2];

		  /* Find the integer cube holding the hit point */
		  int_at[0]=LFLOOR(new_at[0]); /* The macro makes this part a lot faster */
		  int_at[1]=LFLOOR(new_at[1]);
		  int_at[2]=LFLOOR(new_at[2]);

		  /* A simple way to compute the closest neighbors would be to test all
			 boundary cubes exhaustively. This is simple with code like: 
			 {
			   long ii, jj, kk;
			   for (ii=-1; ii<=1; ii++) for (jj=-1; jj<=1; jj++) for (kk=-1; kk<=1; kk++)
			   AddSamples(int_at[0]+ii,int_at[1]+jj,int_at[2]+kk, 
			   max_order, new_at, F, delta, ID);
			 }
			 But this wastes a lot of time working on cubes which are known to be
			 too far away to matter! So we can use a more complex testing method
			 that avoids this needless testing of distant cubes. This doubles the 
			 speed of the algorithm. */

		  /* Test the central cube for closest point(s). */
		  AddSamples(int_at[0], int_at[1], int_at[2], max_order, new_at, F, delta, ID);

		  /* We test if neighbor cubes are even POSSIBLE contributors by examining the
			 combinations of the sum of the squared distances from the cube's lower 
			 or upper corners.*/
		  x2=new_at[0]-int_at[0];
		  y2=new_at[1]-int_at[1];
		  z2=new_at[2]-int_at[2];
		  mx2=(1.0-x2)*(1.0-x2);
		  my2=(1.0-y2)*(1.0-y2);
		  mz2=(1.0-z2)*(1.0-z2);
		  x2*=x2;
		  y2*=y2;
		  z2*=z2;
		  
		  /* Test 6 facing neighbors of center cube. These are closest and most 
			 likely to have a close feature point. */
		  if (x2<F[max_order-1])  AddSamples(int_at[0]-1, int_at[1]  , int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (y2<F[max_order-1])  AddSamples(int_at[0]  , int_at[1]-1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (z2<F[max_order-1])  AddSamples(int_at[0]  , int_at[1]  , int_at[2]-1, 
							 max_order, new_at, F, delta, ID);
		  
		  if (mx2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]  , int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (my2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]+1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (mz2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]  , int_at[2]+1, 
							 max_order, new_at, F, delta, ID);
		  
		  /* Test 12 "edge cube" neighbors if necessary. They're next closest. */
		  if ( x2+ y2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]-1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if ( x2+ z2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]  , int_at[2]-1, 
							 max_order, new_at, F, delta, ID);
		  if ( y2+ z2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]-1, int_at[2]-1, 
							 max_order, new_at, F, delta, ID);  
		  if (mx2+my2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]+1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (mx2+mz2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]  , int_at[2]+1, 
							 max_order, new_at, F, delta, ID);
		  if (my2+mz2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]+1, int_at[2]+1, 
							 max_order, new_at, F, delta, ID);  
		  if ( x2+my2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]+1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if ( x2+mz2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]  , int_at[2]+1, 
							 max_order, new_at, F, delta, ID);
		  if ( y2+mz2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]-1, int_at[2]+1, 
							 max_order, new_at, F, delta, ID);  
		  if (mx2+ y2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]-1, int_at[2]  , 
							 max_order, new_at, F, delta, ID);
		  if (mx2+ z2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]  , int_at[2]-1, 
							 max_order, new_at, F, delta, ID);
		  if (my2+ z2<F[max_order-1]) AddSamples(int_at[0]  , int_at[1]+1, int_at[2]-1, 
							 max_order, new_at, F, delta, ID);  
		  
		  /* Final 8 "corner" cubes */
		  if ( x2+ y2+ z2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]-1, int_at[2]-1, 
								 max_order, new_at, F, delta, ID);
		  if ( x2+ y2+mz2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]-1, int_at[2]+1, 
								 max_order, new_at, F, delta, ID);
		  if ( x2+my2+ z2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]+1, int_at[2]-1, 
								 max_order, new_at, F, delta, ID);
		  if ( x2+my2+mz2<F[max_order-1]) AddSamples(int_at[0]-1, int_at[1]+1, int_at[2]+1, 
								 max_order, new_at, F, delta, ID);
		  if (mx2+ y2+ z2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]-1, int_at[2]-1, 
								 max_order, new_at, F, delta, ID);
		  if (mx2+ y2+mz2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]-1, int_at[2]+1, 
								 max_order, new_at, F, delta, ID);
		  if (mx2+my2+ z2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]+1, int_at[2]-1, 
								 max_order, new_at, F, delta, ID);
		  if (mx2+my2+mz2<F[max_order-1]) AddSamples(int_at[0]+1, int_at[1]+1, int_at[2]+1, 
								 max_order, new_at, F, delta, ID);
		  
		  /* We're done! Convert everything to right size scale */
		  for (i=0; i<max_order; i++)
			{
			  F[i]=sqrt(F[i])*(1.0/DENSITY_ADJUSTMENT);      
			  delta[i][0]*=(1.0/DENSITY_ADJUSTMENT);
			  delta[i][1]*=(1.0/DENSITY_ADJUSTMENT);
			  delta[i][2]*=(1.0/DENSITY_ADJUSTMENT);
			}
		  
		  return;
		}



		static void AddSamples(long xi, long yi, long zi, long max_order,
					   double at[3], double *F,
					   double (*delta)[3], unsigned long *ID)
		{
		  double dx, dy, dz, fx, fy, fz, d2;
		  long count, i, j, index;
		  unsigned long seed, this_id;
		  
		  /* Each cube has a random number seed based on the cube's ID number.
			 The seed might be better if it were a nonlinear hash like Perlin uses
			 for noise but we do very well with this faster simple one.
			 Our LCG uses Knuth-approved constants for maximal periods. */
		  seed=702395077*xi + 915488749*yi + 2120969693*zi;
		  
		  /* How many feature points are in this cube? */
		  count=Poisson_count[seed>>24]; /* 256 element lookup table. Use MSB */

		  seed=1402024253*seed+586950981; /* churn the seed with good Knuth LCG */

		  for (j=0; j<count; j++) /* test and insert each point into our solution */
			{
			  this_id=seed;
			  seed=1402024253*seed+586950981; /* churn */

			  /* compute the 0..1 feature point location's XYZ */
			  fx=(seed+0.5)*(1.0/4294967296.0); 
			  seed=1402024253*seed+586950981; /* churn */
			  fy=(seed+0.5)*(1.0/4294967296.0);
			  seed=1402024253*seed+586950981; /* churn */
			  fz=(seed+0.5)*(1.0/4294967296.0);
			  seed=1402024253*seed+586950981; /* churn */

			  /* delta from feature point to sample location */
			  dx=xi+fx-at[0]; 
			  dy=yi+fy-at[1];
			  dz=zi+fz-at[2];
		      
			  /* Distance computation!  Lots of interesting variations are
			 possible here!
			 Biased "stretched"   A*dx*dx+B*dy*dy+C*dz*dz
			 Manhattan distance   fabs(dx)+fabs(dy)+fabs(dz)
			 Radial Manhattan:    A*fabs(dR)+B*fabs(dTheta)+C*dz
			 Superquadratic:      pow(fabs(dx), A) + pow(fabs(dy), B) + pow(fabs(dz),C)
			 
			 Go ahead and make your own! Remember that you must insure that
			 new distance function causes large deltas in 3D space to map into
			 large deltas in your distance function, so our 3D search can find
			 them! [Alternatively, change the search algorithm for your special
			 cases.]       
			  */

			  d2=dx*dx+dy*dy+dz*dz; /* Euclidian distance, squared */

			  if (d2<F[max_order-1]) /* Is this point close enough to rememember? */
			{
			  /* Insert the information into the output arrays if it's close enough.
				 We use an insertion sort.  No need for a binary search to find
				 the appropriate index.. usually we're dealing with order 2,3,4 so
				 we can just go through the list. If you were computing order 50
				 (wow!!) you could get a speedup with a binary search in the sorted
				 F[] list. */
			  
			  index=max_order;
			  while (index>0 && d2<F[index-1]) index--;

			  /* We insert this new point into slot # <index> */

			  /* Bump down more distant information to make room for this new point. */
			  for (i=max_order-2; i>=index; i--)
				{
				  F[i+1]=F[i];
				  ID[i+1]=ID[i];
				  delta[i+1][0]=delta[i][0];
				  delta[i+1][1]=delta[i][1];
				  delta[i+1][2]=delta[i][2];
				}		
			  /* Insert the new point's information into the list. */
			  F[index]=d2;
			  ID[index]=this_id;
			  delta[index][0]=dx;
			  delta[index][1]=dy;
			  delta[index][2]=dz;
			}
			}
		  
		  return;
		}




	} // math

} // aspect