#ifndef _MATH_NOISE_HPP_
#define _MATH_NOISE_HPP_

namespace aspect
{
	namespace math
	{

		float noise3(float* vec);
		float turbulence(float *v, float freq);
		void cellular(double at[3], long max_order,
			double *F, double (*delta)[3], unsigned long *ID);

	} // math

} // aspect

#endif // _MATH_NOISE_HPP_