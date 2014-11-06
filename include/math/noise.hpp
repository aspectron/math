#ifndef MATH_NOISE_HPP_INCLUDED
#define MATH_NOISE_HPP_INCLUDED

namespace aspect { namespace math {

float noise3(float* vec);
float turbulence(float *v, float freq);
void cellular(double at[3], long max_order, double *F, double (*delta)[3], unsigned long *ID);

}} // aspect::math

#endif // MATH_NOISE_HPP_INCLUDED
