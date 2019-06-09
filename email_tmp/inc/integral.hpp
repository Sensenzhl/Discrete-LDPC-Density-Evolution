#ifndef _INTEGRAL_HPP_
#define _INTEGRAL_HPP_

#define INTEGRAL_ITERATION_SIZE                10000000             // iteration number of the integral
//#define PI                                     3.1415926535897932384626 
#define VARIANCE_DELTA                         0.00000001 
#define QUANTIZE_DELTA                         0.1
#define AWGN_CHANNEL                           true

#define DELTA_AWGN               			   10
#define DELTA_RAYLEIGH           			   10

#define VERBOSE           			           false

const int convolution_range = 20;
const int convolution_length = convolution_range * 2 / QUANTIZE_DELTA + 1;
const int g_nBufferSize = (convolution_range * 2 / QUANTIZE_DELTA + 1) * (convolution_range * 2 / QUANTIZE_DELTA + 1);

int convert_power_of_two(int length, bool verbose);
double pdf_Rayleigh_func(double x, double variance);
double pdf_awgn_func(double x, double variance);
double quantize_number(double w, double delta);
double *quantize_message(double(*f)(double, double), double variance, double delta, double min, double max);
bool normal_convolution(Complex inVec1[], Complex inVec2[], int vecLen, Complex outVec4[], bool verbose);
bool multiply_convolution(Complex const inVec1[], Complex const inVec2[], int const vecLen, double vector[convolution_length][convolution_length], Complex outVec4[], bool verbose);

#endif