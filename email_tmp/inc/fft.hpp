#ifndef _FFT_HPP_
#define _FFT_HPP_
#define MAX_MATRIX_SIZE                4194304             // 2048 * 2048
#define PI                             3.1415926535897932384626
#define MAX_VECTOR_LENGTH              10000             //

typedef struct Complex
{
	double rl;
	double im;
}Complex;
 

class Fft
{
	public:
		Fft();
		~Fft();
	 
	public:
		bool is_power_of_two(int num);
		int  get_computation_layers(int num);         // calculate the layers of computation needed for FFT

		bool fft(Complex const inVec[], int const len, Complex outVec[]);            // 基于蝶形算法的快速傅里叶变换
	    bool ifft(Complex const inVec[], int const len, Complex outVec[]);
	    bool multiply(Complex const inVec1[], Complex const inVec2[], int const vecLen, Complex outVec[]);
};

void OutputDebugStringA(const char *strOutputString);

#endif