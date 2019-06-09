#ifndef _PDF_MATRIX_HPP_
#define _PDF_MATRIX_HPP_

#define MAX_RANGE              20
#define MIN_RANGE              -20
#define QUANTIZE_DELTA         0.1
#define NUMBER_OF_INTERATION   50

typedef struct Complex
{
	double rl;
	double im;
}Complex;

typedef struct PDF_Element_2
{
	int    max = MAX_RANGE;
	int    min = MIN_RANGE;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	//double pdf_rayleigh[length];
	//double pdf_awgn[length];

	double *pdf_rayleigh = new double[length];
	double *pdf_awgn = new double[length];

	//pdf_awgn     = quantize_message(PDF_AWGN_Func, 1, delta, min, max);
	//pdf_rayleigh = quantize_message(PDF_Rayleigh_Func, 1, delta, min, max);

}PDF_Element_2;


typedef struct PDF_Element_RAYLEIGH
{
	int    max = MAX_RANGE;
	int    min = MIN_RANGE;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double *pdf = new double[length];

	//pdf = quantize_message(PDF_Rayleigh_Func, 1, delta, min, max);

}PDF_Element_RAYLEIGH;


typedef struct PDF_Element_AWGN
{
	int    max = MAX_RANGE;
	int    min = MIN_RANGE;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double *pdf = new double[length];
	//double pdf[length];

	//pdf = quantize_message(PDF_AWGN_Func, 1, delta, min, max);

}PDF_Element_AWGN;


typedef struct PDF_Element
{
	int    max = MAX_RANGE;
	int    min = MIN_RANGE;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double *pdf = new double[length];

}PDF_Element;


typedef struct PDF_Element_INITIAL
{
	int    max = MAX_RANGE;
	int    min = MIN_RANGE;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double *pdf = new double[length];

}PDF_Element_INITIAL;


class PDF_MATRIX
{
	private:
		int Row;
		int Col;		
		int matrix_Size;
		bool is_AWGN;

		PDF_Element* m_lpBuf;
		//PDF_Element m_lpBuf[Row * Col];
		//PDF_Element* m_lpBuf;

		//PDF_Element_AWGN* m_lpBuf_AWGN;
		//PDF_Element_RAYLEIGH* m_lpBuf_RAYLEIGH;
	 
	public:
		//bool is_power_of_two(int num);
		//int  get_computation_layers(int num);         // calculate the layers of computation needed for FFT
		//PDF_MATRIX(bool is_AWGN);
		PDF_MATRIX();
		PDF_MATRIX(int row, int col);
		~PDF_MATRIX();

		PDF_Element operator()(int i, int j, int Col, bool is_AWGN);  //Operate like MATLAB MATRIX
		PDF_MATRIX & operator = (PDF_MATRIX& obj);

		int getRow(){return Row;}
		int getCol(){return Col;}

		PDF_Element *getPDF(){return m_lpBuf;}
		//bool getPDF_array(int i, int j);
		PDF_Element get_element(int i, int j);
		double *getPDF_array(int i, int j);
		bool Pop(PDF_Element *in_lpBuf){this->m_lpBuf = in_lpBuf; return true;}
		bool Push(double *pdf, int i, int j, int Row, int Col, int cnt, int length, bool verbose);
		bool clear(bool verbose);
		//PDF_Element_AWGN operator ()(int i, int j);

	    //bool pdf_convolution_RAYLEIGH(PDF_Element_RAYLEIGH const inPDF1[], PDF_Element_RAYLEIGH const inPDF2[], int const vecLen, PDF_Element_RAYLEIGH outPDF[]);
	    //bool pdf_convolution_AWGN(PDF_Element_AWGN const inPDF1[], PDF_Element_AWGN const inPDF2[], int const vecLen, PDF_Element_AWGN outVec[]);
	    //bool pdf_convolution(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN);
	    //bool pdf_convolution_RAYLEIGH(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN);
	    //bool pdf_convolution_AWGN(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN);
		//bool pdf_multiply_convolution(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN);	  

		//bool pdf_convolution(PDF_Element const in_PDF1, PDF_Element in_PDF2, PDF_Element out_PDF, bool verbose);
		PDF_Element pdf_convolution(PDF_Element in_PDF1, PDF_Element in_PDF2, bool verbose);
	    bool pdf_convolution_RAYLEIGH(PDF_Element const in_PDF1, PDF_Element const in_PDF2, PDF_Element out_PDF, bool verbose);
	    bool pdf_convolution_AWGN(PDF_Element const in_PDF1, PDF_Element const in_PDF2, PDF_Element out_PDF, bool verbose); 
		bool pdf_multiply_convolution(PDF_Element const in_PDF1, PDF_Element const in_PDF2, PDF_Element out_PDF, bool verbose);

	    bool matrix_row_update(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose);
	    bool matrix_column_update(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose);
};


void OutputDebugStringA(const char *strOutputString);
bool matrix_row_update(PDF_MATRIX in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose);
bool matrix_column_update(PDF_MATRIX in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose);
bool RCA_predict(PDF_MATRIX HB, double SNR, double p, double R[], double snr_R[], double infor, bool verbose);

#endif