#include <iostream> //standar library
#include <stdio.h>
#include <string.h>
#include "math.h"   //sorce code catalogue
#include "pdf_matrix.hpp"
#include "integral.hpp"

using namespace std;

#define ELEMENT_COUNT               (96 * 180)
//#define ELEMENT_COUNT               (84 * 180)
//#define ELEMENT_COUNT               (60 * 180)

//PDF_MATRIX::PDF_MATRIX(bool isAWGN)
PDF_MATRIX::PDF_MATRIX()
{
	Row = 0;
	Col = 0;
	matrix_Size = Row * Col;
	//matrix_Size = ELEMENT_COUNT;
	is_AWGN = true;

	this->m_lpBuf = new PDF_Element[matrix_Size];
	memset(m_lpBuf,0,Row*Col*sizeof(PDF_Element)); // assign 0 for this memory
}


PDF_MATRIX::PDF_MATRIX(int row, int col)
{
	Row = row;
	Col = col;
	matrix_Size = Row * Col;
	//matrix_Size = ELEMENT_COUNT;
	is_AWGN = true;

	this->m_lpBuf = new PDF_Element[matrix_Size];
	memset(m_lpBuf,0,Row*Col*sizeof(PDF_Element));
}


PDF_MATRIX::~PDF_MATRIX()
{
	Row = 0;
	Col = 0;

	//if((m_lpBuf_AWGN != NULL) || (m_lpBuf_RAYLEIGH != NULL))
	//{
	//	delete[] this->m_lpBuf_AWGN;
	//	delete[] this->m_lpBuf_RAYLEIGH;
	//}
	delete[] this->m_lpBuf;
}


//PDF_Element & PDF_MATRIX::operator ()(int i, int j, int Col, bool is_AWGN)
PDF_Element PDF_MATRIX::operator ()(int i, int j, int Col, bool is_AWGN)
{
    // Unsafe
    //if is_AWGN
    //	return *(m_lpBuf_AWGN + ((i-1) * Col + (j-1)) * sizeof(PDF_Element_AWGN));
    //else
    //	return *(m_lpBuf_RAYLEIGH + ((i-1) * Col + (j-1)) * sizeof(PDF_Element_RAYLEIGH));
    return *(m_lpBuf + ((i-1) * Col + (j-1)) * sizeof(PDF_Element));
}


PDF_MATRIX & PDF_MATRIX::operator = (PDF_MATRIX& obj)
{
	printf(" = Copying0 \n"); 
    // Unsafe
    int matrix_Size = Row * Col;
	double  min = -DELTA_RAYLEIGH;
	double  max = DELTA_RAYLEIGH;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;   

	printf(" matrix_Size is %d and length is %d \n",matrix_Size,length);

    if (m_lpBuf != NULL)
    {
    	/*
    	for (int k = 0; k < matrix_Size; ++k)
    	{
    		delete m_lpBuf[k].pdf; 
    	}
    	*/

    	printf(" = Copying0.3 \n");	
    	
        //delete [] m_lpBuf;
    }
    // 赋值stack内存的值

    printf(" = Copying0.5 \n");	

	printf(" = Copying1 \n"); 

    this-> Row = obj.Row;
    this-> Col = obj.Col;
    this-> matrix_Size = obj.matrix_Size;
    this-> is_AWGN = obj.is_AWGN;

    printf(" = Copying2 \n");

    PDF_Element* m_lpBuf = new PDF_Element[matrix_Size];

    for (int i = 0; i < matrix_Size; ++i)
    {
    	m_lpBuf[i].length = obj.m_lpBuf[i].length;
    	for (int j = 0; j < length; ++j)
    	{
    		m_lpBuf[i].pdf[j] = obj.m_lpBuf[i].pdf[j];
    	}    	
    }
    
    return *this;
}


PDF_Element PDF_MATRIX::get_element(int i, int j)
{
    return m_lpBuf[i * Col + j];
}


double *PDF_MATRIX::getPDF_array(int i, int j)
{
	double  min = -DELTA_RAYLEIGH;
	double  max = DELTA_RAYLEIGH;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;

	double  *tmp_pdf = new double [length];

    //return true;
    return m_lpBuf[i * 60 + j].pdf;
}


bool PDF_MATRIX::Push(double element_pdf[], int i, int j, int Row, int Col, int index, int length, bool verbose)
{
	int idx;

	if(!index)
		idx = i * Col + j;
	else
		idx = index;

	if(verbose)
	{
		for (int l = 0; l < length; ++l)	
		{
			printf("element_pdf is %f\n",element_pdf[l]);
		}

		printf("length is %d\n",length);
	}

	this->m_lpBuf[idx].length = length;

	//printf("/// 1 /// \n");

	this->m_lpBuf[idx].pdf = element_pdf;
		/* code */
	//printf("/// 2 /// \n");

	return true;
}


bool PDF_MATRIX::clear(bool verbose)
{
	Row = 0;
	Col = 0;	

	double  min = -DELTA_RAYLEIGH;
	double  max = DELTA_RAYLEIGH;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;	

	for (int i = 0; i < matrix_Size; ++i)
	{
		for (int j = 0; j < length; ++j)
		{
			this->m_lpBuf[i].pdf[j] = 0;
			/* code */
		}
	}

	matrix_Size = 0;

	if(verbose)
	{
		printf("clear done %f\n");
	}

	return true;
}


bool PDF_MATRIX::pdf_convolution_RAYLEIGH(PDF_Element in_PDF1, PDF_Element in_PDF2, PDF_Element out_PDF, bool verbose)
{
	double  min = -DELTA_RAYLEIGH;
	double  max = DELTA_RAYLEIGH;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	int  length_power_of_two;

	length_power_of_two = convert_power_of_two(length, verbose);

	double *pdf1 = new double[length_power_of_two];
	double *pdf2 = new double[length_power_of_two];	
	double *pdf_out = new double[length_power_of_two];	

	pdf1 = in_PDF1.pdf;
	pdf2 = in_PDF2.pdf;

	pdf2 = quantize_message(pdf_Rayleigh_func, 1, delta, min, max);

	Complex *inVec1 = new Complex[length_power_of_two];
	Complex *inVec2 = new Complex[length_power_of_two];
	Complex *outVec4 = new Complex[length_power_of_two];

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec1[i].rl = pdf2[i];
		else
			inVec1[i].rl = 0;
	}

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec2[i].rl = pdf2[i];
		else
			inVec2[i].rl = 0;
	}

	normal_convolution(inVec1, inVec2, length_power_of_two, outVec4, verbose);

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			pdf_out[i] = outVec4[i].rl;
		else
			pdf_out[i] = 0;
	}

	delete [] pdf1;
	delete [] pdf2;
	delete [] pdf_out;

	delete [] inVec1;
	delete [] inVec2;
	delete [] outVec4;

	return true;
}


bool PDF_MATRIX::pdf_convolution_AWGN(PDF_Element in_PDF1, PDF_Element in_PDF2, PDF_Element out_PDF, bool verbose)
{
	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;

	int  length_power_of_two;

	length_power_of_two = convert_power_of_two(length, verbose);

	double *pdf1 = new double[length_power_of_two];
	double *pdf2 = new double[length_power_of_two];	
	double *pdf_out = new double[length_power_of_two];	

	pdf1 = in_PDF1.pdf;
	pdf2 = in_PDF2.pdf;

	pdf2 = quantize_message(pdf_awgn_func, 1, delta, min, max);

	//Complex inVec1[];
	//Complex inVec2[];
	//Complex outVec4[];

	Complex *inVec1 = new Complex[length_power_of_two];
	Complex *inVec2 = new Complex[length_power_of_two];
	Complex *outVec4 = new Complex[length_power_of_two];	

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec1[i].rl = pdf2[i];
		else
			inVec1[i].rl = 0;
	}

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec2[i].rl = pdf2[i];
		else
			inVec2[i].rl = 0;
	}

	normal_convolution(inVec1, inVec2, length_power_of_two, outVec4, verbose);

	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			pdf_out[i] = outVec4[i].rl;
		else
			pdf_out[i] = 0;
	}

	delete [] pdf1;
	delete [] pdf2;
	delete [] pdf_out;
	delete [] inVec1;
	delete [] inVec2;
	delete [] outVec4;

	return true;
}

//bool PDF_MATRIX::pdf_convolution(PDF_Element in_PDF1, PDF_Element in_PDF2, PDF_Element out_PDF, bool verbose)
PDF_Element PDF_MATRIX::pdf_convolution(PDF_Element in_PDF1, PDF_Element in_PDF2, bool verbose)
{
	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;

	int  length_power_of_two;

	length_power_of_two = convert_power_of_two(length, verbose);

	double *pdf1 = new double[length];
	double *pdf2 = new double[length];
	double *pdf_out = new double[length];

	PDF_Element out_PDF;
	out_PDF.max = MAX_RANGE;
	out_PDF.min = MIN_RANGE;
	out_PDF.length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	out_PDF.pdf = new double[length];

	if (verbose)
	{
		printf("///// assigning local pdf ///// \n");

		for (int i = 0; i < length; i++)
		{
			printf("///// value of in_PDF1.pdf[%d] is %f ///// \n",i,in_PDF1.pdf[i]);
			//pdf1[i] = in_PDF1.pdf[i];
			//pdf2[i] = in_PDF2.pdf[i];
		}

		printf("///// assign local pdf done ///// \n");
	}

	Complex *inVec1 = new Complex[length];
	Complex *inVec2 = new Complex[length];
	Complex *outVec4 = new Complex[length];	

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			inVec1[i].rl = pdf2[i];
		else
			inVec1[i].rl = 0;
	}

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			inVec2[i].rl = pdf2[i];
		else
			inVec2[i].rl = 0;
	}

	if(verbose)
		printf("///// Begin normal_convolution ///// \n");

	normal_convolution(inVec1, inVec2, length, outVec4, verbose);

	if(verbose)
		printf("///// Copying to pdf_out///// \n");

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			pdf_out[i] = outVec4[i].rl;
		else
			pdf_out[i] = 0;
	}

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			out_PDF.pdf[i] = outVec4[i].rl;
		else
			out_PDF.pdf[i] = 0;
	}

	//if(verbose)
	//printf("///// pdf_out assignment done ///// \n");

	delete [] pdf1;
	delete [] pdf2;
	delete [] pdf_out;

	//delete [] inVec1;
	//delete [] inVec2;
	//delete [] outVec4;

	//return true;
	return out_PDF;
}


bool PDF_MATRIX::pdf_multiply_convolution(PDF_Element const in_PDF1, PDF_Element const in_PDF2, PDF_Element out_PDF, bool verbose)
{
	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;

	int  length_power_of_two;

	length_power_of_two = convert_power_of_two(length, verbose);

	double *pdf1 = new double[length];
	double *pdf2 = new double[length];	
	double *pdf_out = new double[length];	

	for (int i = 0; i < length; i++)
	{
		pdf1[i] = in_PDF1.pdf[i];
		pdf2[i] = in_PDF2.pdf[i];
	}

	Complex *inVec1 = new Complex[length];
	Complex *inVec2 = new Complex[length];
	Complex *outVec4 = new Complex[length];

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			inVec1[i].rl = pdf2[i];
		else
			inVec1[i].rl = 0;
	}

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			inVec2[i].rl = pdf2[i];
		else
			inVec2[i].rl = 0;
	}

	// 这部分后续放在main中初始化查询表
	const int convolution_length = (convolution_range * 2 / QUANTIZE_DELTA + 1);

	double (*vector)[convolution_length] = new double[convolution_length][convolution_length];

	for(int i = 0; i < (convolution_range * 2 / QUANTIZE_DELTA + 1); i++)
	{
		for (int j = 0; j < (convolution_range * 2 / QUANTIZE_DELTA + 1); j++)
		{
			vector[i][j] = quantize_number((2 * atanh(tanh(0.5 * i) * tanh(0.5 * j))), delta) / delta;
		}
	}

	//multiply_convolution(inVec1, inVec2, length_power_of_two, outVec4);
	multiply_convolution(inVec1, inVec2, length, vector, outVec4, verbose);

	for (int i = 0; i < length; i++)
	{
		if(i < length)
			pdf_out[i] = outVec4[i].rl;
		else
			pdf_out[i] = 0;
	}

	delete [] pdf1;
	delete [] pdf2;
	delete [] pdf_out;

	return true;
}


double minimize_array(double *array, int length)
{
	double min_number = *(array);
	double tmp;

	for (int i = 1; i < length; ++i)
	{
		tmp = *(array + sizeof(int) * i);
		if (tmp < min_number)
			min_number = tmp;
	}

	return min_number;
}


bool Judge_llr(PDF_MATRIX in_MATRIX, bool is_AWGN, bool verbose)
{
	PDF_Element in_PDF_AWGN;
	PDF_Element in_PDF_RAYLEIGH;

	PDF_Element* m_lpBuf = in_MATRIX.getPDF();
	int row = in_MATRIX.getRow();
	int col = in_MATRIX.getCol();
	double *pdf_sum = new double[col];

	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;	

	bool above_threshold = false;

	if(verbose)
	{
		printf("/// the value of row is %d /// \n", row);
		printf("/// the value of col is %d /// \n", col);
	}

	PDF_Element pdf_element_tmp;
	PDF_Element pdf_element_iter;
	
	for (int i = 0; i < col; ++i)
	{
		//pdf_element_tmp = *(m_lpBuf + i * sizeof(PDF_Element));
		pdf_element_tmp = m_lpBuf[i];
		
		for (int l = 0; l < row - 1; l++)
		{		
			printf("/// the value of row is %d /// \n", row);
		}

		// Convolution all the elements of the jth column except the ith row
		for (int l = 0; l < row - 1; l++)
		{
			//out_MATRIX.pdf_convolution_AWGN(PDF_Element_AWGN const in_PDF1, PDF_Element_AWGN const PDF_Element_AWGN, PDF_Element out_PDF);
			printf("/// pdf_convolution times: %d /// \n", l);
			if (l != i)
			{
				pdf_element_iter = m_lpBuf[l * col + i];
				//in_MATRIX.pdf_convolution(pdf_element_tmp, pdf_element_iter, pdf_element_tmp, verbose);
				pdf_element_tmp = in_MATRIX.pdf_convolution(pdf_element_tmp, pdf_element_iter, verbose);
			}
		}

		if(is_AWGN)
		{
			//in_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN, pdf_element_tmp, verbose);
			pdf_element_tmp = in_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN, verbose);
		}

		else
		{
			//in_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH, pdf_element_tmp, verbose);
			pdf_element_tmp = in_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH, verbose);
		}

		for (int m = 0; m < length; m++)
		{
			//printf("The %d th element is %f \n", m, pdf_element_tmp.m_lpBuf[m]);
			pdf_sum[i] += pdf_element_tmp.pdf[m];
		}
	}

	printf("******************************************** Judge Done ****************************************** \n");
	printf("*************************************** Printing the pdf_sum ************************************* \n");

	for (int n = 0; n < col; n++)
	{
		printf("***************************** Pdf_sum of the %dth column is %f *************************** \n", n, pdf_sum[n]);
	}

	if (pdf_sum > 0)
		above_threshold = true;
	else
		above_threshold = false;

	printf("*************************************** Printing the end ***************************************** \n");

	delete [] pdf_sum;

	return above_threshold;
}


bool RCA_predict(PDF_MATRIX HB, double SNR, double p, double R[], double snr_R[], double infor, bool verbose)
{
	double s = 2 * pow(10,SNR/10);
	
	int iter_num = NUMBER_OF_INTERATION;
	int mvout_threshold = 100;

	int row = HB.getRow();
	int col = HB.getCol();
	int matrix_Size = row * col;

	double mu0[col];
	double mvout[col];
	double muout[col];

	bool flag = false;
	bool is_AWGN = true;
	bool flag_row = false;
	bool flag_col = false;
	bool push_seccess = false;
	bool clear_seccess = false;

	PDF_MATRIX tmp_MATRIX; 
	tmp_MATRIX = HB;

	printf("*************************************** row and col of tmp_MATRIX are %d and %d ************************************* \n",tmp_MATRIX.getRow(),tmp_MATRIX.getCol());

	PDF_MATRIX tmp_MATRIX2;
	PDF_MATRIX out_MATRIX;

	printf("*************************************** Line 620 Reached ************************************* \n");

	PDF_Element *row_update_Buf;
	PDF_Element *colmumn_update_Buf;

	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;		

	printf("*************************************** The row and col in RCA_predict of HB is %d and %d ***************************************** \n",row, col);

	for(int k = 0; k < col; k++)
	{
		if(k <= p)
		{
			mu0[k] = 0;
			mvout[k] = 0;
		}

		else
		{
			mu0[k] = s;
			mvout[k] = s;
		}
	}

	for(int k = 0; k < row; k++)
	{
		muout[k] = 0;
	}

	printf("******************************* Print the value of pdf in PDF_MATRIX ************************************\n");

	PDF_Element *Buf_tmp = tmp_MATRIX.getPDF();
	double *pdf_tmp = Buf_tmp[0].pdf;

	printf("/// get pdf_tmp OK ///\n");
	if(verbose)
	{
		for(int cnt = 0; cnt < length; cnt++)
		{
			printf("/// pdf in tmp_MATRIX is %f ;///\n", pdf_tmp[cnt]);
		}
	}

	printf("*********************************** Begin to Iterate ************************************\n");

	for(int l = 0; l < iter_num; l++)
	{
		for(int m = 0; m < row; m++)
		{
			printf("///////// Updating Row Times: %d //////////\n", m);
			push_seccess = false;
			flag_row = matrix_row_update(tmp_MATRIX, tmp_MATRIX2, is_AWGN, verbose);
			printf("///////// matrix_row_update 'flag_row' is %d //////////\n", flag_row);
			printf("///////// matrix_Size is %d with row: %d and col: %d //////////\n", matrix_Size, row, col);
			//tmp_MATRIX = tmp_MATRIX2;
			
			// write tmp_MATRIX2 into tmp_MATRIX
			//for(int k = 0; k < matrix_Size; k = k + 1)
			for(int k = 0; k < 2; k = k + 1)
			{
				row_update_Buf = tmp_MATRIX2.getPDF();
				push_seccess = tmp_MATRIX.Push(row_update_Buf[k].pdf, 0, 0, row, col, k, length, verbose);  //
				printf("///////// push is seccess or not : %d //////////\n", push_seccess);
			}//

			push_seccess = false;

			// Clear tmp_MATRIX2
			clear_seccess = tmp_MATRIX2.clear(verbose);		

			printf("///////// clear is seccess or not : %d //////////\n", clear_seccess);
			clear_seccess = false;

			printf("///////// assignment done //////////\n");
			//free tmp_MATRIX2;
			//PDF_MATRIX tmp_MATRIX2;
		}
		
		for(int n = 0; n < col; n++)
		{
			printf("///////// Updating Col Times: %d /////////\n", n);
			push_seccess = false;
			flag_col = matrix_column_update(tmp_MATRIX, tmp_MATRIX2, is_AWGN, verbose);

			printf("///////// matrix_row_update 'flag_col' is %d //////////\n", flag_col);
			printf("///////// matrix_Size is %d with row: %d and col: %d //////////\n", matrix_Size, row, col);
			//tmp_MATRIX = tmp_MATRIX2;
			
			// write tmp_MATRIX2 into tmp_MATRIX
			//for(int k = 0; k < matrix_Size; k = k + 1)
			for(int k = 0; k < 2; k = k + 1)
			{
				row_update_Buf = tmp_MATRIX2.getPDF();
				push_seccess = tmp_MATRIX.Push(row_update_Buf[k].pdf, 0, 0, row, col, k, length, verbose);  //
				printf("///////// push is seccess or not : %d //////////\n", push_seccess);
			}//

			push_seccess = false;

			// Clear tmp_MATRIX2
			clear_seccess = tmp_MATRIX2.clear(verbose);
			clear_seccess = false;

			printf("///////// clear is seccess or not : %d //////////\n", clear_seccess);
			printf("///////// assignment done //////////\n");
		}
	}

	printf("*********************************** Iteration Done ************************************\n");

	flag = Judge_llr(tmp_MATRIX, is_AWGN, verbose);

	if(flag)
		printf(" The LDPC code has passed the threshhold. \n");
	else
		printf(" The LDPC code failed to pass the threshhold. \n");

	return flag;
}


bool matrix_row_update(PDF_MATRIX in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose)
{
	int row = in_MATRIX.getRow();
	int col = in_MATRIX.getCol();

	PDF_Element *local_Buf = in_MATRIX.getPDF();
	PDF_Element in_PDF_AWGN;
	PDF_Element in_PDF_RAYLEIGH;

	PDF_Element pdf_element_tmp;
	PDF_Element pdf_element_iter;

	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;	

	in_PDF_AWGN.max = MAX_RANGE;
	in_PDF_AWGN.min = MIN_RANGE;
	in_PDF_AWGN.length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	in_PDF_AWGN.pdf = quantize_message(pdf_awgn_func, 1, delta, min, max);

	in_PDF_RAYLEIGH.max = MAX_RANGE;
	in_PDF_RAYLEIGH.min = MIN_RANGE;
	in_PDF_RAYLEIGH.length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	in_PDF_RAYLEIGH.pdf = quantize_message(pdf_Rayleigh_func, 1, delta, min, max);

	//for(int i = 0; i < row; i++)
	for(int i = 0; i < 1; i++)
	{
		printf("///// Entering row iteration %d, while row number is %d \n", i, row);

		for(int j = 0; j < col; j++)
		{
			//pdf_element_tmp = *(m_lpBuf + ((i * col + j) * sizeof(PDF_Element)));
			pdf_element_tmp = local_Buf[i * col + j];

			if(verbose)
			{
				for (int k = 0; k < length; k++)
				{
					printf("///// value of pdf_element_tmp.pdf[%d] is %f ///// \n", k ,pdf_element_tmp.pdf[k]);
				}
			}

			for (int l = 0; l < row - 1; l++)
			{
				//out_MATRIX.pdf_convolution_AWGN(PDF_Element_AWGN const in_PDF1, PDF_Element_AWGN const PDF_Element_AWGN, PDF_Element out_PDF);
				//if(verbose)
				printf("///// convolution in row update: %d. \n", l);

				if(l != i)
				{
					pdf_element_iter = local_Buf[l * col + j];
					//out_MATRIX.pdf_convolution(pdf_element_tmp, pdf_element_iter, pdf_element_tmp, verbose);
					pdf_element_tmp = out_MATRIX.pdf_convolution(pdf_element_tmp, pdf_element_iter, verbose);
				}
			}

			//out_MATRIX.pdf_convolution_AWGN(PDF_Element_AWGN const in_PDF1, PDF_Element_AWGN const out_PDF, PDF_Element_AWGN out_PDF);
			if(is_AWGN)
				//out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN , pdf_element_tmp, verbose);
				pdf_element_tmp = out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN, verbose);
			else
				//out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH , pdf_element_tmp, verbose);
				pdf_element_tmp = out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH, verbose);

			//out_MATRIX.m_lpBuf[i * col + j] = pdf_element_tmp;
			bool push_seccess = out_MATRIX.Push(pdf_element_tmp.pdf, i, j, row, col, 0, length, verbose);
		}
	}

	printf("///// before delete \n");

	//delete [] local_Buf;

	printf("///// after delete \n");

	return true;
}


//bool PDF_MATRIX::matrix_column_update(PDF_MATRIX const in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN)
bool matrix_column_update(PDF_MATRIX in_MATRIX, PDF_MATRIX out_MATRIX, bool is_AWGN, bool verbose)
{
	int row = in_MATRIX.getRow();
	int col = in_MATRIX.getCol();
	//for(int i = 0; i < in_MATRIX.Row; i++)
	PDF_Element* local_Buf = in_MATRIX.getPDF();
	PDF_Element in_PDF_AWGN;
	PDF_Element in_PDF_RAYLEIGH;	

	PDF_Element pdf_element_tmp;
	PDF_Element pdf_element_iter;	

	double  min = -DELTA_AWGN;
	double  max = DELTA_AWGN;
	double  delta = QUANTIZE_DELTA;

	int  length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;		

	for(int i = 0; i < row; i++)
	{
		printf(" // Entering column iteration %d, while row number is %d \n", i, row);

		for(int j = 0; j < col; j++)
		{
			//p[k] = sum(pa*pb)		
			pdf_element_tmp = local_Buf[i * col + j];

			if(verbose)
			{
				for (int k = 0; k < length; k++)
				{
					printf("///// value of pdf_element_tmp.pdf[%d] is %f ///// \n",k,pdf_element_tmp.pdf[k]);
				}
			}			

			//mam
			for (int l = 0; l < col - 1; l++)
			{
				//out_MATRIX.pdf_convolution_AWGN(PDF_Element_AWGN const in_PDF1, PDF_Element_AWGN const PDF_Element_AWGN, PDF_Element out_PDF);
				if(l != j)
				{
					pdf_element_iter = local_Buf[l * col + j];
					out_MATRIX.pdf_multiply_convolution(pdf_element_tmp, pdf_element_iter, pdf_element_tmp, verbose);
				}
			}

			if(is_AWGN)
				//out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN , pdf_element_tmp, verbose);
				pdf_element_tmp = out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_AWGN, verbose);
			else
				//out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH, pdf_element_tmp, verbose);
				pdf_element_tmp = out_MATRIX.pdf_convolution(pdf_element_tmp, in_PDF_RAYLEIGH, verbose);

			//bool push_seccess = out_MATRIX.m_lpBuf[i * col + j] = pdf_element_tmp;
			bool push_seccess = out_MATRIX.Push(pdf_element_tmp.pdf, i, j, row, col, 0, length, verbose);
		}
	}

	delete [] local_Buf;
	return true;
}


