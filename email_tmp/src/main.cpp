/*#include "defines.hpp"*/
#include <iostream> //standar library
#include <stdio.h>
#include <string.h>
//#include </usr/local/MATLAB/R2016b/extern/include/mat.h>
#include <mat.h>
#include <Eigen/Dense>
#include <vector>
//#include "El.hpp"

//#include "stdafx.h"  precompile which is only used in Visual Studio
#include "pdf_matrix.hpp"
#include "integral.hpp"

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	bool        verbose = VERBOSE;
	double		score;
//	IplImage	*img;
//	CvMat		temp;
//	Mat			image_draw;
//	Mat			image_temp;
	float 		*ptr;
	float 		*ptr_temp;
	float 		*ptr_thresholded;
	float 		*ptr_scored;
	float 		*ptr_removed_duplicate;
	int         flag_box_exist = 0;

	double		SNR = -11;
	double      p = 2;
	double      R[5001];
	double      snr_R[5001];
	double		infor = 22;

	double 		*HB_TMP; //array for storing data in a matrix

	vector<double> v;
	PDF_Element *m_lpBuf;

	MATFile 	*pmatFile = NULL;
	mxArray 	*pMxArray = NULL;

	// get the length of pmf
	int     	length_power_of_two;
	double  	min = -DELTA_AWGN;
	double  	max = DELTA_AWGN;
	double  	delta = QUANTIZE_DELTA;	
	int     	length = floor(max/delta) - ceil(min/delta) + 1;
	double  	*pdf_tmp = new double[length];

	printf("**************************************** Initialization ************************************ \n");
	length_power_of_two = convert_power_of_two(length, verbose);

	for (int i = 0; i < 5001; ++i)
	{
		R[i] = 0;
		snr_R[i] = 0;
	}

	// Get matrix data from .mat file
	char *file = "HB_64800_10_15.mat";

	printf("***************************************** Opening File ************************************* \n");
	pmatFile = matOpen(file,"r");
	if (pmatFile == NULL)
	{
		printf("Error opening file %s\n", file);
		return(1);
	}

	printf("**************************************** Opening Matrix ************************************ \n");

	// Judge whether the file is empty
	pMxArray = matGetVariable(pmatFile, "B");

	if (pMxArray != NULL && mxIsDouble(pMxArray) && !mxIsEmpty(pMxArray))
	{
		//copy data
		printf("************************************ Opening Success ************************************ \n");
        mwSize num = mxGetNumberOfElements(pMxArray);
        double *pr = mxGetPr(pMxArray);
        if (pr != NULL) {
            v.reserve(num); //is faster than resize :-)
            v.assign(pr, pr+num);
        }
	}
	
	HB_TMP = (double*) mxGetData(pMxArray);

	int M = mxGetM(pMxArray);   //Get the row number
	int N = mxGetN(pMxArray);   //Get the col number

	printf("///////////////// The row number is %d and the column number is %d///////////////// \n",M,N);
	if(verbose)
	{
		printf("///////////////// Verifying the Matrix ///////////////// \n");
		printf("///////////////// The row number is %d and the column number is %d///////////////// \n",M,N);
		printf("HB_TMP \n");
	
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				printf("/// The Matrix of %d row and %d col is %f//// \n",i,j,HB_TMP[M*j+i]); 
			}
		}
	}

	printf("************************************ Copying into PDF_MATRIX ************************************ \n");
	// Copy data and convert them into PDF_MATRIX HB, thus the onrigal matrix is tranformed into a PDF matrix
	PDF_MATRIX HB(M,N);

	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			printf("/// HB_TMP[%d,%d] is%f//// \n",i,j,HB_TMP[M*j+i]);
			if(HB_TMP[M*j+i] == 1)
			{
				for(int k = 0; k < length; k++)
				{
					//if ((k == (length / 2 - 1)) || (k == (length / 2)))
					if (k == (length - 1) / 2)
						pdf_tmp[k] = 1;
					else
						pdf_tmp[k] = 0;

					//printf("/// pdf_tmp[%d]  is  %f//// \n",k, pdf_tmp[k]);
				}			
			}

			else
			{
				for(int k = 0; k < length; k++)
				{
					pdf_tmp[k] = 0;
					//printf("/// pdf_tmp[%d]  is  %f//// \n",k, pdf_tmp[k]);
				}	
			}

			if(verbose)
				printf("/// Pdf_tmp assignment Done //// \n");

			HB.Push(pdf_tmp, i, j, M, N, 0, length, verbose);
		}
	}

	printf("/// Closing matfile //// \n");

	matClose(pmatFile);
	//mxFree(HB_TMP);
	//Free(HB_TMP);

	printf("/// Closing Done //// \n");

	printf("*************************************  Opening File Done ************************************ \n");


	printf("****************************  Preparation Done, Checking PDF Data ************************************ \n");
	// 我在这里想打印出类里的pdf数据，用了两种方法。一种是类里的成员函数，另一中是指针。好像都不行。
	printf("length is %d ",length);

	std::cout <<"element length: "<< HB.get_element(61, 61).length << endl;

	for (int k = (length - 1) / 2; k < length; ++k)
	{
		//std::cout <<"element pdf is: "<< HB.get_element(1, 1).pdf[k] << endl;
		//printf("element pdf %d is %f \n",k,HB.get_element(10, 3).pdf[k]);
	}
	
	HB.getPDF_array(1, 1);
	int matrix_size = M*N;
	PDF_Element *tmp_Element = new PDF_Element[matrix_size];
	tmp_Element = HB.getPDF();
	double *pdf_tmp2 = new double[length];

	printf("***********************************  Checking PDF Data Done  ************************************ \n");	

	// Perform Density Evolution of the matrix
	RCA_predict(HB, SNR, p, R, snr_R, infor, verbose);

	printf("**************************************  Predict Done  ************************************ \n");

    while(1)
    {
    	//123
    }

    return 0;
}
