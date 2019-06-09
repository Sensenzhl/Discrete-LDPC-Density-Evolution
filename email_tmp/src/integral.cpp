#include <iostream>
#include <stdio.h>
#include <string.h>
#include <Eigen/Dense>

#include "math.h"
#include "fft.hpp"
#include "integral.hpp"

using namespace std;

//Matrix<double, Dynamic, Dynamic> C;
typedef struct R_mapping
{
	int index_i;
	int index_j;
	double index_k;
}R_mapping;

R_mapping table[g_nBufferSize];

//R_mapping table[] =
//{
//	{ -20,     -20,  13},  { -20,   -19.9,  13},  { -20,   -19.8,  13},  { -20,   -19.7,  13},  { -20,   -19.6,  13},  
//	{ -20,   -19.5,  13},  { -20,   -19.4,  13},  { -20,   -19.3,  13},  { -20,   -19.2,  13},  { -20,   -19.1,  13},  
//	{ -20,     -19,  13},  { -20,   -18.9,  13},  { -20,   -18.8,  13},  { -20,   -18.7,  13},  { -20,   -18.6,  13},  
//	{ -20,   -18.5,  13},  { -20,   -18.4,  12},  { -20,   -18.3,  12},  { -20,   -18.2,  12},  { -20,   -18.1,  12},  
//	{ -20,     -18,  12},  { -20,   -17.9,  12},  { -20,   -17.8,  12},  { -20,   -17.7,  12},  { -20,   -17.6,  12}, 
//	{ -20,   -17.5,  12},  { -20,   -17.4,  12},  { -20,   -17.3,  12},  { -20,   -17.2,  12},  { -20,   -17.1,  12},  
//	{ -20,     -17,  12},  { -20,   -16.9,  12},  { -20,   -16.8,  12},  { -20,   -16.7,  12},  { -20,   -16.6,  11}, 
//    { -20,   -16.5,  11},  { -20,   -16.4,  11},  { -20,   -16.3,  11},  { -20,   -16.2,  11},  { -20,   -16.1,  11},  
//    { -20,     -16,  11},  { -20,   -15.9,  11},  { -20,   -15.8,  11},  { -20,   -15.7,  11},  { -20,   -15.6,  11},  
//    { -20,   -15.5,  11},  { -20,   -15.4,  11},  { -20,   -15.3,  11},  { -20,   -15.2,  11},  { -20,   -15.1,  11},  
//    { -20,     -15,  11},  { -20,   -14.9,  10},  { -20,   -14.8,  10},  { -20,   -14.7,  10},  { -20,   -14.6,  10},  
//    { -20,   -14.5,  10},  { -20,   -14.4,  10},  { -20,   -14.3,  10},  { -20,   -14.2,  10},  { -20,   -14.1,  10},  
//    { -20,     -14,  10},  { -20,   -13.9,  10},  { -20,   -13.8,  10},  { -20,   -13.7,  10},  { -20,   -13.6,  10}, 
//    { -20,   -13.5,  10},  { -20,   -13.4,  10},  { -20,   -13.3,  10},  { -20,   -13.2,   9},  { -20,   -13.1,   9}, 
//    { -20,     -13,   9},  { -20,   -12.9,   9},  { -20,   -12.8,   9},  { -20,   -12.7,   9},  { -20,   -12.6,   9},  
//    { -20,   -12.5,   9},  { -20,   -12.4,   9},  { -20,   -12.3,   9},  { -20,   -12.2,   9},  { -20,   -12.1,   9},  
//    { -20,     -12,   9},  { -20,   -11.9,   9},  { -20,   -11.8,   9},  { -20,   -11.7,   8},  { -20,   -11.6,   8},  
//    { -20,   -11.5,   8},  { -20,   -11.4,   8},  { -20,   -11.3,   8},  { -20,   -11.2,   8},  { -20,   -11.1,   8},  
//    { -20,     -11,   8},  { -20,   -10.9,   8},  { -20,   -10.8,   8},  { -20,   -10.7,   8},  { -20,   -10.6,   8},  
//    { -20,   -10.5,   8},  { -20,   -10.4,   8},  { -20,   -10.3,   8},  { -20,   -10.2,   7},  { -20,   -10.1,   7},  
//    { -20,     -10,   7},  { -20,    -9.9,   7},  { -20,    -9.8,   7},  { -20,    -9.7,   7},  { -20,    -9.6,   7}, 
//    { -20,    -9.5,   7},  { -20,    -9.4,   7},  { -20,    -9.3,   7},  { -20,    -9.2,   7},  { -20,    -9.1,   7},  
//    { -20,      -9,   7},  { -20,    -8.9,   7},  { -20,    -8.8,   7},  { -20,    -8.7,   6},  { -20,    -8.6,   6},  
//    { -20,    -8.5,   6},  { -20,    -8.4,   6},  { -20,    -8.3,   6},  { -20,    -8.2,   6},  { -20,    -8.1,   6},  
//    { -20,      -8,   6},  { -20,    -7.9,   6},  { -20,    -7.8,   6},  { -20,    -7.7,   6},  { -20,    -7.6,   6},  
//    { -20,    -7.5,   6},  { -20,    -7.4,   6},  { -20,    -7.3,   5},  { -20,    -7.2,   5},  { -20,    -7.1,   5},  
//    { -20,      -7,   5},  { -20,    -6.9,   5},  { -20,    -6.8,   5},  { -20,    -6.7,   5},  { -20,    -6.6,   5},  
//    { -20,    -6.5,   5},  { -20,    -6.4,   5},  { -20,    -6.3,   5},  { -20,    -6.2,   5},  { -20,    -6.1,   5},     
//    { -20,      -6,   5},  { -20,    -5.9,   4},  { -20,    -5.8,   4},  { -20,    -5.7,   4},  { -20,    -5.6,   4},  
//    { -20,    -5.5,   4},  { -20,    -5.4,   4},  { -20,    -5.3,   4},  { -20,    -5.2,   4},  { -20,    -5.1,   4},  
//    { -20,      -5,   4},  { -20,    -4.9,   4},  { -20,    -4.8,   4},  { -20,    -4.7,   4},  { -20,    -4.6,   3},  
//    { -20,    -4.5,   3},  { -20,    -4.4,   3},  { -20,    -4.3,   3},  { -20,    -4.2,   3},  { -20,    -4.1,   3},  
//    { -20,      -4,   3},  { -20,    -3.9,   3},  { -20,    -3.8,   3},  { -20,    -3.7,   3},  { -20,    -3.6,   3},  
//    { -20,    -3.5,   3},  { -20,    -3.4,   3},  { -20,    -3.3,   3},  { -20,    -3.2,   2},  { -20,    -3.1,   2},  
//    { -20,      -3,   2},  { -20,    -2.9,   2},  { -20,    -2.8,   2},  { -20,    -2.7,   2},  { -20,    -2.6,   2}, 
//    { -20,    -2.5,   2},  { -20,    -2.4,   2},  { -20,    -2.3,   2},  { -20,    -2.2,   2},  { -20,    -2.1,   2}, 
//    { -20,      -2,   2},  { -20,    -1.9,   1},  { -20,    -1.8,   1},  { -20,    -1.7,   1},  { -20,    -1.6,   1},  
//    { -20,    -1.5,   1},  { -20,    -1.4,   1},  { -20,    -1.3,   1},  { -20,    -1.2,   1},  { -20,    -1.1,   1},  
//    { -20,      -1,   1},  { -20,    -0.9,   1},  { -20,    -0.8,   1},  { -20,    -0.7,   1},  { -20,    -0.6,   0},  
//    { -20,    -0.5,   0},  { -20,    -0.4,   0},  { -20,    -0.3,   0},  { -20,    -0.2,   0},  { -20,    -0.1,   0},  
//    { -20,       0,   0},  { -20,     0.1,   0},  { -20,     0.2,   0},  { -20,     0.3,   0},  { -20,     0.4,   0},
//    { -20,     0.5,   0},  { -20,     0.6,   0},  { -20,     0.7,  -1},  { -20,     0.8,  -1},  { -20,     0.9,  -1},
//    { -20,       1,  -1},  { -20,     1.1,  -1},  { -20,     1.2,  -1},  { -20,     1.3,  -1},  { -20,     1.4,  -1},
//    { -20,     1.5,  -1},  { -20,     1.6,  -1},  { -20,     1.7,  -1},  { -20,     1.8,  -1},  { -20,     1.9,  -1},
//    { -20,       2,  -1},  { -20,     1.1,  -1},  { -20,     1.2,  -1},  { -20,     1.3,  -1},  { -20,     1.4,  -1},
//    { -20,     1.5,  -1},  { -20,     1.6,  -1},  { -20,     1.7,  -1},  { -20,     1.8,  -1},  { -20,     1.9,  -1},
//    { -20,       2,  -2},  { -20,     2.1,  -2},  { -20,     2.2,  -2},  { -20,     2.3,  -2},  { -20,     2.4,  -2},
//    { -20,     2.5,  -2},  { -20,     2.6,  -2},  { -20,     2.7,  -2},  { -20,     2.8,  -2},  { -20,     2.9,  -2},
//    { -20,       3,  -2},  { -20,     3.1,  -2},  { -20,     3.2,  -2},  { -20,     3.3,  -3},  { -20,     3.4,  -3},
//    { -20,     3.5,  -3},  { -20,     3.6,  -3},  { -20,     3.7,  -3},  { -20,     3.8,  -3},  { -20,     3.9,  -3},
//    { -20,       4,  -3},  { -20,     4.1,  -3},  { -20,     4.2,  -3},  { -20,     4.3,  -3},  { -20,     4.4,  -3},
//    { -20,     4.5,  -3},  { -20,     4.6,  -3},  { -20,     4.7,  -4},  { -20,     4.8,  -4},  { -20,     4.9,  -4},
//    { -20,       5,  -4},  { -20,     5.1,  -4},  { -20,     5.2,  -4},  { -20,     5.3,  -4},  { -20,     5.4,  -4},
//    { -20,     5.5,  -4},  { -20,     5.6,  -4},  { -20,     5.7,  -4},  { -20,     5.8,  -4},  { -20,     5.9,  -4},
//    { -20,       6,  -5},  { -20,     6.1,  -5},  { -20,     6.2,  -5},  { -20,     6.3,  -5},  { -20,     6.4,  -5},
//    { -20,     6.5,  -5},  { -20,     6.6,  -5},  { -20,     6.7,  -5},  { -20,     6.8,  -5},  { -20,     6.9,  -5},
//    { -20,       7,  -5},  { -20,     7.1,  -5},  { -20,     7.2,  -5},  { -20,     7.3,  -5},  { -20,     7.4,  -6},
//    { -20,     7.5,  -6},  { -20,     7.6,  -6},  { -20,     7.7,  -6},  { -20,     7.8,  -6},  { -20,     7.9,  -6},
//    { -20,       8,  -6},  { -20,     8.1,  -6},  { -20,     8.2,  -6},  { -20,     8.3,  -6},  { -20,     8.4,  -6},
//    { -20,     8.5,  -6},  { -20,     8.6,  -6},  { -20,     8.7,  -6},  { -20,     8.8,  -7},  { -20,     8.9,  -7},
//    { -20,       9,  -7},  { -20,     9.1,  -7},  { -20,     9.2,  -7},  { -20,     9.3,  -7},  { -20,     9.4,  -7},
//    { -20,     9.5,  -7},  { -20,     9.6,  -7},  { -20,     9.7,  -7},  { -20,     9.8,  -7},  { -20,     9.9,  -7},
//    { -20,      10,  -7},  { -20,    10.1,  -7},  { -20,    10.2,  -7},  { -20,    10.3,  -8},  { -20,    10.4,  -8},
//    { -20,    10.5,  -8},  { -20,    10.6,  -8},  { -20,    10.7,  -8},  { -20,    10.8,  -8},  { -20,    10.9,  -8},
//    { -20,      11,  -8},  { -20,    11.1,  -8},  { -20,    11.2,  -8},  { -20,    11.3,  -8},  { -20,    11.4,  -8},
//    { -20,    11.5,  -8},  { -20,    11.6,  -8},  { -20,    11.7,  -8},  { -20,    11.8,  -9},  { -20,    11.9,  -9},
//    { -20,      12,  -9},  { -20,    12.1,  -9},  { -20,    12.2,  -9},  { -20,    12.3,  -9},  { -20,    12.4,  -9},
//    { -20,    12.5,  -9},  { -20,    12.6,  -9},  { -20,    12.7,  -9},  { -20,    12.8,  -9},  { -20,    12.9,  -9},
//    { -20,      13,  -9},  { -20,    13.1,  -9},  { -20,    13.2,  -9},  { -20,    13.3, -10},  { -20,    13.4, -10},
//    { -20,    13.5, -10},  { -20,    13.6, -10},  { -20,    13.7, -10},  { -20,    13.8, -10},  { -20,    13.9, -10},
//    { -20,      14, -10},  { -20,    14.1, -10},  { -20,    14.2, -10},  { -20,    14.3, -10},  { -20,    14.4, -10},
//    { -20,    14.5, -10},  { -20,    14.6, -10},  { -20,    14.7, -10},  { -20,    14.8, -10},  { -20,    14.9, -10},
//    { -20,      15, -11},  { -20,    15.1, -11},  { -20,    15.2, -11},  { -20,    15.3, -11},  { -20,    15.4, -11},
//    { -20,    15.5, -11},  { -20,    15.6, -11},  { -20,    15.7, -11},  { -20,    15.8, -11},  { -20,    15.9, -11},
//    { -20,      16, -11},  { -20,    16.1, -11},  { -20,    16.2, -11},  { -20,    16.3, -11},  { -20,    16.4, -11},
//    { -20,    16.5, -11},  { -20,    16.6, -11},  { -20,    16.7, -12},  { -20,    16.8, -12},  { -20,    16.9, -12},
//    { -20,      17, -12},  { -20,    17.1, -12},  { -20,    17.2, -12},  { -20,    17.3, -12},  { -20,    17.4, -12},
//    { -20,    17.5, -12},  { -20,    17.6, -12},  { -20,    17.7, -12},  { -20,    17.8, -12},  { -20,    17.9, -12},
//    { -20,      18, -12},  { -20,    18.1, -12},  { -20,    18.2, -12},  { -20,    18.3, -12},  { -20,    18.4, -12},
//    { -20,    18.4, -12},  { -20,    18.5, -13},  { -20,    18.6, -13},  { -20,    18.7, -13},  { -20,    18.8, -13},
//    { -20,    18.9, -13},  { -20,      19, -13},  { -20,    19.1, -13},  { -20,    19.2, -13},  { -20,    19.3, -13},  
//    { -20,    19.4, -13},  { -20,    19.5, -13},  { -20,    19.6, -13},  { -20,    19.7, -13},  { -20,    19.8, -13},  
//    { -20,    19.9, -13},  { -20,    19.9, -13},  { -20,      20, -13},
//
//	{ -19.9,     -20,  13},  { -19.9,   -19.9,  13},  { -19.9,   -19.8,  13},  { -19.9,   -19.7,  13},  { -19.9,   -19.6,  13},  
//	{ -19.9,   -19.5,  13},  { -19.9,   -19.4,  13},  { -19.9,   -19.3,  13},  { -19.9,   -19.2,  13},  { -19.9,   -19.1,  13},  
//	{ -19.9,     -19,  13},  { -19.9,   -18.9,  13},  { -19.9,   -18.8,  13},  { -19.9,   -18.7,  13},  { -19.9,   -18.6,  13},  
//	{ -19.9,   -18.5,  12},  { -19.9,   -18.4,  12},  { -19.9,   -18.3,  12},  { -19.9,   -18.2,  12},  { -19.9,   -18.1,  12},  
//	{ -19.9,     -18,  12},  { -19.9,   -17.9,  12},  { -19.9,   -17.8,  12},  { -19.9,   -17.7,  12},  { -19.9,   -17.6,  12}, 
//	{ -19.9,   -17.5,  12},  { -19.9,   -17.4,  12},  { -19.9,   -17.3,  12},  { -19.9,   -17.2,  12},  { -19.9,   -17.1,  12},  
//	{ -19.9,     -17,  12},  { -19.9,   -16.9,  12},  { -19.9,   -16.8,  12},  { -19.9,   -16.7,  11},  { -19.9,   -16.6,  11}, 
//    { -19.9,   -16.5,  11},  { -19.9,   -16.4,  11},  { -19.9,   -16.3,  11},  { -19.9,   -16.2,  11},  { -19.9,   -16.1,  11},  
//    { -19.9,     -16,  11},  { -19.9,   -15.9,  11},  { -19.9,   -15.8,  11},  { -19.9,   -15.7,  11},  { -19.9,   -15.6,  11},  
//    { -19.9,   -15.5,  11},  { -19.9,   -15.4,  11},  { -19.9,   -15.3,  11},  { -19.9,   -15.2,  11},  { -19.9,   -15.1,  11},  
//    { -19.9,     -15,  11},  { -19.9,   -14.9,  10},  { -19.9,   -14.8,  10},  { -19.9,   -14.7,  10},  { -19.9,   -14.6,  10},  
//    { -19.9,   -14.5,  10},  { -19.9,   -14.4,  10},  { -19.9,   -14.3,  10},  { -19.9,   -14.2,  10},  { -19.9,   -14.1,  10},  
//    { -19.9,     -14,  10},  { -19.9,   -13.9,  10},  { -19.9,   -13.8,  10},  { -19.9,   -13.7,  10},  { -19.9,   -13.6,  10}, 
//    { -19.9,   -13.5,  10},  { -19.9,   -13.4,  10},  { -19.9,   -13.3,   9},  { -19.9,   -13.2,   9},  { -19.9,   -13.1,   9}, 
//    { -19.9,     -13,   9},  { -19.9,   -12.9,   9},  { -19.9,   -12.8,   9},  { -19.9,   -12.7,   9},  { -19.9,   -12.6,   9},  
//    { -19.9,   -12.5,   9},  { -19.9,   -12.4,   9},  { -19.9,   -12.3,   9},  { -19.9,   -12.2,   9},  { -19.9,   -12.1,   9},  
//    { -19.9,     -12,   9},  { -19.9,   -11.9,   9},  { -19.9,   -11.8,   9},  { -19.9,   -11.7,   8},  { -19.9,   -11.6,   8},  
//    { -19.9,   -11.5,   8},  { -19.9,   -11.4,   8},  { -19.9,   -11.3,   8},  { -19.9,   -11.2,   8},  { -19.9,   -11.1,   8},  
//    { -19.9,     -11,   8},  { -19.9,   -10.9,   8},  { -19.9,   -10.8,   8},  { -19.9,   -10.7,   8},  { -19.9,   -10.6,   8},  
//    { -19.9,   -10.5,   8},  { -19.9,   -10.4,   8},  { -19.9,   -10.3,   8},  { -19.9,   -10.2,   7},  { -19.9,   -10.1,   7},  
//    { -19.9,     -10,   7},  { -19.9,    -9.9,   7},  { -19.9,    -9.8,   7},  { -19.9,    -9.7,   7},  { -19.9,    -9.6,   7}, 
//    { -19.9,    -9.5,   7},  { -19.9,    -9.4,   7},  { -19.9,    -9.3,   7},  { -19.9,    -9.2,   7},  { -19.9,    -9.1,   7},  
//    { -19.9,      -9,   7},  { -19.9,    -8.9,   7},  { -19.9,    -8.8,   7},  { -19.9,    -8.7,   6},  { -19.9,    -8.6,   6},  
//    { -19.9,    -8.5,   6},  { -19.9,    -8.4,   6},  { -19.9,    -8.3,   6},  { -19.9,    -8.2,   6},  { -19.9,    -8.1,   6},  
//    { -19.9,      -8,   6},  { -19.9,    -7.9,   6},  { -19.9,    -7.8,   6},  { -19.9,    -7.7,   6},  { -19.9,    -7.6,   6},  
//    { -19.9,    -7.5,   6},  { -19.9,    -7.4,   6},  { -19.9,    -7.3,   5},  { -19.9,    -7.2,   5},  { -19.9,    -7.1,   5},  
//    { -19.9,      -7,   5},  { -19.9,    -6.9,   5},  { -19.9,    -6.8,   5},  { -19.9,    -6.7,   5},  { -19.9,    -6.6,   5},  
//    { -19.9,    -6.5,   5},  { -19.9,    -6.4,   5},  { -19.9,    -6.3,   5},  { -19.9,    -6.2,   5},  { -19.9,    -6.1,   5},     
//    { -19.9,      -6,   4},  { -19.9,    -5.9,   4},  { -19.9,    -5.8,   4},  { -19.9,    -5.7,   4},  { -19.9,    -5.6,   4},  
//    { -19.9,    -5.5,   4},  { -19.9,    -5.4,   4},  { -19.9,    -5.3,   4},  { -19.9,    -5.2,   4},  { -19.9,    -5.1,   4},  
//    { -19.9,      -5,   4},  { -19.9,    -4.9,   4},  { -19.9,    -4.8,   4},  { -19.9,    -4.7,   4},  { -19.9,    -4.6,   3},  
//    { -19.9,    -4.5,   3},  { -19.9,    -4.4,   3},  { -19.9,    -4.3,   3},  { -19.9,    -4.2,   3},  { -19.9,    -4.1,   3},  
//    { -19.9,      -4,   3},  { -19.9,    -3.9,   3},  { -19.9,    -3.8,   3},  { -19.9,    -3.7,   3},  { -19.9,    -3.6,   3},  
//    { -19.9,    -3.5,   3},  { -19.9,    -3.4,   3},  { -19.9,    -3.3,   2},  { -19.9,    -3.2,   2},  { -19.9,    -3.1,   2},  
//    { -19.9,      -3,   2},  { -19.9,    -2.9,   2},  { -19.9,    -2.8,   2},  { -19.9,    -2.7,   2},  { -19.9,    -2.6,   2}, 
//    { -19.9,    -2.5,   2},  { -19.9,    -2.4,   2},  { -19.9,    -2.3,   2},  { -19.9,    -2.2,   2},  { -19.9,    -2.1,   2}, 
//    { -19.9,      -2,   2},  { -19.9,    -1.9,   1},  { -19.9,    -1.8,   1},  { -19.9,    -1.7,   1},  { -19.9,    -1.6,   1},  
//    { -19.9,    -1.5,   1},  { -19.9,    -1.4,   1},  { -19.9,    -1.3,   1},  { -19.9,    -1.2,   1},  { -19.9,    -1.1,   1},  
//    { -19.9,      -1,   1},  { -19.9,    -0.9,   1},  { -19.9,    -0.8,   1},  { -19.9,    -0.7,   1},  { -19.9,    -0.6,   0},  
//    { -19.9,    -0.5,   0},  { -19.9,    -0.4,   0},  { -19.9,    -0.3,   0},  { -19.9,    -0.2,   0},  { -19.9,    -0.1,   0},  
//    { -19.9,       0,   0},  { -19.9,     0.1,   0},  { -19.9,     0.2,   0},  { -19.9,     0.3,   0},  { -19.9,     0.4,   0},
//    { -19.9,     0.5,   0},  { -19.9,     0.6,   0},  { -19.9,     0.7,  -1},  { -19.9,     0.8,  -1},  { -19.9,     0.9,  -1},
//    { -19.9,       1,  -1},  { -19.9,     1.1,  -1},  { -19.9,     1.2,  -1},  { -19.9,     1.3,  -1},  { -19.9,     1.4,  -1},
//    { -19.9,     1.5,  -1},  { -19.9,     1.6,  -1},  { -19.9,     1.7,  -1},  { -19.9,     1.8,  -1},  { -19.9,     1.9,  -1},
//    { -19.9,       2,  -1},  { -19.9,     1.1,  -1},  { -19.9,     1.2,  -1},  { -19.9,     1.3,  -1},  { -19.9,     1.4,  -1},
//    { -19.9,     1.5,  -1},  { -19.9,     1.6,  -1},  { -19.9,     1.7,  -1},  { -19.9,     1.8,  -1},  { -19.9,     1.9,  -1},
//    { -19.9,       2,  -2},  { -19.9,     2.1,  -2},  { -19.9,     2.2,  -2},  { -19.9,     2.3,  -2},  { -19.9,     2.4,  -2},
//    { -19.9,     2.5,  -2},  { -19.9,     2.6,  -2},  { -19.9,     2.7,  -2},  { -19.9,     2.8,  -2},  { -19.9,     2.9,  -2},
//    { -19.9,       3,  -2},  { -19.9,     3.1,  -2},  { -19.9,     3.2,  -2},  { -19.9,     3.3,  -2},  { -19.9,     3.4,  -3},
//    { -19.9,     3.5,  -3},  { -19.9,     3.6,  -3},  { -19.9,     3.7,  -3},  { -19.9,     3.8,  -3},  { -19.9,     3.9,  -3},
//    { -19.9,       4,  -3},  { -19.9,     4.1,  -3},  { -19.9,     4.2,  -3},  { -19.9,     4.3,  -3},  { -19.9,     4.4,  -3},
//    { -19.9,     4.5,  -3},  { -19.9,     4.6,  -3},  { -19.9,     4.7,  -4},  { -19.9,     4.8,  -4},  { -19.9,     4.9,  -4},
//    { -19.9,       5,  -4},  { -19.9,     5.1,  -4},  { -19.9,     5.2,  -4},  { -19.9,     5.3,  -4},  { -19.9,     5.4,  -4},
//    { -19.9,     5.5,  -4},  { -19.9,     5.6,  -4},  { -19.9,     5.7,  -4},  { -19.9,     5.8,  -4},  { -19.9,     5.9,  -4},
//    { -19.9,       6,  -4},  { -19.9,     6.1,  -5},  { -19.9,     6.2,  -5},  { -19.9,     6.3,  -5},  { -19.9,     6.4,  -5},
//    { -19.9,     6.5,  -5},  { -19.9,     6.6,  -5},  { -19.9,     6.7,  -5},  { -19.9,     6.8,  -5},  { -19.9,     6.9,  -5},
//    { -19.9,       7,  -5},  { -19.9,     7.1,  -5},  { -19.9,     7.2,  -5},  { -19.9,     7.3,  -5},  { -19.9,     7.4,  -6},
//    { -19.9,     7.5,  -6},  { -19.9,     7.6,  -6},  { -19.9,     7.7,  -6},  { -19.9,     7.8,  -6},  { -19.9,     7.9,  -6},
//    { -19.9,       8,  -6},  { -19.9,     8.1,  -6},  { -19.9,     8.2,  -6},  { -19.9,     8.3,  -6},  { -19.9,     8.4,  -6},
//    { -19.9,     8.5,  -6},  { -19.9,     8.6,  -6},  { -19.9,     8.7,  -6},  { -19.9,     8.8,  -7},  { -19.9,     8.9,  -7},
//    { -19.9,       9,  -7},  { -19.9,     9.1,  -7},  { -19.9,     9.2,  -7},  { -19.9,     9.3,  -7},  { -19.9,     9.4,  -7},
//    { -19.9,     9.5,  -7},  { -19.9,     9.6,  -7},  { -19.9,     9.7,  -7},  { -19.9,     9.8,  -7},  { -19.9,     9.9,  -7},
//    { -19.9,      10,  -7},  { -19.9,    10.1,  -7},  { -19.9,    10.2,  -7},  { -19.9,    10.3,  -8},  { -19.9,    10.4,  -8},
//    { -19.9,    10.5,  -8},  { -19.9,    10.6,  -8},  { -19.9,    10.7,  -8},  { -19.9,    10.8,  -8},  { -19.9,    10.9,  -8},
//    { -19.9,      11,  -8},  { -19.9,    11.1,  -8},  { -19.9,    11.2,  -8},  { -19.9,    11.3,  -8},  { -19.9,    11.4,  -8},
//    { -19.9,    11.5,  -8},  { -19.9,    11.6,  -8},  { -19.9,    11.7,  -8},  { -19.9,    11.8,  -9},  { -19.9,    11.9,  -9},
//    { -19.9,      12,  -9},  { -19.9,    12.1,  -9},  { -19.9,    12.2,  -9},  { -19.9,    12.3,  -9},  { -19.9,    12.4,  -9},
//    { -19.9,    12.5,  -9},  { -19.9,    12.6,  -9},  { -19.9,    12.7,  -9},  { -19.9,    12.8,  -9},  { -19.9,    12.9,  -9},
//    { -19.9,      13,  -9},  { -19.9,    13.1,  -9},  { -19.9,    13.2,  -9},  { -19.9,    13.3,  -9},  { -19.9,    13.4, -10},
//    { -19.9,    13.5, -10},  { -19.9,    13.6, -10},  { -19.9,    13.7, -10},  { -19.9,    13.8, -10},  { -19.9,    13.9, -10},
//    { -19.9,      14, -10},  { -19.9,    14.1, -10},  { -19.9,    14.2, -10},  { -19.9,    14.3, -10},  { -19.9,    14.4, -10},
//    { -19.9,    14.5, -10},  { -19.9,    14.6, -10},  { -19.9,    14.7, -10},  { -19.9,    14.8, -10},  { -19.9,    14.9, -10},
//    { -19.9,      15, -11},  { -19.9,    15.1, -11},  { -19.9,    15.2, -11},  { -19.9,    15.3, -11},  { -19.9,    15.4, -11},
//    { -19.9,    15.5, -11},  { -19.9,    15.6, -11},  { -19.9,    15.7, -11},  { -19.9,    15.8, -11},  { -19.9,    15.9, -11},
//    { -19.9,      16, -11},  { -19.9,    16.1, -11},  { -19.9,    16.2, -11},  { -19.9,    16.3, -11},  { -19.9,    16.4, -11},
//    { -19.9,    16.5, -11},  { -19.9,    16.6, -11},  { -19.9,    16.7, -11},  { -19.9,    16.8, -12},  { -19.9,    16.9, -12},
//    { -19.9,      17, -12},  { -19.9,    17.1, -12},  { -19.9,    17.2, -12},  { -19.9,    17.3, -12},  { -19.9,    17.4, -12},
//    { -19.9,    17.5, -12},  { -19.9,    17.6, -12},  { -19.9,    17.7, -12},  { -19.9,    17.8, -12},  { -19.9,    17.9, -12},
//    { -19.9,      18, -12},  { -19.9,    18.1, -12},  { -19.9,    18.2, -12},  { -19.9,    18.3, -12},  { -19.9,    18.4, -12},
//    { -19.9,    18.4, -12},  { -19.9,    18.5, -12},  { -19.9,    18.6, -13},  { -19.9,    18.7, -13},  { -19.9,    18.8, -13},
//    { -19.9,    18.9, -13},  { -19.9,      19, -13},  { -19.9,    19.1, -13},  { -19.9,    19.2, -13},  { -19.9,    19.3, -13},  
//    { -19.9,    19.4, -13},  { -19.9,    19.5, -13},  { -19.9,    19.6, -13},  { -19.9,    19.7, -13},  { -19.9,    19.8, -13},  
//    { -19.9,    19.9, -13},  { -19.9,    19.9, -13},  { -19.9,      20, -13},


//vector<vector<int>> array;
double quantize_number(double w, double delta)
{
	double Q_value;
	double Q_value_quantized;
	//double w;

	if(w >= delta / 2)
		Q_value = delta * floor(w/delta + 1/2);
	else if(w <= -delta / 2)
		Q_value = delta * ceil(w/delta - 1/2);
	else
		Q_value = 0;

	Q_value_quantized = Q_value;	

	return Q_value_quantized;
}

double (*vector_initialization(double delta))[convolution_length]{

	const int convolution_length = (convolution_range * 2 / QUANTIZE_DELTA + 1);

	double (*vector)[convolution_length] = new double[convolution_length][convolution_length];

	//vector[convolution_length][convolution_length];
	
	for(int i = 0; i < (convolution_range * 2 / QUANTIZE_DELTA + 1); i++)
	{
		for (int j = 0; j < (convolution_range * 2 / QUANTIZE_DELTA + 1); j++)
		{
			vector[i][j] = quantize_number((2 * atanh(tanh(0.5 * i) * tanh(0.5 * j))), delta) / delta;
		}
	}

	return vector;
}



double func1(double x){
	return x;
}

double func2(double x){
	return x*x;
}

double pdf_awgn_func(double x, double variance){

	//double pi = 4 * atan(1.0); 

	//printf("parameters are %f %f  %f  %f \n", variance, M_PI, pow(variance,2), (2 / pow(variance,2)));

	return sqrt(pow(variance,2) / (8 * M_PI)) * exp(- pow((x - (2 / pow(variance,2))),2) * pow(variance,2) / 8.0f);
	//return sqrt(1 / (8 * 3.14)) * exp(- pow((x - 2),2) / 8.0f);
}

double pdf_Rayleigh_func(double x, double variance){

	//double c,d,e;

	//printf("the pdf_Rayleigh_func result of x is %f\n", pow(variance,2) * exp(d) / (2 * sqrt(2 * pow(variance,2) + 1)));
	return pow(variance,2) * exp(-1 / 2.0f * (sqrt(2 * pow(variance,2) + 1) * abs(x) - x)) / (2 * sqrt(2 * pow(variance,2) + 1));
	//return pow(variance,2) * exp(d) / (2 * sqrt(2 * pow(variance,2) + 1));
}

//Traditional integral using rectangle method
double tanh_rule(double a, double b)
{
	double result;

	result = 2 * atanh(tanh(a / 2.0f) * tanh(b / 2.0f));

	return result;
}


double integral(double(*f)(double, double), double min, double max, double variance)
{

	double result = 0;
	const int N = INTEGRAL_ITERATION_SIZE; //define integral step

	double delta = (max - min) / N;

	for (double i = min+delta; i < max; i+=delta)
	{
		result += f(i,variance)*delta;
	}
	return result;
}


double *quantize_message(double(*f)(double, double), double variance, double delta, double min, double max)
{
	double Q_value;
	int    length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double *Q_value_quantized = new double[length];
	int    j = 0;
	double i,w;

	for (double k = ceil(min/QUANTIZE_DELTA); k * QUANTIZE_DELTA < max; k++)
	{
		i = k * QUANTIZE_DELTA;
		w = f(i, variance);

		if(w >= delta / 2)
			Q_value = delta * floor(w/delta + 1/2);
		else if(w <= -delta / 2)
			Q_value = delta * ceil(w/delta - 1/2);
		else
			Q_value = 0;

		Q_value_quantized[j] = Q_value;		
		j++;
	}

	return Q_value_quantized;
}

int convert_power_of_two(int length, bool verbose)
{
	int length_power_of_two = 1;
	double k = length;

	while(k > 1)
	{
		k = k / 2;
		length_power_of_two = length_power_of_two * 2;
	}

	if(verbose)
		printf("length_power_of_two is %d \n", length_power_of_two);

	return length_power_of_two;
}


double operator_r(double a, double b, double delta)
{
	double r;
	double tmp;

	tmp = 2 * atanh(tanh(0.5 * a) * tanh(0.5 * b));
	r = quantize_number(tmp,delta);
 
	return r;
}


bool normal_convolution(Complex inVec1[], Complex inVec2[], int vecLen, Complex outVec4[], bool verbose)
{
	Fft    t,t2,t3;
	int    length_power_of_two;

	length_power_of_two = convert_power_of_two(vecLen, verbose);

	Complex *outVec   = new Complex[length_power_of_two];
	Complex *outVec2  = new Complex[length_power_of_two];
	Complex *outVec3  = new Complex[length_power_of_two];
	//Complex *outVec4  = new Complex[length_power_of_two];

	if(verbose)
		printf("///// Before FFT ///// \n");

	t.fft(inVec1, vecLen, outVec);
	t2.fft(inVec2, vecLen, outVec2);
	t3.multiply(outVec,outVec2,length_power_of_two, outVec3);

	t.ifft(outVec, length_power_of_two, outVec4);

	if(verbose)
		printf("///// FFT Finished///// \n");

	return true;
}


bool multiply_convolution(Complex const inVec1[], Complex const inVec2[], int const vecLen, double vector[convolution_length][convolution_length], Complex outVec4[], bool verbose)
{
	Fft     t,t2,t3;
	int     length_power_of_two;

	length_power_of_two = convert_power_of_two(vecLen, verbose);

	double  pdf1[length_power_of_two];
	double  pdf2[length_power_of_two];
	double  pdf_out[length_power_of_two];
	double  sum[length_power_of_two];
	double  tmp = 0;
	double  k = 0;

	Complex *outVec   = new Complex[length_power_of_two];
	Complex *outVec2  = new Complex[length_power_of_two];
	Complex *outVec3  = new Complex[length_power_of_two];
	//Complex *outVec4  = new Complex[length_power_of_two];

	//t.fft(inVec, vecLen, outVec);
	//t2.fft(inVec2, vecLen, outVec2);
	//t3.multiply(outVec,outVec2,length_power_of_two, outVec3);
	//t.ifft(outVec, length_power_of_two, outVec4);

	for(int i = 0; i < length_power_of_two; i++)
	{
		pdf1[i] = inVec1[i].rl;
		pdf2[i] = inVec2[i].rl;
		sum[i] = 0;
	}

	for(int i = 0; i < length_power_of_two; i++)
	{
		for(int j = 0; j < length_power_of_two; j++)
		{
			k = vector[i][j];

			sum[int(k)] = sum[int(k)] + pdf1[i] * pdf2[j];
		}
	}

	for(int i = 0; i < length_power_of_two; i++)
	{
		outVec4[i].rl = sum[i];
		outVec4[i].im = 0;
	}

	//quantize_number	

	return true;
}



double *new_convolution_LUT_initialization(double delta, double range)
{
	int length = range * 2.0 / delta;
	int LUT_length = length * length;
	double a[length];
	double b[length];
	double LUT_R[LUT_length];
	double R_before_quantize;
	double Q_value;
	int i = 0;

	for(int a = -range; a < range; a = a + delta)
	{
		for(int b = -range; b < range; b = b + delta)
		{
			R_before_quantize = 2 * atanh(tanh(0.5 * a) * tanh(0.5 * b));
			if(R_before_quantize >= delta / 2)                  //Is this delta equals to the delta above?
				Q_value = delta * floor(R_before_quantize/delta + 1/2);
			else if(R_before_quantize <= -delta / 2)
				Q_value = delta * ceil(R_before_quantize/delta - 1/2);
			else
				Q_value = 0;

			i = i + 1;
		}
	}

	return LUT_R;
}




/*
int main()
{
	double  result = 0;

	double  delta = QUANTIZE_DELTA;
	double  delta_awgn = DELTA_AWGN;
	double  delta_rayleigh = DELTA_RAYLEIGH;

	double  min = AWGN_CHANNEL? -DELTA_AWGN : -DELTA_RAYLEIGH;
	double  max = AWGN_CHANNEL?  DELTA_AWGN : DELTA_RAYLEIGH;

	int     length = floor(max/QUANTIZE_DELTA) - ceil(min/QUANTIZE_DELTA) + 1;
	double  *pdf1 = new double[length];
	double  *pdf2 = new double[length];
	int     k = 0;
	int     length_power_of_two;

	// convert the length into a number which is the power of two for fft
	length_power_of_two = convert_power_of_two(length, verbose);

	result = integral(pdf_awgn_func, -delta_awgn, delta_awgn, 1);
	printf("integral of pdf_awgn_func is %f\n", result);
	result = integral(pdf_Rayleigh_func, -delta_rayleigh, delta_rayleigh, 1);
	printf("integral of pdf_Rayleigh_func is %f\n", result);
	//
	Complex *inVec  = new Complex[length_power_of_two];  //Complex data for pdf1
	Complex *outVec = new Complex[length_power_of_two];
	Complex *invert = new Complex[length_power_of_two];

	Complex *inVec3  = new Complex[length_power_of_two]; //Complex data for pdf2
	Complex *outVec3 = new Complex[length_power_of_two];
	Complex *invert3 = new Complex[length_power_of_two];

	//multiply result
	Complex *inVec4  = new Complex[length_power_of_two];
	Complex *outVec4 = new Complex[length_power_of_two];
	Complex *invert4 = new Complex[length_power_of_two];		

	pdf1 = quantize_message(pdf_Rayleigh_func, 1, delta, min, max); //Complex data obey Rayleigh distribution
	pdf2 = quantize_message(pdf_Rayleigh_func, 1, delta, min, max);

	printf("length is %d\n", length_power_of_two);

	memset(inVec, 0, length_power_of_two * sizeof(Complex));
	memset(inVec3, 0, length_power_of_two * sizeof(Complex));

	//append data into input vector of fft
	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec[i].rl = pdf1[i];
		else
			inVec[i].rl = 0;
	}

	//append data into input vector of fft
	for (int i = 0; i < length_power_of_two; i++)
	{
		if(i < length)
			inVec3[i].rl = pdf2[i];
		else
			inVec3[i].rl = 0 ;
	}	

	//display input vector
	for(int l = 0; l < length_power_of_two; l++)
	{
		printf("pdf1 [%d] is %f \n",l,pdf1[l]);
	}	
	
	//double vec[] = {15, 32, 9, 222, 118, 151, 5, 7, 56, 233, 56, 121, 235, 89, 98, 111};
	double vec[] = {15, 32, 9, 222, 118, 151, 5, 7, 56, 233, 56, 121, 235, 89, 0, 0};
	int len = sizeof(vec) / sizeof(double);

	Complex *inVec2  = new Complex[len];
	Complex *outVec2 = new Complex[len];
	Complex *invert2 = new Complex[len];

 	// Append vector function
	memset(inVec2, 0, len * sizeof(Complex));
	for (int i = 0; i < len; i++)
		inVec2[i].rl = vec[i];

	// Fourier transformation
	Fft t,t2,t3;
	t.fft(inVec, length_power_of_two, outVec);
	t2.fft(inVec2, len, outVec2);
	t3.fft(inVec3, length_power_of_two, outVec3);
	
	// print result
	OutputDebugStringA("快速傅里叶变换结果为：\n");
	for (int i = 0; i < length_power_of_two; i++) {
		if (outVec[i].im < 0)
			printf("result[%d]: %lf - %lfi \n", i + 1, outVec[i].rl, -outVec[i].im);
		else
			printf("result[%d]: %lf + %lfi \n", i + 1, outVec[i].rl, outVec[i].im);
	}

	// multiply pmf and pmf2
	//for (int i = 0; i < length_power_of_two; i++) {
	//	outVec4[i].rl = outVec[i].rl * outVec3[i].rl - outVec[i].im * outVec3[i].im;
	//	outVec4[i].im = outVec[i].rl * outVec3[i].im + outVec[i].im * outVec3[i].rl;
	//}
	t3.multiply(outVec,outVec3,length_power_of_two, outVec4);

	// print result
	OutputDebugStringA("乘积为：\n");
	for (int i = 0; i < length_power_of_two; i++) {
		if (outVec4[i].im < 0)
			printf("result[%d]: %lf - %lfi \n", i + 1, outVec4[i].rl, -outVec4[i].im);
		else
			printf("result[%d]: %lf + %lfi \n", i + 1, outVec4[i].rl, outVec4[i].im);
	}
	
	OutputDebugStringA("逆变换结果为：");
	t.ifft(outVec, length_power_of_two, invert);

	for (int i = 0; i < length_power_of_two; i++) {
		printf("ifft[%d]: %lf \n", i + 1, invert[i].rl);
	}

	OutputDebugStringA("乘积的逆变换结果为：");
	t3.ifft(outVec4, length_power_of_two, invert4);

	for (int i = 0; i < length_power_of_two; i++) {
		//printf("ifft[%d]: %lf \n", i + 1, invert4[i].rl);
		if (invert4[i].im < 0)
			printf("ifft[%d]: %lf - %lfi \n", i + 1, invert4[i].rl, -invert4[i].im);
		else
			printf("ifft[%d]: %lf + %lfi \n", i + 1, invert4[i].rl, invert4[i].im);
	}
	
	//print result2
	OutputDebugStringA("快速傅里叶变换结果为：\n");
	for (int i = 0; i < len; i++) {
		if (outVec2[i].im < 0)
			printf("result[%d]: %lf - %lfi \n", i + 1, outVec2[i].rl, -outVec2[i].im);
		else
			printf("result[%d]: %lf + %lfi \n", i + 1, outVec2[i].rl, outVec2[i].im);
	}
	
	OutputDebugStringA("逆变换结果为：");
	t.ifft(outVec2, len, invert);
	for (int i = 0; i < len; i++) {
		printf("ifft[%d]: %lf \n", i + 1, invert[i].rl);
	}

	delete [] inVec;
	delete [] outVec;
	delete [] invert;
	delete [] inVec2;
	delete [] outVec2;
	delete [] invert2;	
	delete [] pdf1;
	delete [] pdf2;
	//system("pause");
	return 0;
}



//check integral of adjustment to be 1

result = integral(pdf_awgn_func, -delta_awgn, delta_awgn, 0, 1);
while(abs(result- 1) > VARIANCE_DELTA)
{
	result = integral(pdf_awgn_func, -delta_awgn, delta_awgn, 0, 1);
	cout << result << endl;
	if (result > 1)
		delta_awgn = delta_awgn - 0.01 / (k+1);
	else if (result == 1)
		break;
	else
		delta_awgn = delta_awgn + 0.01 / (k+1);		

	k = k + 1;
}

printf("delta_awgn adjustment end \n");

k = 0;
while(abs(result- 1) > VARIANCE_DELTA)
{
	result = integral(pdf_Rayleigh_func, -delta_rayleigh, delta_rayleigh, 0, 1);
	cout << result << endl;
	if (result > 1)
		delta_rayleigh = delta_rayleigh - 0.01 / (k+1);
	else if (result == 1)
		break;
	else
		delta_rayleigh = delta_rayleigh + 0.01 / (k+1);

	k = k + 1;
}

printf("delta_rayleigh adjustment end \n");

cout << delta_awgn << endl;
cout << delta_rayleigh << endl;


bool RCA_predict(HB,SNR,p,R,snr_R,infor)
bool RCA_predict(Matrix<double, Dynamic , Dynamic> HB, double pdf,  pdf, double SNR, double p, double R[5001], double snr_R[5001], double infor)
{
	bool predict_flag;

	double s;
	//123123;
	int iter_num = 50;

	int mvout_threshold = 100;

	int row = HB.rows();
	int col = HB.cols();

	s = 2 * pow(10,SNR / 10);

	for (int k = 0; k < iter_num; ++k)
	{

		for (int m = 0; i < m; ++m)
		{
			12313
		}


		for (int n = 0; j < n; ++n)
		{
			123123
			
		}
	}
	
	return predict_flag;
}

*/


