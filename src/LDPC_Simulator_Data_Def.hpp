// Copyright 2012 Xishuo Liu, Stark C. Draper, Benjamin Recht
//
// This program is distributed under the terms of the GNU General Public License.
//
// This file is part of ADMM Decoder.
//
//    ADMM Decoder is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ADMM Decoder is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ADMM Decoder.  If not, see <http://www.gnu.org/licenses/>.
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Project:				ADMM Decoder
// Files:				LDPC_Simulator_Example.cpp, LDPC_Class.cpp,
//						LDPC_Simulator_Data_Def.h, LDPC_Class.h,
//						MersenneTwister.h
// Date:				8.30.2012
//
// Author:				Xishuo Liu, xliu94@wisc.edu
// Thanks to:			S. Barman, S. Draper and B. Recht.
//
// Papers:				1. S. Barman, X. Liu, S. Draper and B. Recht,  
//						"Decomposition Methods for Large Scale LP Decoding"
//						http://arxiv.org/abs/1204.0556
//						2. X. Liu, S. Draper and B. Recht,
//						"Suppressing Pseudocodewords by Penalizing the Objective of LP Decoding"
//						IEEE Information Theory Workshop (ITW), 2012. Lausanne: Switzerland, 2012
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
This file contains basic data types and functions for the ADMM decoder classes. 
*/
#ifndef BIG_HEADER_HPP_
#define BIG_HEADER_HPP_

#include<iostream>
#include<cmath>
#include<ctime>
#include<fstream>
#include<string>
#include<iomanip>
#include<cstdlib>
#include<algorithm>
#include<sstream>

#include "MersenneTwister.hpp"
#include "MESSAGE.hpp"

using namespace std;


//!  Calculate erf function
/*!
  erf(z)     = 2/sqrt(pi) * Integral(0..x) exp( -t^2) dt
  erf(0.01)  = 0.0112834772 erf(3.7) = 0.9999998325
  Abramowitz/Stegun: p299, |erf(z)-erf| <= 1.5*10^(-7) 
*/
inline double myerf(double x)
{  
 double y = 1.0 / ( 1.0 + 0.3275911 * x);   
 return 1 - (((((
        + 1.061405429  * y
        - 1.453152027) * y
        + 1.421413741) * y
        - 0.284496736) * y 
        + 0.254829592) * y) 
        * exp (-x * x);      
}

//! True if x = Nan
inline int my_isnan(double x)
{
return x != x;
}

//! True if x = Inf
inline int my_isinf(double x)
{
if ((x == x) && ((x - x) != 0.0)) return (x < 0.0 ? -1 : 1);
else return 0;
}

//! generate random string
/*!
	\param length Length of string.
*/
inline string randomStrGen(int length) {
    static string charset = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
    string result;
    result.resize(length);

    for (int i = 0; i < length; i++)
        result[i] = charset[rand() % charset.length()];

    return result;
}

//!  A struct for sparse matrix structure
typedef struct Triple{
	int col;    /*!< column index. */
	int row;    /*!< row index.  */
	int value;  /*!< value index, reserved for non-binary codes*/
}TRIPLE;

//!  A struct for keeping sorting indicies 
/*!
	The idea is to obtain the indicies for sorting not just the final results for sorting.
*/
typedef struct  
{
	int    index;
	double value;
}NODE;

//! Sorting definition for struct NODE
bool sort_compare_de(const NODE & a, const NODE & b);

#endif
