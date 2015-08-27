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
This file contains declarations for the ADMM decoder classes. 
*/
#include "Channel.hpp"
#include "DecodedSeq.hpp"
#include "LDPC_Simulator_Data_Def.hpp"
#include "NoisySeq.hpp"
#include "Simulator.hpp"
