/*
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 *
 * stable version: june 24th, 2015
 * still under modification
 */
#include <x86intrin.h>
#include <xmmintrin.h>

#ifndef DECODER_HPP_
#define DECODER_HPP_

//#define PROFILE_ON

#include "MatriceProjection.hpp"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


#define _type_ float

class Decoder{
public:
	//_type_* vector_before_proj;
	#define NEW_PROJECTION
	MatriceProjection<_type_, 32> mp;

	int get_deg_vn(int pos){
		return VariableDegree[pos];
	}

private:

	string DecoderName;
	string ChannelName;
	string mRealDecoderName;

	double* mProjectionPPd   (NODE v[], int length);
	double* mProjectionPPd_d6(NODE v[]);
	double* mProjectionPPd_d7(NODE v[]);

    long long int t_vn = 0;
    long long int t_cn = 0;
    long long int t_pj = 0;
    long long int t_ex = 0;

	// parameters for decoders
	string mPCMatrixFileName;  /*!< parity check matrix file */
	int mNChecks, mBlocklength, mPCheckMapSize;
	int *CheckDegree;
	int *VariableDegree;
	double RateOfCode; /*!< code rate, k/n */

	TRIPLE *mPCheckMap; /*!< parity check mapping */
	TRIPLE *mParityCheckMatrix; /*!< parity check matrix */

	unsigned int* t_col;
	unsigned int* t_row;
	unsigned int* t_col1;
	unsigned int* t_row1;

	float *Lambda;
	float *zReplica;
	float *latestProjVector;

	// algorithm parameter
	double para_end_feas;   /*!< ADMM alg end tolerance, typically 1e-5 */
	double para_mu;         /*!< ADMM alg \mu, typically 5.5 */
	double para_rho;        /*!< ADMM over relaxation para, typically 1.8, 1.9 */
	int maxIteration;       /*!< max iteration */

	//! Learn the degree of the parity check matrix
	/*!
	  The parity check matrix should be in correct format. This function is useful for irregular LDPC codes.
	*/
	void mLearnParityCheckMatrix();
	//! Read and store the parity check matrix
	void mReadParityCheckMatrix();
	void mSetDefaultParameters();


	//! Calculate log likelihood ratio using the output from the channel
	/*!
	  \param channeloutput Output from channel.
	  \sa NoisySeq
	*/
	void mGenerateLLR(NoisySeq &channeloutput);

	//Decoders:
	void ADMMDecoder();
	void ADMMDecoder_float();
	void ADMMDecoder_float_coeffs();

	//
	//
	//

	// data
	float *_LogLikelihoodRatio; /*!< log-likelihood ratio from received vector */
	float *OutputFromDecoder; /*!<soft information from decoder. could be pseudocodeword */

	double alpha;       /*!< the constant for penalty term */

	float  o_value[32]; /*!< the constant for penalty term */
	bool   coeff_mode;

	//decoder stats
	unsigned long int mExeTime; /*!< exevution time */
	bool mAlgorithmConverge; /*!< true if the decoder converges*/
	int mIteration; /*!< number of iterations used for decoding */
	bool mValidCodeword; /*!< true if the output is a valid codeword */

public:
	//! Set decoder parameters.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	*/
	void SetParameters(int mxIt, double feas, double p_mu, double p_rho);
	void SetParameters(int mxIt, double feas, float p_mu, float p_value[], int p_degs[], int n_degs);
	void SetDecoder(string decoder);
	void SetChannel(string channel){ChannelName = channel;}
	void SetPenaltyConstant(double alp){alpha = alp;}
	bool ValidateCodeword();
	//! See if the input sequence is valid codeword
	/*!
		\param input Input sequence array.
	*/
	bool ValidateCodeword(int input[]);
	//! Invoke decode function
	/*!
		\param channeloutput Noisy sequence from channel.
		\param decoderesult Decode results. Contains decoded sequence and other information.
		\sa DecodedSeq
		\sa NoisySeq
	*/
	double Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult);

	//! Constructor.
	/*!
		\param FileName File for parity check matrix.
		\param blocklength Blocklength of code (i.e. n)
		\param nChecks Number of checks (i.e. row number of parity check matrix or (n-k))
	*/
	Decoder(string FileName,int nChecks, int BlockLength);
	//! Destructor
	~Decoder();
};


// Decoder Class
Decoder::Decoder(string FileName,int nChecks, int BlockLength) : mp()
{
#ifdef PROFILE_ON
    t_vn = 0;
    t_cn = 0;
    t_pj = 0;
    t_ex = 0;
#endif

    coeff_mode = false; // BLG

	// Learn and read parity check matrix
	mPCMatrixFileName = FileName;
	mNChecks = nChecks;
	mBlocklength = BlockLength;
	mPCheckMapSize = 0;
	CheckDegree         = (int*)_mm_malloc(mNChecks * sizeof(int), 64);
	VariableDegree      = (int*)_mm_malloc(mBlocklength * sizeof(int), 64);
	OutputFromDecoder   = (float*)_mm_malloc(mBlocklength * sizeof(float), 64);
	_LogLikelihoodRatio = (float*)_mm_malloc(mBlocklength * sizeof(float), 64);

	mSetDefaultParameters();
	mLearnParityCheckMatrix();
	mPCheckMap = new TRIPLE[mPCheckMapSize];
	mParityCheckMatrix = new TRIPLE[mPCheckMapSize];
	mReadParityCheckMatrix();
	RateOfCode = (double)(mBlocklength -mNChecks)/mBlocklength;
	//vector_before_proj = (float*)_mm_malloc(64 * sizeof(float), 64);

	Lambda   = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);
	zReplica = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);
	latestProjVector = (float*)_mm_malloc(mPCheckMapSize * sizeof(float), 64);

	t_col              = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mParityCheckMatrix[i].row == k )
			{
				t_col[ptr] = mParityCheckMatrix[i].col;
				ptr += 1;
			}

		}
	}

	t_col1             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int k = 0; k < mPCheckMapSize; k++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].row == k )
			{
				t_col1[ptr] = mPCheckMap[i].col;
				ptr += 1;
                                break; // on en a fini avec le msg_k
			}

		}
	}

	t_row             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int j = 0; j < mBlocklength; j++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].col == j )
			{
				t_row[ptr] = mPCheckMap[i].row;
				ptr += 1;
			}
		}
	}
	t_row1             = (unsigned int*)_mm_malloc(mPCheckMapSize * sizeof(unsigned int), 64);
        ptr = 0;
	for(int j = 0; j < mBlocklength; j++)
	{
		for(int i = 0; i < mPCheckMapSize; i++)
		{
			if( mPCheckMap[i].col == j )
			{
				t_row1[ptr] = mPCheckMap[i].row;
 				ptr += 1;
			}
		}
	}
}

Decoder::~Decoder()
{
#ifdef PROFILE_ON
    if( t_ex != 0 ){
        t_cn /= t_ex;
        t_vn /= t_ex;
        t_pj /= t_ex;
        long long int sum = t_cn + t_vn;
        double CN = (100.0 * t_cn) / (double)sum;
        double VN = (100.0 * t_vn) / (double)sum;
        double PJ = (100.0 * t_pj) / (double)t_cn;
        printf("VN : %1.3f - CN = %1.3f (PJ = %1.3f)\n", VN, CN, PJ);
    }else{
    	printf("(WW) PROFILE WAS SWITWED ON AND NO DATA ARE AVAILABLE !\n");
    }
#endif

	_mm_free ( CheckDegree );
	_mm_free ( VariableDegree );
	delete [] mPCheckMap;
	delete [] mParityCheckMatrix;
	_mm_free ( OutputFromDecoder );
	_mm_free (_LogLikelihoodRatio );
	//_mm_free ( vector_before_proj );
	_mm_free ( t_col );
	_mm_free ( t_row );
	_mm_free ( t_col1 );
	_mm_free ( t_row1 );
	_mm_free ( Lambda );
	_mm_free ( zReplica );
	_mm_free( latestProjVector );
}

//Decodersettings
void Decoder::SetParameters(int mxIt, double feas,  double p_mu, double p_rho)
{
	maxIteration  = mxIt;
	para_end_feas = feas;
	para_mu       = p_mu;
	para_rho      = p_rho;
}

void Decoder::SetParameters(int mxIt, double feas, float p_mu, float p_value[], int p_degs[], int n_degs)
{
	maxIteration  = mxIt;
	para_end_feas = feas;
	para_mu       = p_mu;
	para_rho      = 1.90f; // BLG
	coeff_mode    = true;

	for(int i=0; i<32; i++)
	{
		o_value[ i ] = 0.0f;
	}
	for(int i=0; i<n_degs; i++)
	{
		o_value[ p_degs[i] ] = p_value[i] / p_mu;
	}
}



void Decoder::SetDecoder(string decoder)
{
	mRealDecoderName = decoder;
        DecoderName  = "ADMM";
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mSetDefaultParameters()
{
	// decoder settings
	maxIteration = 1000;
	para_end_feas = 1e-6;
	para_mu = 5;
	para_rho= 1.8;

	alpha = 1;
	for(int i = 0; i < mBlocklength; i++)
	{
		VariableDegree[i] = 0;
	}
	for(int i = 0; i < mNChecks; i++)
	{
		CheckDegree[i] = 0;
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mLearnParityCheckMatrix()// learn the degree of the parity check matrix
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	string line;
	if (myfile.is_open())
	{
		while(getline(myfile, line))
		{
		   istringstream is(line);
		   int curr_degree = 0;
		   while( is >> tempvalue )
		   {
			   VariableDegree[tempvalue]++;
			   curr_degree++;
			   mPCheckMapSize++;
		   }
		   CheckDegree[i] = curr_degree;
		   i++;
		}
		if (myfile.is_open())
			myfile.close();
	}
	else
	{
		cout << "Unable to open file";
		exit(0);
	}
//	cout <<"mPCheckMapSize= "<< mPCheckMapSize<<endl;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mReadParityCheckMatrix()// read the parity check matrix, initialize message passing
{
	ifstream myfile (mPCMatrixFileName.c_str());
	int tempvalue = 0;
	int i = 0;
	int count = 0;
	string line;
	if (myfile.is_open())
	{
		while(getline(myfile, line))
		{
			istringstream is(line);
			while( is >> tempvalue )
			{
				mParityCheckMatrix[count].col = tempvalue;
				mParityCheckMatrix[count].row = i;

				mPCheckMap[count].col = tempvalue;
				mPCheckMap[count].row = count;
				count++;
			}
			i++;
		}
	}
	if (myfile.is_open())
		myfile.close();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Decoder::mGenerateLLR(NoisySeq &channeloutput)
{
	if(ChannelName == "AWGN")
	{
		double std = channeloutput.ChannelParameter;
		for(int i = 0; i < mBlocklength; i++)
		{
			double receivedbit = channeloutput.NoisySequence[i];
			_LogLikelihoodRatio[i] = ((receivedbit - 1)*(receivedbit - 1) - (receivedbit + 1)*(receivedbit + 1))/2.0/std/std;
		}
	}
	else if(ChannelName == "BSC")
	{
		double CrossOverProb = channeloutput.ChannelParameter;
		for(int i = 0; i < mBlocklength; i++)
		{
			int receivedbit = floor(channeloutput.NoisySequence[i] + 0.5);
			if(receivedbit == 0)
				_LogLikelihoodRatio[i] = - log(CrossOverProb/(1 - CrossOverProb));
			else
				_LogLikelihoodRatio[i] = log(CrossOverProb/(1 - CrossOverProb));
		}
	}
	else
	{
		cout<<"channel not supported"<<endl;
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double Decoder::Decode(NoisySeq &channeloutput, DecodedSeq &decoderesult)
{
	mExeTime           = 0;
	mAlgorithmConverge = false;
	mValidCodeword     = false;

	mGenerateLLR(channeloutput);

	clock_t tStart = clock();
    auto start     = chrono::steady_clock::now();

    //
    //
    //
    if( coeff_mode == true ){
        ADMMDecoder_float_coeffs( );

    }else{
            ADMMDecoder_float();
    }

    auto end                       = chrono::steady_clock::now();
	mExeTime                       = clock() - tStart;
	decoderesult.ExeTime           = mExeTime;
	decoderesult.Iteration         = mIteration;
	decoderesult.ValidCodeword     = mValidCodeword;
	decoderesult.AlgorithmConverge = mAlgorithmConverge;
	decoderesult.Iteration         = mIteration;

	for(int i = 0; i < mBlocklength; i++)
	{
		if(my_isnan(OutputFromDecoder[i]))
			decoderesult.IsNan = true;
		if(my_isinf(OutputFromDecoder[i]))
			decoderesult.IsInf = true;
		decoderesult.SoftInfo[i] = OutputFromDecoder[i];
		decoderesult.HardDecision[i] = floor(OutputFromDecoder[i] + 0.5);
	}
    auto diff      = end - start;
    return chrono::duration <double, milli> (diff).count();
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void show6(float* a)
{
	printf("__m256 : %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n", a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show8(float* a)
{
	printf("__m256 : %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f\n", a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show6(int* a)
{
	printf("__m256 : %d %d %d %d %d %d\n", a[5], a[4], a[3], a[2], a[1], a[0]);
}

void show8(int* a)
{
	printf("__m256 : %d %d %d %d %d %d %d %d\n", a[7], a[6], a[5], a[4], a[3], a[2], a[1], a[0]);
}

#include "decoders/decoder_float.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool Decoder::ValidateCodeword()
{
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		int syndrom = 0;
		int deg     = CheckDegree[k];
		for(int i = 0; i < deg; i++)
		{
			syndrom += (int) floor(OutputFromDecoder[ t_col[ptr] ] + 0.5);
			ptr += 1;
		}
		if( syndrom & 0x01 ) return false;
	}
	return true;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool Decoder::ValidateCodeword(int input[])
{
	int ptr = 0;
	for(int k = 0; k < mNChecks; k++)
	{
		int syndrom = 0;
		int deg     = CheckDegree[k];
		for(int i = 0; i < deg; i++)
		{
			syndrom += (int) floor(input[ t_col[ptr] ]);
			ptr += 1;
		}
		if( syndrom & 0x01 ) return false;
	}
	return true;
}

#endif /* DECODER_HPP_ */
