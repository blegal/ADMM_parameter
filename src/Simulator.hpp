/*
 * Simulator.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef SIMULATOR_HPP_
#define SIMULATOR_HPP_

#include "Decoder.hpp"
#include <chrono>

using namespace std;

//!  Class for simulation
/*!
Run simulations.
Take input sequence, pass it through channel to get noisy sequence. Pass noisy sequence to decoder and then compare decoded sequence with original input.
*/
class Simulator{
private:
	typedef unsigned long uint32;
	int *mInput;

    const int max_threads = 8;
    int num_threads = 1;
	NoisySeq*   mNoisySequence  [8];
	Channel*    mChannel        [8];
	DecodedSeq* mDecodedSequence[8];

public:
	Decoder*    mDecoder        [8];

	string mChannelName;
	string mDecoderName;

	string mParityCheckFile;
	string mCodewordFile;

	int    mBlocklength;
	int    mNChecks;
	double mRateOfCode;
	double EbN0;
	double mSTD;
	double mCrossOverProb;
	double mChannelParameter;
	bool   mErrorFlag;

    long long tick = 0;

    //
	// Simulation control
    //
	uint32   mOutputFileInterval;
	string   mOutputFileName;
	ofstream mOutputFileStream;

    //
    //
    //
	uint32 mOutputCommandInterval;
	uint32 mOutputCommandDetailLevel;

    //
	// simulation stats
    //
	uint32 mTargetErrors;
	uint32 mTargetSims;
	uint32 mTargetExeTime;

	double mTotalExeTime;
	double mTotalSimTime;

	uint32 nbBitErrors;//imen BER calc

	uint32 mTotalErrors;
	uint32 mTotalSims;
	uint32 mTotalDecodingTime;
	uint32 mTotalIterations;

	uint32 mTotalCorrectIterations;
	uint32 mTotalWrongIterations;

	uint32 mTotalCorrectDecodingTime;
	uint32 mTotalWrongDecodingTime;

	uint32 mTotalMLErrors;
	uint32 mTotalUndetectedErrors;



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	//! The input is Eb/N0 (dB). This should be converted to std in AWGN or crossover probability in BSC
	void mTranslateEbN0()
	{
        mSTD           = (double)sqrt(1.0 / pow(10.0, EbN0/10) / 2.0 / mRateOfCode);
        double temp    = sqrt(2 * mRateOfCode * pow(10.0 , EbN0/10));
        mCrossOverProb = 0.5 - 0.5 * myerf(temp/sqrt(2.0));
        if(mChannelName == "AWGN")
            mChannelParameter = mSTD;
        else
            mChannelParameter = mCrossOverProb;
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	void mUpdateStats(int qq);



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	void mClear()
	{
        nbBitErrors=0;//imen BER calc
		mTotalErrors              = 0;
		mTotalSims                = 0;
		mTotalExeTime             = 0;
		mTotalSimTime             = 0;
		mTotalDecodingTime        = 0;
		mTotalIterations          = 0;
		mTotalCorrectIterations   = 0;
		mTotalWrongIterations     = 0;
		mTotalCorrectDecodingTime = 0;
		mTotalWrongDecodingTime   = 0;
		mTotalMLErrors            = 0;
		mTotalUndetectedErrors    = 0;
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	void ShowTime(unsigned long secondes)
	{
	    int ss = secondes % 60;
	    int mn = (int)round(secondes / 60.0f) % 60;
	    int hh = (int)round(secondes / 3600.0f);
	    printf("%2.2dh%2.2d'%2.2d", hh, mn, ss);
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	void mOutputCommandFlush( )
	{
//		if (mTotalSims % mOutputCommandInterval == 0)
//		{
			if(mOutputCommandDetailLevel == 3)
			{
			    cout<<"  TotalSims = "       <<mTotalSims
				<<"  Time = "            << (int)mTotalExeTime << "s  Total error = "<<mTotalErrors
				<<"  Avg # iteration = " << mTotalIterations/(double)(mTotalSims)
				<<"  Undetected error = "<<mTotalUndetectedErrors
				<<"  ML error = "        <<mTotalMLErrors
				<<"  Avg decode time = " <<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
				<<endl;
			}
			else if(mOutputCommandDetailLevel == 2)
			{
			    cout<<"  TotalSims = "               << mTotalSims
				<<"  Time = "                    <<(int)mTotalExeTime
				<<"  Total error = "             <<mTotalErrors
				<<"  Undetected error = "        << mTotalUndetectedErrors
				<<"  ML error = "                <<mTotalMLErrors
				<<"  Avg # iteration = "         <<mTotalIterations/(double)(mTotalSims)
				<<"  Avg # iteration correct = " <<mTotalCorrectIterations/(double)(mTotalSims - mTotalErrors)
				<<"  Avg # iteration wrong = "   <<mTotalWrongIterations/(double)(mTotalErrors)
				<<"  Avg decode time = "         <<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
				<<"  Avg decode time correct = " <<(double)mTotalCorrectDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims - mTotalErrors) <<"(s)"
				<<"  Avg decode time wrong = "   <<(double)mTotalWrongDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalErrors) <<"(s)"
				<<endl;
			}
			else if(mOutputCommandDetailLevel != -1)
			{
				//		}
				//			printf("mTotalExeTime = %f milli-sec.\n", mTotalExeTime);
				// imen : consider all bits for BER calc
				double tBER = ((double)nbBitErrors) / (mBlocklength * (double)mTotalSims);
				// imen
				//double tBER = (mBlocklength * (double)mTotalErrors) / (mBlocklength * (double)mTotalSims);
				double tFER          = (double)mTotalErrors / (double)mTotalSims;
				double temps         = (double)mTotalExeTime;
				double sec_par_trame = (temps/1000.f) / mTotalSims;
				double trame_par_sec = 1.0f / sec_par_trame;
				double        bps    = (trame_par_sec * mBlocklength) / 1000.0 / 1000.0;
				//		    double be_par_fe     = mTotalUndetectedErrors;
				double avgIters      = (double)mTotalIterations/(double)(mTotalSims);
				printf("SNR = %.2f | BER = %1.3e | FER = %1.3e | BPS = %2.2f | MATRICES = %8ld | FE = %d | ERR = %d | MlErr = %d | Avg#Iter = %1.2f | ", EbN0, tBER, tFER, bps, mTotalSims, (int) mTotalErrors, (int) mTotalUndetectedErrors, (int)mTotalUndetectedErrors, avgIters);
				ShowTime( temps/1000.f );
				printf(" | ");
				ShowTime( (double)mTotalSimTime/1000.f );
				printf("\n");
				fflush(stdout);
			}
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	void mOutputCommand()
	{
		if (mTotalSims % mOutputCommandInterval == 0)
		{
			if(mOutputCommandDetailLevel == /*1*/ 3)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"*"<<(mTotalDecodingTime*100)/mTotalExeTime<<"%(s)  Total error = "<<mTotalErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<endl;
			}
			else if(mOutputCommandDetailLevel == 2)
			{
				cout<<"TotalSims = "<<mTotalSims
					<<"  Time = "<<(double)mTotalExeTime/CLOCKS_PER_SEC<<"*"<<(mTotalDecodingTime*100)/mTotalExeTime<<"%(s)"
					<<"  Total error = "<<mTotalErrors
					<<"  Undetected error = "<<mTotalUndetectedErrors
					<<"  ML error = "<<mTotalMLErrors
					<<"  Avg # iteration = "<<mTotalIterations/(double)(mTotalSims)
					<<"  Avg # iteration correct = "<<mTotalCorrectIterations/(double)(mTotalSims - mTotalErrors)
					<<"  Avg # iteration wrong = "<<mTotalWrongIterations/(double)(mTotalErrors)
					<<"  Avg decode time = "<<(double)mTotalDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims) <<"(s)"
					<<"  Avg decode time correct = "<<(double)mTotalCorrectDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalSims - mTotalErrors) <<"(s)"
					<<"  Avg decode time wrong = "<<(double)mTotalWrongDecodingTime/CLOCKS_PER_SEC/ (double)(mTotalErrors) <<"(s)"
					<<endl;
			}
			else if(mOutputCommandDetailLevel == -1)
			{
				/* */
			}
            else
            {
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                double temps = (double)mTotalExeTime;
                if( mTotalErrors == 0)
                {
                    double tFER          = (double)1.0f / (double)mTotalSims;
                    double tBER          = ((double)1.0f) / (mBlocklength * (double)mTotalSims);
                    double sec_par_trame = (temps/1000.f) / mTotalSims;
                    double trame_par_sec = 1.0f / sec_par_trame;
                    double        bps    = (trame_par_sec * mBlocklength) / 1000.0 / 1000.0;
                    double avgIters      = (double)mTotalIterations/(double)(mTotalSims);
                    printf("?NR = %.2f | BER < %1.3e | FER < %2.2e | BPS = %2.2f | MATRICES = %8ld | FE = %3d | #Iter = %1.2f | fERR = %d | MlErr = %d | ", EbN0, tBER, tFER, bps, mTotalSims, (int) mTotalErrors, avgIters, (int) mTotalUndetectedErrors, (int)mTotalUndetectedErrors);

                }else{
                    double tFER          = (double)mTotalErrors / (double)mTotalSims;
                    double tBER          = ((double)nbBitErrors) / (mBlocklength * (double)mTotalSims);
                    double sec_par_trame = (temps/1000.f) / mTotalSims;
                    double trame_par_sec = 1.0f / sec_par_trame;
                    double        bps    = (trame_par_sec * mBlocklength) / 1000.0 / 1000.0;
                    double avgIters      = (double)mTotalIterations/(double)(mTotalSims);
                    printf("?NR = %.2f | BER = %1.3e | FER = %2.2e | BPS = %2.2f | MATRICES = %8ld | FE = %3d | #Iter = %1.2f | fERR = %d | MlErr = %d | ", EbN0, tBER, tFER, bps, mTotalSims, (int) mTotalErrors, avgIters, (int) mTotalUndetectedErrors, (int)mTotalUndetectedErrors);

                }
                ShowTime( temps/1000.f );
    		    printf(" | ");
    		    ShowTime( (double)mTotalSimTime/1000.f );
                if( tick%4 == 0 ) printf("   ");
                if( tick%4 == 1 ) printf(".  ");
                if( tick%4 == 2 ) printf(".. ");
                if( tick%4 == 3 ) printf("...");
                tick += 1;
                fflush( stdout );
                printf("\r");
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
		}
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



public:
    double GetLatestFerValue( ){
    	return (double)GetTotalErrors() / (double)GetTotalSims();
    }

	//! Return total number of errors
	uint32 GetTotalErrors(){return mTotalErrors;}
	//! Return total number of simulations
	uint32 GetTotalSims(){return mTotalSims;}
	//! Return time used for decoding
	uint32 GetTotalDecodingTime(){return mTotalDecodingTime;}
	//! Return total number of iterations
	uint32 GetTotalIterations(){return mTotalIterations;}
	//! Return total number of iterations for correct decoding events
	uint32 GetTotalCorrectIterations(){return mTotalCorrectIterations;}
	//! Return total number of iterations for erroneous decoding events
	uint32 GetTotalWrongIterations(){return mTotalWrongIterations;}
	//! Return time used for correct decoding events
	uint32 GetTotalCorrectDecodingTime(){return mTotalCorrectDecodingTime;}
	//! Return time used for erroneous decoding events
	uint32 GetTotalWrongDecodingTime(){return mTotalWrongDecodingTime;}
	//! Return total number of errors that are sure to be ML errors.
	/*!
	The ML error is calculated whenever there is a decoded codeword but is not the ML error. This is a lower bound for number of ML errors.
	*/
	uint32 GetTotalMLErrors(){return mTotalMLErrors;}
	//! Return total number of undetected errors. These are valid codeword but not the transmitted codeword.
	uint32 GetTotalUndetectedErrors(){return mTotalUndetectedErrors;}
	//! Return time for simulation. This include decoding time and time for channels, etc.
	uint32 GetTotalExeTime(){return mTotalExeTime;}
	//! Set channel's random number seed
	/*!
	\param seed Seed for channel.
	*/
	void SetChannelSeed(int seed){
        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mChannel[qq]->SetSeed(seed + qq);
    }
	//! Set simulation targets
	/*!
		\param tError Target errors. Stop simulation if more than this number of errors are collected.
		\param tSim Target simulations. Stop simulation if more than this number of simulations are performed.
		\param tExeTime Target time. Stop simulation if more than this amount of time is consumed.
	*/
	void SetTargets(uint32 tError, uint32 tSim, uint32 tExeTime){mTargetErrors = tError; mTargetSims = tSim; mTargetExeTime = tExeTime;}
	//! Set codeword.
	/*!
		Currently there is no encoding function. But you can choose a valid codeword to simulation if you have one other than the all zero codeword.
		\param filename Codeword file. The file should contain codeword. It should be 0,1 sequence, broken by space. e.g. 0 1 0 0 1 0
	*/
	void SetCodeword(string filename);
	//! Set command line output appearance.
	/*!
		\param simInterval Interval for each command line output.
		\param detaillevel Detail level for output. Only takes values 1 and 2.
	*/
	void SetCommandLineOutput(int simInterval, int detaillevel){mOutputCommandInterval = simInterval; mOutputCommandDetailLevel = detaillevel;}

	//! Set decoder name.
	/*!
	\param decoder Decoder name.
	 \sa Decoder
	*/
	void SetDecoder(string decoder){
        mDecoderName = decoder;
        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]->SetDecoder(mDecoderName);
        }
    }

	//! Set decoder parameters. Sufficient for ADMM LP.
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho)
	{
        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]->SetParameters(mxIt, feas, p_mu, p_rho);
        }
	}
	//! Set decoder parameters. Sufficient for ADMM L1 and ADMM L2
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp)
	{
        for (int qq = 0; qq<num_threads; qq+=1){
            mDecoder[qq]->SetParameters(mxIt,feas,p_mu,p_rho);
        }

        for (int qq = 0; qq<num_threads; qq+=1){
        	// BLG
            mDecoder[qq]->SetPenaltyConstant(alp);
        }
	}

	void SetDecoderParameters(int mxIt, double feas, float p_mu, float p_value[], int p_degs[], int n_degs)
	{
        for (int qq = 0; qq<num_threads; qq+=1){
            mDecoder[qq]->SetParameters(mxIt, feas, p_mu, p_value, p_degs, n_degs);
        }

        for (int qq = 0; qq<num_threads; qq+=1){
        	// BLG
            mDecoder[qq]->SetPenaltyConstant(0.00f);
        }
	}

	//! Set decoder parameters. Sufficient for all ADMM algorithms
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param alp Penalty constant.
	 \param th Penalty thershold for entropy penalty
	 \sa Decoder
	*/
/*	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, double alp, double th)
	{
		mDecoder.SetParameters(mxIt,feas,p_mu,p_rho); mDecoder.SetPenaltyConstant(alp); mDecoder.SetEntropyParameter(th);
	}
	*/
	//! Set decoder parameters. Sufficient for BP decoding
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param nansat True if non-saturating BP is used
	 \param value Saturation value is original BP is used
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool nonsat, double value)
	{
        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]->SetParameters(mxIt, feas, p_mu, p_rho); //mDecoder.SetBPSaturation(nonsat, value);
        }
	}

    //! Set decoder parameters. Setting for normalized LLR for BSC
	/*!
	 \param mxIt Maximum number of iterations. Default: 1000
	 \param feas ADMM stopping parameter. Default: 1e-6
	 \param p_mu ADMM step parameter mu. Default: 5
	 \param p_rho ADMM over relaxation parameter. Default: 1.8
	 \param normal True if the LLRs from BSC are normalized to +1 or -1.
	 \sa Decoder
	*/
	void SetDecoderParameters(int mxIt, double feas,  double p_mu, double p_rho, bool normal)
	{
        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]->SetParameters(mxIt, feas, p_mu, p_rho); //mDecoder.SetBSCnormalized(normal);
        }
	}

    //! Run simulations
	void RunSim();

	//! Constructor
	/*!
		\param blocklength Blocklength of the code
		\param nchecks Number of checks
		\param pcfilename Parity check matrix file name.
		\param channelname Channel name.
		\param eb_over_nzero Eb/N0 for the channel.
	*/
	Simulator(int blocklength, int nchecks, string pcfilename, string channelname, double eb_over_nzero, int _num_threads)
//		: mNoisySequence(blocklength)
//		, mDecodedSequence(blocklength)
//		, mDecoder(pcfilename, nchecks, blocklength)
//		, mChannel(blocklength)
	{
        num_threads = _num_threads;

        /* ON INSTANCIE TOUS LES DECODEURS */
        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mNoisySequence[qq]   = new NoisySeq(blocklength);

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mDecodedSequence[qq] = new DecodedSeq(blocklength);

        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]         = new  Decoder(pcfilename, nchecks, blocklength);
        }

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mChannel[qq]         = new Channel(blocklength);


		mChannelName     = channelname;
		mBlocklength     = blocklength;
		mNChecks         = nchecks;
		mParityCheckFile = pcfilename;
		mRateOfCode      = (double)(mBlocklength - mNChecks)/mBlocklength;
		EbN0             = eb_over_nzero;
		mTranslateEbN0();

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mNoisySequence[qq]->SetProperties(mChannelName, mChannelParameter);

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            mChannel[qq]->SetChannel(mBlocklength, mChannelName, mChannelParameter);

        for (int qq = 0; qq<num_threads; qq+=1)
        {
            mDecoder[qq]->SetChannel(mChannelName);
        }

        mInput = new int[mBlocklength];
		for(int i = 0; i < mBlocklength; i++)
			mInput[i] = 0;
	}



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	~Simulator()
	{
        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            delete mNoisySequence[qq];

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            delete mDecodedSequence[qq];

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
        {
            delete mDecoder[qq];
        }

        for (int qq = 0; qq<num_threads; qq+=1) // BLG
            delete mChannel[qq];

        delete [] mInput;
	}

};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Simulator::SetCodeword(string codewordfile)// use user defined codeword
{
	ifstream cwfile(codewordfile.c_str());
	int temp_input = 0;
	int i = 0;
	if (cwfile.is_open())
	{
		while (cwfile >> temp_input)
		{
			mInput[i] = temp_input;
			i++;
			if(i >= mBlocklength)
				break;
		}
	}

    // BLG : on utilise le decoder 0
    if(!mDecoder[0]->ValidateCodeword(mInput))
	{
		cout<<"invalid c.w., use all zero instead"<<endl;
		for(int i = 0; i < mBlocklength; i++)
			mInput[i] = 0;
	}
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Simulator::mUpdateStats(int qq)
{

	mErrorFlag = false;
    // imen
    // if the code is systematic then consider only the information bits
    //for(int i = 0; i < (mBlocklength/2)+1; i++)
    for(int i = 0; i < mBlocklength; i++)
    {
        if(mInput[i] != mDecodedSequence[qq]->HardDecision[i])
        {
            mErrorFlag = true;
            break;
        } // decoded seq does not agree with the input
        if(mDecodedSequence[qq]->IsNan)
        {
            mErrorFlag = true;
            break;
        }// decoded seq is nan
        if(mDecodedSequence[qq]->IsInf){
            mErrorFlag = true;
            break;
        }// decoded seq is inf
    }
    // imen : calcul de BER
	for(int i = 0; i < mBlocklength; i++)
	{
		if(mInput[i] != mDecodedSequence[qq]->HardDecision[i])
        {
            nbBitErrors += 1;
        }
	}
    // imen end
	if(mErrorFlag)
	{
		mTotalErrors++;
		mTotalWrongIterations   += mDecodedSequence[qq]->Iteration;
		mTotalWrongDecodingTime += mDecodedSequence[qq]->ExeTime;
		if( mDecodedSequence[qq]->ValidCodeword )
		{
			mTotalUndetectedErrors++;
			if(mChannelName == "AWGN")
			{
				double dis_orig = 0;
				double dis_decode = 0;

                // imen
                // if the code is systematic then consider only the information bits

				//for(int i = 0; i < (mBlocklength/2)+1; i++)
				for(int i = 0; i < mBlocklength; i++)
				{
					dis_orig   += ((double)mInput[i] * 2.0 - 1 - mNoisySequence[qq]->NoisySequence[i])*((double)mInput[i]*2.0 - 1 - mNoisySequence[qq]->NoisySequence[i]);
					dis_decode += ((double)floor(mDecodedSequence[qq]->HardDecision[i] + 0.5)*2.0 - 1 - mNoisySequence[qq]->NoisySequence[i]) * ((double)floor(mDecodedSequence[qq]->HardDecision[i] + 0.5)*2.0 - 1 - mNoisySequence[qq]->NoisySequence[i]);
				}
				if (dis_orig > dis_decode)
				{
					mTotalMLErrors++;
				}
			}
			else if(mChannelName == "BSC")
			{
				double hmdis_orig = 0;
				double hmdis_decoder = 0;
				for(int i = 0; i < mBlocklength; i++)
				{
					hmdis_orig += abs(mInput[i] -  mNoisySequence[qq]->NoisySequence[i]);
					hmdis_decoder += abs(mDecodedSequence[qq]->HardDecision[i] -  mNoisySequence[qq]->NoisySequence[i]);
				}
				if (hmdis_orig > hmdis_decoder)
				{
					mTotalMLErrors++;
				}
			}
		}
	}
	else
	{
		mTotalCorrectIterations   += mDecodedSequence[qq]->Iteration;
		mTotalCorrectDecodingTime += mDecodedSequence[qq]->ExeTime;
	}
	mTotalSims         += 1;
	mTotalDecodingTime += mDecodedSequence[qq]->ExeTime;
	mTotalIterations   += mDecodedSequence[qq]->Iteration;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Simulator::RunSim( )
{
	mTotalExeTime      = 0;
	mTotalSimTime      = 0;
    mTotalSims         = 0;
    mTotalDecodingTime = 0;
    mTotalIterations   = 0;

    auto start     = chrono::steady_clock::now();
    mClear(); // clear stats
	while(true)
	{

        //for (int qq = 0; qq<num_threads; qq+=1) // BLG

        double t_exec[max_threads];
        for (int qq = 0; qq<num_threads; qq+=1) t_exec[qq] = 0.0f;

//        auto mesure_debut     = chrono::steady_clock::now();

        #pragma omp parallel for num_threads(num_threads)
        for (int qq = 0; qq<num_threads; qq+=1) // BLG
        {
            for(int p=0; p<100; p++){
                mChannel[qq]->GenerateOutput(mInput, *mNoisySequence[qq]);
                double eTime = mDecoder[qq]->Decode(*mNoisySequence[qq], *mDecodedSequence[qq]);
                t_exec[qq] += eTime;
                #pragma omp critical // section critique: ne peut être executee que par un thread a la fois
                {
                    mUpdateStats(qq);
                }
            }
        }

//        auto mesure_fin = chrono::steady_clock::now();
//        auto diff       = mesure_fin - mesure_debut;
//        float measure  = chrono::duration <double, milli> (diff).count();

        int max_value = t_exec[0];
        for (int qq = 1; qq<num_threads; qq+=1)
            max_value = std::fmax(max_value, t_exec[qq]);

        mTotalExeTime += max_value;

//        printf("measure 1 = %d and Sum = %f\n", measure, max_value);

//        mTotalExeTime += _time_;//chrono::duration <double, milli> (diff).count();

        for (int qq = 0; qq<num_threads; qq+=1) // BLG

        mOutputCommand();

        mTotalSimTime = chrono::duration <double, milli> (chrono::steady_clock::now() - start).count();

        if(mTotalErrors >= mTargetErrors && mTargetErrors > 0)
        {
            break;
        }
		if(mTotalSims >= mTargetSims && mTargetSims > 0)
        {
            break;
        }
		if(mTotalExeTime >= mTargetExeTime && mTargetExeTime > 0)
        {
            break;
        }
	}

    if(mOutputCommandDetailLevel != -1){
    	mOutputCommandFlush( );
    }
}



#endif /* SIMULATOR_HPP_ */
