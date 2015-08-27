/*
 * NoisySeq2.h
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef NOISYSEQ2_H_
#define NOISYSEQ2_H_

#include "LDPC_Simulator_Data_Def.hpp"

class NoisySeq{
public:
	int Blocklength; /*!< blocklength of code */
	double *NoisySequence; /*!< noisy sequence vector, can be R^n for AWGN or {0,1}^n for BSC */
	string ChannelName; /*!< channel type: AWGN or BSC */
	double ChannelParameter;
	//! Print function. Intended for debugging.
	void print()
	{
		for(int i = 0; i < Blocklength; i++)
		{
			cout<<NoisySequence[i]<<" ";
		}
		cout<<endl;
	}

	//! Constructor.
    /*!
      \param blocklength Blocklength of code
	  \param channel Channel name, can only be "AWGN" or "BSC"
    */
	NoisySeq(int blocklength);
	void SetProperties(string channel, double channel_para){ChannelName = channel; ChannelParameter = channel_para;}
	//! destructor.
	~NoisySeq();
};


//Noisy sequence Class
NoisySeq::NoisySeq(int blocklength)
{
	Blocklength = blocklength;
	NoisySequence = new double[blocklength];
}

NoisySeq::~NoisySeq()
{
#ifdef DEBUG
cout<<"~NoisySeq()"<<endl;
#endif
	delete [] NoisySequence;
}

#endif /* NOISYSEQ2_H_ */
