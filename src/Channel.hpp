/*
 * Channel.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef CHANNEL_HPP_
#define CHANNEL_HPP_

#include "MersenneTwister.hpp"
#include "NoisySeq.hpp"

#include<iostream>
#include<cmath>
#include<ctime>
#include<fstream>
#include<string>
#include<iomanip>
#include<cstdlib>
#include<algorithm>
#include<sstream>

using namespace std;

//!  Class for noisy sequence.
/*!
  Could be replaced by just an array ... anyway.
*/

//!  Class for channel.
/*!
  Takes input and generate output. There are two types of channel supported: AWGN and BSC.
  The random number generator uses the Mersenne Twister random number generator.
*/
class Channel{
private:

	int mBlocklength;  /*!< blocklength of code */
	string mChannelName; /*!< channel type: AWGN or BSC */

	MTRand channelrand; /*!< random number generator*/
	double mCrossoverProb;/*!< cross over probability, is used only when channel is BSC*/
	double mSTD;/*!< std, is used only when channel is AWGN */

	void mGenerateAWGN(int input[], NoisySeq &ChannelOutput); // AWGN by std
	void mGenerateBitFlip(int input[], NoisySeq &ChannelOutput); // bit flip channel by num of bit flips
public:
	//! Set Channel
	/*!
      \param blocklength Blocklength of code
	  \param channel channel name, can only be "AWGN" or "BSC"
	  \param parameter Parameter for the channel. std for AWGN or crossover prob for BSC
    */
	void SetChannel(int blocklength, string channel, double parameter);

	//! Generate output noisy sequence.
    /*!
      \param input Vector for input. Should contains only 0s and 1s.
	  \param ChannelOutput Object for output.
	  \sa NoisySeq
    */
	void GenerateOutput(int input[], NoisySeq &ChannelOutput);
	//void SetAWGNSTD(double std){mSTD = std;}
	//void SetBSCCrossoverProb(double p){mCrossoverProb = p;}
	//! Constructor.
	Channel(int blocklength){mBlocklength = blocklength;}
	//! destructor.
	~Channel()
	{
#ifdef DEBUG
	cout<<"~Channel()"<<endl;
#endif
	}

	//! Set fixed seed.
    /*!
      \param seed. Seed
    */
	void SetSeed(int seed);
};



//Channel Class
void Channel::SetChannel(int blocklength, string channel, double parameter)
{
	if(channel == "AWGN")
	{
		mSTD = parameter;
	}
	else if (channel == "BSC")
	{
		mCrossoverProb = parameter;
	}
	else
	{
		cout<<"channel not supported"<<endl;
		exit(0);
	}

	mBlocklength = blocklength;
	mChannelName = channel;
}
void Channel::GenerateOutput(int input[], NoisySeq &ChannelOutput)
{

	if(mChannelName == "AWGN")
	{
		ChannelOutput.ChannelParameter = mSTD;
		mGenerateAWGN(input, ChannelOutput);
	}
	else if (mChannelName == "BSC")
	{
		ChannelOutput.ChannelParameter = mCrossoverProb;
		mGenerateBitFlip(input, ChannelOutput);
	}
	else
	{
		cout<<"channel not supported"<<endl;
		exit(0);
	}
}
void Channel::mGenerateAWGN(int input[], NoisySeq &ChannelOutput)// generate noisy sequence for AWGN
{
	double std_check = 0;
	for (int i = 0; i < mBlocklength; i++)
	{
		int transmittedbit = 0;
		transmittedbit = 2 * input[i] - 1;

		double receivedbit = channelrand.randNorm(transmittedbit, mSTD); // assume all zero codeword is sent
		ChannelOutput.NoisySequence[i] = receivedbit;
		//_LogLikelihoodRatio[i] = ((receivedbit - 1)*(receivedbit - 1) - (receivedbit + 1)*(receivedbit + 1))/2.0/std/std;
		std_check += (receivedbit - transmittedbit)*(receivedbit - transmittedbit);
	}
#ifdef DEBUG
	cout<<"std = "<<sqrt(std_check/(double)mBlocklength)<<" real std = "<<mSTD<<endl;
#endif
}
void Channel::mGenerateBitFlip(int input[], NoisySeq &ChannelOutput)// generate noisy sequence for BSC
{
	int total_bit_flips = 0;
	for (int i = 0; i < mBlocklength; i++)
	{
		double rand_num = channelrand.rand();
		if(rand_num > mCrossoverProb)
		{
			ChannelOutput.NoisySequence[i] = input[i];
		}
		else
		{
			ChannelOutput.NoisySequence[i] = 1 - input[i];
			total_bit_flips++;
		}
	}
#ifdef DEBUG
	cout<<"p = "<<(double)total_bit_flips/mBlocklength<<" real p = "<<mCrossoverProb<<endl;
#endif
}

void Channel::SetSeed(int seed){
	channelrand.seed(seed);
}


#endif /* CHANNEL_HPP_ */
