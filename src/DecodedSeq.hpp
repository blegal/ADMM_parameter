/*
 * DecodedSeq.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */

#ifndef DECODEDSEQ_HPP_
#define DECODEDSEQ_HPP_

//!  Class for decoded sequence.
/*!
  Stores the sequence from decoder, also contains information about decoder status.
*/
class DecodedSeq{
public:
	int Blocklength; /*!< blocklength of code */

	double *SoftInfo; /*!< soft information */
	int *HardDecision; /*!< 0/1 decision */
	bool AlgorithmConverge; /*!< true if the decoder converged */
	int Iteration; /*!< number of iteration used */
	bool ValidCodeword; /*!< true if the codeword is a valid codeword */
	bool LPIntegerSol; /*!< true if the codeword is an integral solution from ADMM-LP decoder */
	unsigned long int ExeTime; /*!< decoding time */
	bool IsNan, IsInf;
	//! Constructor.
    /*!
      \param blocklength Blocklength of code
	  \param channel Channel name, can only be "AWGN" or "BSC"
    */
	DecodedSeq(int blocklength)
	{
		Blocklength  = blocklength;
		HardDecision = new int[blocklength];
		SoftInfo     = new double[blocklength];
		IsNan        = false;
		IsInf        = false;
	}

	//! destructor.
	~DecodedSeq()
	{
#ifdef DEBUG
	cout<<"~DecodedSeq()"<<endl;
#endif
		delete [] HardDecision;
		delete [] SoftInfo;
	}
};




#endif /* DECODEDSEQ_HPP_ */
