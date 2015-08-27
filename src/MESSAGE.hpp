/*
 * MESSAGE.hpp
 *
 *  Created on: 17 avr. 2015
 *      Author: legal
 */
#ifndef MESSAGE_HPP_
#define MESSAGE_HPP_

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

//!  A class template for message passing values
/*!
  The template thing is original intended to test error floors with float versus double.
*/
template <class DataType>
class MESSAGE {
 public:
	DataType Mprob; /*!< LLRs, can be of different data type, say float, double.  */
	int Mrow, Mcol; /*!< column, row index. */
	int degree;     /*!< degree of the node */
	int *i;         /*!< store nodes (in the bipartite graph) that are linked to this node */

	//! Set Degree, initialize the array
    /*!
      \param deg the degree of the current node
    */
	void SetDeg(int deg)
	{
		degree = deg;
		i      = new int[deg];
	}
	//! A constructor.
	MESSAGE()
	{
		i = NULL;
	}
	//! A destructor.
	~MESSAGE()
	{
		delete [] i;
	}
};

#endif /* MESSAGE_HPP_ */
