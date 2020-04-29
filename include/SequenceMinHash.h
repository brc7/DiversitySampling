#pragma once 

#include <limits>
#include <vector>
#include <string>

#include <iostream>

#include "MurmurHash.h"


class SequenceMinHash {
private: 
	int _numhashes; 
public: 
	SequenceMinHash(int number_of_hashes); 
	void getHash(size_t k, const std::string& sequence, int* hashes); 
	unsigned int internalHash(int input, int seed); 
}; 



