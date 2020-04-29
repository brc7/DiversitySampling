#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "MurmurHash.h"

bool SequenceFeatures(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat); 


bool SequenceFeaturesI(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat); 
bool SequenceFeaturesPE(std::istream& in1, std::istream& in2, std::string& sequence, std::string& chunk1, std::string& chunk2, std::string fastWhat); 
bool SequenceFeaturesSE(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat); 



