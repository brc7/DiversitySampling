#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "MurmurHash.h"

bool SequenceFeatures(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat); 

