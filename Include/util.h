#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <math.h>
#include <algorithm>

#include "MurmurHash.h"

// mathematic and hash utilities, collected in one place
void rehash(int* input_hashes, int* output_hashes, int nhashes, int values_per_set); 

