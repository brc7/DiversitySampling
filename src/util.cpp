#include "util.h"

/*
Copyright 2019, Benjamin Coleman, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/


void rehash(int* input_hashes, int* output_hashes, int nhashes, int values_per_set){
    #pragma omp parallel for 
    for (size_t i = 0; i < nhashes; i++){
       output_hashes[i] = MurmurHash(input_hashes + values_per_set*i, sizeof(int)*values_per_set, 42);
    }
}