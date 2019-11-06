#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <string>
#include <cstring>
#include <algorithm>

/*
Copyright 2019, Benjamin Coleman, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/

int main(int argc, char **argv){

    if (argc < 4){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"samplerace <tau> <input> <output> [--range race_range] [--reps race_reps] [--hashes n_minhashes] [-k kmer_size]"<<std::endl; 
        std::clog<<"Positional arguments: "<<std::endl; 
        std::clog<<"tau: floating point RACE sampling threshold. Roughly determines how many samples you will store"<<std::endl; 
        std::clog<<"input: path to input data file (.fastq or .fasta extension)"<<std::endl; 
        std::clog<<"output: path to output sample file (same extension as input)"<<std::endl; 
        std::clog<<"Optional arguments: "<<std::endl; 
        std::clog<<"[--range race_range]: (Optional, default 10000) Hash range for each ACE (B)"<<std::endl;
        std::clog<<"[--reps race_reps]: (Optional, default 10) Number of ACE repetitions (R)"<<std::endl;
        std::clog<<"[--hashes n_minhashes]: (Optional, default 1) Number of MinHashes for each ACE (n)"<<std::endl;
        std::clog<<"[--k kmer_size]: (Optional, default 16) Size of each MinHash k-mer (k)"<<std::endl;
        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<"samplerace 15.0 data/input.fastq data/output.fastq --range 100 --reps 50 --hashes 3 --k 5"<<std::endl; 
        std::clog<<"samplerace 0.1 data/input.fasta data/output.fasta --range 500 --k 20"<<std::endl; 
        return -1; 
    }

    double tau = std::stod(argv[1]);
    std::ifstream datastream(argv[2]);
    std::ofstream samplestream(argv[3]);

    int race_range = 1000;
    int race_repetitions = 10;
    int hash_power = 1;
    int kmer_k = 16;

    for (int i = 0; i < argc; ++i){
        if (std::strcmp("--range",argv[i]) == 0){
            if ((i+1) < argc){
                race_range = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --range"<<std::endl; 
                return -1;
            }
        }
        if (std::strcmp("--reps",argv[i]) == 0){
            if ((i+1) < argc){
                race_repetitions = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --reps"<<std::endl; 
                return -1;
            }
        }
        if (std::strcmp("--hashes",argv[i]) == 0){
            if ((i+1) < argc){
                hash_power = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --hashes"<<std::endl; 
                return -1;
            }
        }
        if (std::strcmp("--k",argv[i]) == 0){
            if ((i+1) < argc){
                kmer_k = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --k"<<std::endl; 
                return -1;
            }
        }
    }

    // determine file extension
    std::string filename(argv[2]); 
    std::string file_extension = "";
    size_t idx = filename.rfind('.',filename.length()); 
    if (idx != std::string::npos){
        file_extension = filename.substr(idx+1, filename.length() - idx); 
    } else {
        std::cerr<<"Input file does not appear to have any file extension."<<std::endl; 
        return -1; 
    }
    if (file_extension != "fasta" && file_extension != "fastq"){
        std::cerr<<"Unknown file extension: "<<file_extension<<std::endl; 
        std::cerr<<"Please specify either a .fasta or a .fastq file"<<std::endl; 
        return -1; 
    }

    // done parsing information. Begin RACE algorithm: 

    // buffer for sequences and fasta/fastq chunks
    std::string sequence;
    std::string chunk;

    // set up the hash function that will be used to hash input sequences
    SequenceMinHash hash = SequenceMinHash(race_repetitions*hash_power);
    int* raw_hashes = new int[race_repetitions*hash_power]; 
    int* rehashes = new int[race_repetitions];

    RACE sketch = RACE(race_repetitions,race_range); 

    do{
        bool success = SequenceFeatures(datastream, sequence, chunk, file_extension);
        if (!success) continue;

        hash.getHash(kmer_k, sequence, raw_hashes); 
        // now that we have the sequence and label
        // feed the sequence into the RACE structure
        // first rehash so that the arrays can fit into RACE
        rehash(raw_hashes, rehashes, race_repetitions, hash_power);
        // then simultaneously query and add 
        double KDE = sketch.query_and_add(rehashes); 
        // note: KDE is on a scale from [0,N] not the normalized interval [0,1]
        if (KDE < tau){
            // then keep this sample
            samplestream<<chunk;
        }
    }
    while(datastream);
}
