#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <numeric>

#include <chrono>

// #define DEBUG

/*
Copyright 2022, Benjamin Coleman, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author
*/


/*
Three types of reads: paired, interleaved and single

For single reads just do normally 
For interleaved reads just do normally but save "chunks"
For paired reads, assume they're in order (if they're not,you can use fastq-pair) 
and "rescue" saved reads
*/


typedef struct read_info_SEI {
    std::streampos pos1;
    float score;
} read_info_SEI;

typedef struct read_info_PE {
    std::streampos pos1;
    std::streampos pos2;
    float score;
} read_info_PE;


int main(int argc, char **argv){

    if (argc < 4){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"permuterace <format> <input> <output>"; 
        std::clog<<" [--range race_range] [--reps race_reps] [--hashes n_minhashes] [-k kmer_size]"<<std::endl; 
        std::clog<<"Positional arguments: "<<std::endl; 
        std::clog<<"format: Either PE, SE, or I for paired-end, single-end, and interleaved paired reads"<<std::endl; 
        std::clog<<"input: path to input data file (.fastq or .fasta extension). For PE format, specify two files."<<std::endl; 
        std::clog<<"output: path to output sample file (same extension as input). For PE format, specify two files."<<std::endl; 
        
        std::clog<<"Optional arguments: "<<std::endl; 
        std::clog<<"[--range race_range]: (Optional, default 10000) Hash range for each ACE (B)"<<std::endl;
        std::clog<<"[--reps race_reps]: (Optional, default 10) Number of ACE repetitions (R)"<<std::endl;
        std::clog<<"[--hashes n_minhashes]: (Optional, default 1) Number of MinHashes for each ACE (n)"<<std::endl;
        std::clog<<"[--k kmer_size]: (Optional, default 16) Size of each MinHash k-mer (k)"<<std::endl;
        std::clog<<"[--chunksize chunk_size]: (Optional, default 100000) Number of reads to write to disk at once during permutations. Larger = faster but more RAM."<<std::endl;
        std::clog<<"[--scoretype score_type]: (Optional, default R) Type of score to use. Either R (for running KDE score), N (for normalized KDE score) or F (for full KDE score). Note that F will require ~50% more computation time."<<std::endl;

        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<"permuterace PE data/input-1.fastq data/input-2.fastq data/output-1.fastq data/output-2.fastq --range 100 --reps 50 --hashes 3 --k 5"<<std::endl; 
        std::clog<<"permuterace SE data/input.fastq data/output.fastq --range 100 --reps 5 --hashes 1 --k 33"<<std::endl; 
        std::clog<<"permuterace SE data/input.fasta data/output.fasta --range 100000 --k 20"<<std::endl; 
        return -1; 
    }


    // POSITIONAL ARGUMENTS
    int format; // ENUM: 1 = unpaired, 2 = interleaved, 3 = paired
    if (std::strcmp("SE",argv[1]) == 0){
        format = 1;
    } else if (std::strcmp("I",argv[1]) == 0){
        format = 2; 
    } else if (std::strcmp("PE",argv[1]) == 0){
        format = 3; 
        if (argc < 6){
            std::cerr<<"For paired-end reads, please specify the input and output files as:"<<std::endl; 
            std::cerr<<"input1.fastq input2.fastq output1.fastq output2.fastq"<<std::endl; 
            return -1; 
        }
    } else {
        std::cerr<<"Invalid format, please specify either SE, PE, or I"<<std::endl; 
        return -1;
    }

    // open the correct file streams given the format
    std::ifstream datastream1;
    std::ofstream samplestream1;
    std::ifstream datastream2;
    std::ofstream samplestream2;

    if (format != 3){
        datastream1.open(argv[2]);
        samplestream1.open(argv[3]);
    } else {
        datastream1.open(argv[2]);
        datastream2.open(argv[3]);
        samplestream1.open(argv[4]);
        samplestream2.open(argv[5]);
    }

    // determine file extension
    std::string filename(argv[2]); 
    std::string file_extension = "";
    size_t idx = filename.rfind('.',filename.length()); 
    if (file_extension == "fq"){
        file_extension = "fastq"; 
    }
    if (idx != std::string::npos){
        file_extension = filename.substr(idx+1, filename.length() - idx); 
    } else {
        std::cerr<<"Input file does not appear to have any file extension."<<std::endl; 
        return -1; 
    }
    if (file_extension != "fasta" && file_extension != "fastq"){
        std::cerr<<"Unknown file extension: "<<file_extension<<std::endl; 
        std::cerr<<"Please specify either a file with the .fasta or .fastq extension."<<std::endl; 
        return -1; 
    }

    // OPTIONAL ARGUMENTS
    int race_range = 10000;
    int race_repetitions = 10;
    int hash_power = 1;
    int kmer_k = 16;
    int batch_size = 100000;
    int score_type = 1; // ENUM: 1 = running, 2 = normalized, 3 = full

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
        if (std::strcmp("--chunksize",argv[i]) == 0){
            if ((i+1) < argc){
                batch_size = std::stoi(argv[i+1]);
            } else {
                std::cerr<<"Invalid argument for optional parameter --chunksize"<<std::endl; 
                return -1;
            }
        }
        if (std::strcmp("--scoretype",argv[i]) == 0){
            if ((i+1) < argc){
                if (std::strcmp("R",argv[i+1]) == 0){
                    score_type = 1;
                } else if (std::strcmp("N",argv[i+1]) == 0){
                    score_type = 2; 
                } else if (std::strcmp("F",argv[i+1]) == 0){
                    score_type = 3; 
                } else {
                    std::cerr<<"Invalid scoretype, please specify either R, N or F"<<std::endl; 
                    return -1;
                }
            } else {
                std::cerr<<"Invalid argument for optional parameter --scoretype"<<std::endl; 
                return -1;
            }
        }
    }

    // Check if arguments are valid
    if (race_range <= 0){ std::cerr<<"Invalid value for optional parameter --range"<<std::endl; return -1; }
    if (race_repetitions <= 0){ std::cerr<<"Invalid value for optional parameter --reps"<<std::endl; return -1; }
    if (hash_power <= 0){ std::cerr<<"Invalid value for optional parameter --hashes"<<std::endl; return -1; }
    if (kmer_k <= 0){ std::cerr<<"Invalid value for optional parameter --k"<<std::endl; return -1; }
    if (batch_size <= 0){ std::cerr<<"Invalid value for optional parameter --chunksize"<<std::endl; return -1; }

    // done parsing information. Begin permutation algorithm: 
    // buffer for sequences and fasta/fastq chunks
    std::string sequence;
    std::string chunk1;
    std::string chunk2;

    // set up the hash function that will be used to hash input sequences
    SequenceMinHash hash = SequenceMinHash(race_repetitions*hash_power);
    int* raw_hashes = new int[race_repetitions*hash_power]; 
    int* rehashes = new int[race_repetitions];

    RACE sketch = RACE(race_repetitions,race_range); 

    // vector of scores
    std::vector<read_info_PE> read_index_PE;
    std::vector<read_info_SEI> read_index_SEI;

    auto score_start = std::chrono::high_resolution_clock::now();

    do{
        bool success = false; 
        int c = datastream1.peek(); 
        if (c == EOF) {
            if (datastream1.eof()){
                continue; 
            }
        }

        std::streampos pos1 = datastream1.tellg();
        std::streampos pos2 = datastream2.tellg();

        switch(format){
            case 1: // 1 = unpaired
            success = SequenceFeaturesSE(datastream1, sequence, chunk1, file_extension);
            break; 
            case 2: // 2 = interleaved
            success = SequenceFeaturesI(datastream1, sequence, chunk1, file_extension); 
            break; 
            case 3: // 3 = paired
            success = SequenceFeaturesPE(datastream1, datastream2, sequence, chunk1, chunk2, file_extension);
            break; 
        }
        if (!success) continue;

        hash.getHash(kmer_k, sequence, raw_hashes); 
        // now that we have the sequence and label
        // feed the sequence into the RACE structure
        // first rehash so that the arrays can fit into RACE
        rehash(raw_hashes, rehashes, race_repetitions, hash_power);
        // then simultaneously query and add 
        double KDE = sketch.query_and_add(rehashes); 
        // note: KDE is on a scale from [0,N] not the normalized interval [0,1]
        switch(format){
            case 1: // 1 = unpaired
            case 2: // 2 = interleaved
            read_index_SEI.push_back({pos1, (float)(KDE)});
            break;
            case 3: // 3 = paired
            read_index_PE.push_back({pos1, pos2, (float)(KDE)});
            break;
        }
    }
    while(datastream1);

    // handle alternative score types
    if (score_type == 2){
        // normalize scores by the number of elements seen so far
        for(size_t data_id = 0; data_id < read_index_SEI.size(); data_id++){
            read_index_SEI[data_id].score = read_index_SEI[data_id].score / (data_id + 1);
        }
        for(size_t data_id = 0; data_id < read_index_PE.size(); data_id++){
            read_index_PE[data_id].score = read_index_PE[data_id].score / (data_id + 1);
        }
    } else if (score_type == 3){
        // re-compute scores using the total KDE
        datastream1.clear();
        datastream1.seekg(0);
        datastream2.clear();
        datastream2.seekg(0);
        size_t data_id = 0;
        do{
            bool success = false; 
            int c = datastream1.peek(); 
            if (c == EOF) {
                if (datastream1.eof()){
                    continue; 
                }
            }
            switch(format){
                case 1: // 1 = unpaired
                success = SequenceFeaturesSE(datastream1, sequence, chunk1, file_extension);
                break; 
                case 2: // 2 = interleaved
                success = SequenceFeaturesI(datastream1, sequence, chunk1, file_extension); 
                break; 
                case 3: // 3 = paired
                success = SequenceFeaturesPE(datastream1, datastream2, sequence, chunk1, chunk2, file_extension);
                break; 
            }
            if (!success) continue;

            hash.getHash(kmer_k, sequence, raw_hashes); 
            // now that we have the sequence and label
            // feed the sequence into the RACE structure
            // first rehash so that the arrays can fit into RACE
            rehash(raw_hashes, rehashes, race_repetitions, hash_power);
            // then simultaneously query and add 
            double KDE = sketch.query(rehashes); 
            // note: KDE is on a scale from [0,N] not the normalized interval [0,1]
            switch(format){
                case 1: // 1 = unpaired
                case 2: // 2 = interleaved
                read_index_SEI[data_id].score = (float)(KDE);
                break;
                case 3: // 3 = paired
                read_index_PE[data_id].score = (float)(KDE);
                break;
            }
            data_id++;
        }
        while(datastream1);
    }

    // sort reads by score to get the permutation function
    std::sort(read_index_SEI.begin(), read_index_SEI.end(),
        [](const read_info_SEI &left, const read_info_SEI &right){
            return left.score < right.score;});
    std::sort(read_index_PE.begin(), read_index_PE.end(),
        [](const read_info_PE &left, const read_info_PE &right){
            return left.score < right.score;});

    std::chrono::duration<double, std::milli> score_duration = std::chrono::high_resolution_clock::now() - score_start;
    std::cout<<"Scoring took "<<score_duration.count()/1000.0<<" seconds."<<std::endl;    

    // collect the reads into properly-ordered chunks and write each block to the output file
    std::vector<size_t> sorted_idx(batch_size); // for sorting by index
    std::vector<std::string> batch_reads1(batch_size);
    std::vector<std::string> batch_reads2(batch_size);
    size_t num_reads = std::max(read_index_PE.size(), read_index_SEI.size());

    auto write_start = std::chrono::high_resolution_clock::now();

    for (size_t num_written = 0; num_written < num_reads; num_written += batch_size){
        // auto batch_start = std::chrono::high_resolution_clock::now();
        size_t num_available = (num_reads - num_written < batch_size) ? num_reads - num_written : batch_size;
        switch(format){
            case 1:
                // sort the reads in order of their appearance in the file.
                // sorted_idx[n] is the location inside the batch of the
                // nth read to appear in the datastream.
                #ifdef DEBUG
                std::cout<<"Sorting batch of "<<num_available<<" from "<<num_written<<" to "<<num_written+num_available<<std::endl;
                #endif
                std::iota(sorted_idx.begin(), sorted_idx.begin() + num_available, num_written);
                // This lambda is a bit tricky.
                std::sort(sorted_idx.begin(), sorted_idx.begin() + num_available,
                    [&read_index_SEI](size_t idx_left, size_t idx_right){
                        return read_index_SEI[idx_left].pos1 < read_index_SEI[idx_right].pos1;});
                #ifdef DEBUG
                std::cout<<"Batch of items in input-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = sorted_idx[n];
                    std::cout<<" ("<<i<<" ["<<read_index_SEI[i].score<<"] @ "<<read_index_SEI[i].pos1<<") ";
                }
                std::cout<<std::endl;

                std::cout<<"Batch of items in output-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = num_written + n;
                    std::cout<<" ("<<i<<"["<<read_index_SEI[i].score<<"] @ "<<read_index_SEI[i].pos1<<") ";
                }
                std::cout<<std::endl;
                #endif

                for (size_t n = 0; n < num_available; n++){
                    // n is the order of occurence in the input file
                    // i is the order of occurence in the output file
                    size_t i = sorted_idx[n];
                    size_t batch_idx = i - num_written;
                    datastream1.clear();
                    datastream1.seekg(read_index_SEI[i].pos1);
                    // this is messy, we should handle the multiple formats more elegantly
                    bool success = SequenceFeaturesSE(datastream1, sequence, chunk1, file_extension);
                    batch_reads1[batch_idx] = chunk1;
                    #ifdef DEBUG
                    std::cout<<"\tRead "<<i<<" from position "<<read_index_SEI[i].pos1<<" into batch index "<<batch_idx<<"("<<n<<"/"<<num_available<<")"<<std::endl;
                    #endif
                }
                for (size_t n = 0; n < num_available; n++){
                    samplestream1 << batch_reads1[n];
                }
            break;
            case 2:
                // sort the reads in order of their appearance in the file.
                // sorted_idx[n] is the location inside the batch of the
                // nth read to appear in the datastream.
                #ifdef DEBUG
                std::cout<<"Sorting batch of "<<num_available<<" from "<<num_written<<" to "<<num_written+num_available<<std::endl;
                #endif
                std::iota(sorted_idx.begin(), sorted_idx.begin() + num_available, num_written);
                // This lambda is a bit tricky.
                std::sort(sorted_idx.begin(), sorted_idx.begin() + num_available,
                    [&read_index_SEI](size_t idx_left, size_t idx_right){
                        return read_index_SEI[idx_left].pos1 < read_index_SEI[idx_right].pos1;});

                #ifdef DEBUG
                std::cout<<"Batch of items in input-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = sorted_idx[n];
                    std::cout<<" ("<<i<<" ["<<read_index_SEI[i].score<<"] @ "<<read_index_SEI[i].pos1<<") ";
                }
                std::cout<<std::endl;

                std::cout<<"Batch of items in output-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = num_written + n;
                    std::cout<<" ("<<i<<"["<<read_index_SEI[i].score<<"] @ "<<read_index_SEI[i].pos1<<") ";
                }
                std::cout<<std::endl;
                #endif

                for (size_t n = 0; n < num_available; n++){
                    // n is the order of occurence in the input file
                    // i is the order of occurence in the output file
                    size_t i = sorted_idx[n];
                    size_t batch_idx = i - num_written;
                    datastream1.clear();
                    datastream1.seekg(read_index_SEI[i].pos1);
                    // this is messy, we should handle the multiple formats more elegantly
                    bool success = SequenceFeaturesSE(datastream1, sequence, chunk1, file_extension);
                    batch_reads1[batch_idx] = chunk1;
                    #ifdef DEBUG
                    std::cout<<"\tRead "<<i<<" from position "<<read_index_SEI[i].pos1<<" into batch index "<<batch_idx<<"("<<n<<"/"<<num_available<<")"<<std::endl;
                    #endif
                }
                for (size_t n = 0; n < num_available; n++){
                    samplestream1 << batch_reads1[n];
                }
            break;
            case 3:
                // sort the reads in order of their appearance in the file.
                // sorted_idx[n] is the location inside the batch of the
                // nth read to appear in the datastream.
                #ifdef DEBUG
                std::cout<<"Sorting batch of "<<num_available<<" from "<<num_written<<" to "<<num_written+num_available<<std::endl;
                #endif
                std::iota(sorted_idx.begin(), sorted_idx.begin() + num_available, num_written);
                // This lambda is a bit tricky.
                std::sort(sorted_idx.begin(), sorted_idx.begin() + num_available,
                    [&read_index_PE](size_t idx_left, size_t idx_right){
                        return read_index_PE[idx_left].pos1 < read_index_PE[idx_right].pos1;});

                #ifdef DEBUG
                std::cout<<"Batch of items in input-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = sorted_idx[n];
                    std::cout<<" ("<<i<<" ["<<read_index_PE[i].score<<"] @ "<<read_index_PE[i].pos1<<","<<read_index_PE[i].pos2<<") ";
                }
                std::cout<<std::endl;

                std::cout<<"Batch of items in output-file order:"<<std::endl;
                for (size_t n = 0; n < num_available; n++){
                    size_t i = num_written + n;
                    std::cout<<" ("<<i<<"["<<read_index_PE[i].score<<"] @ "<<read_index_PE[i].pos1<<","<<read_index_PE[i].pos2<<") ";
                }
                std::cout<<std::endl;
                #endif

                for (size_t n = 0; n < num_available; n++){
                    // n is the order of occurence in the input file
                    // i is the order of occurence in the output file
                    size_t i = sorted_idx[n];
                    size_t batch_idx = i - num_written;
                    datastream1.clear();
                    datastream1.seekg(read_index_PE[i].pos1);
                    datastream2.clear();
                    datastream2.seekg(read_index_PE[i].pos2);
                    // this is messy, we should handle the multiple formats more elegantly
                    bool success = success = SequenceFeaturesPE(datastream1, datastream2, sequence, chunk1, chunk2, file_extension);
                    batch_reads1[batch_idx] = chunk1;
                    batch_reads2[batch_idx] = chunk2;
                    #ifdef DEBUG
                    std::cout<<"\tRead "<<i<<" from positions "<<read_index_PE[i].pos1<<","<<read_index_PE[i].pos2<<" into batch index "<<batch_idx<<"("<<n<<"/"<<num_available<<")"<<std::endl;
                    #endif
                }
                for (size_t n = 0; n < num_available; n++){
                    samplestream1 << batch_reads1[n];
                    samplestream2 << batch_reads2[n];
                }
            break;
        }
        // std::chrono::duration<double, std::milli> duration = std::chrono::high_resolution_clock::now() - batch_start;
        // std::cout<<"Batch write time: "<<duration.count()<<" ms. Time/read: "<<duration.count() / (float)(num_available)<<" ms."<<std::endl;
    }
    std::chrono::duration<double, std::milli> duration = std::chrono::high_resolution_clock::now() - write_start;
    std::cout<<"Writing took "<<duration.count()/1000.0<<" seconds."<<std::endl;
}
