#include "io.h"
#include "SequenceMinHash.h"
#include "RACE.h"
#include "util.h"

#include <chrono>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>

/*
Copyright 2019, Benjamin Coleman, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/


/*
Three types of reads: paired, interleaved and single

For single reads just do normally 
For interleaved reads just do normally but save "chunks"
For paired reads, assume they're in order (if they're not,you can use fastq-pair) 
and "rescue" saved reads


Desired interface:
samplerace tau SE input output <flags>
samplerace tau PE input1 input2 output1 output2 <flags>
samplerace tau I input output <flags>

*/

int main(int argc, char **argv){

    if (argc < 4){
        std::clog<<"Usage: "<<std::endl; 
        std::clog<<"samplerace <tau> <format> <input> <output>";
        std::clog<<" [--range race_range] [--reps race_reps] [--hashes n_minhashes] [-k kmer_size]"<<std::endl; 
        std::clog<<"Positional arguments: "<<std::endl; 
        std::clog<<"tau: csv list of one or more floating point RACE sampling thresholds. Roughly determines how many samples you will store. You may specify this in scientific notation (i.e. 10e-6). Note: Don't include spaces between the elements."<<std::endl;
        std::clog<<"format: Either PE, SE, or I for paired-end, single-end, and interleaved paired reads"<<std::endl;
		std::clog<<"input: path to text file containing working directory of read files followed by names of the read files to be sampled (.fastq or .fasta extension). For PE format, specify two filenames, one containing \"_1\" and the other \"_2\", for each read."<<std::endl;
        std::clog<<"output: path to output sample file (same extension as input). For PE format, specify two files."<<std::endl; 
        
        std::clog<<"Optional arguments: "<<std::endl; 
        std::clog<<"[--range race_range]: (Optional, default 10000) Hash range for each ACE (B)"<<std::endl;
        std::clog<<"[--reps race_reps]: (Optional, default 10) Number of ACE repetitions (R)"<<std::endl;
        std::clog<<"[--hashes n_minhashes]: (Optional, default 1) Number of MinHashes for each ACE (n)"<<std::endl;
        std::clog<<"[--k kmer_size]: (Optional, default 16) Size of each MinHash k-mer (k)"<<std::endl;

        std::clog<<std::endl<<"Example usage:"<<std::endl; 
        std::clog<<"samplerace 15.0 PE data/input-1.fastq data/input-2.fastq data/output-1.fastq data/output-2.fastq --range 100 --reps 50 --hashes 3 --k 5"<<std::endl; 
        std::clog<<"samplerace 10e-6 SE data/input.fastq data/output.fastq --range 100 --reps 5 --hashes 1 --k 33"<<std::endl; 
        std::clog<<"samplerace 0.1 SE data/input.fasta data/output.fasta --range 100000 --k 20"<<std::endl; 
        return -1; 
    }


    // POSITIONAL ARGUMENTS

	// get set of taus
	std::vector<double> taus;
	std::stringstream ss(argv[1]); // get set of taus
	double tau = 0;
	while(ss >> tau){
		if (tau <= 0){ std::cerr<<"Invalid value for parameter <tau>"<<std::endl; return -1; }
		taus.push_back(tau);
		if (ss.peek() == ',')
			ss.ignore();
	}
	ss.clear();

    int format; // ENUM: 1 = unpaired, 2 = interleaved, 3 = paired
    if (std::strcmp("SE",argv[2]) == 0){
        format = 1;
    } else if (std::strcmp("I",argv[2]) == 0){
        format = 2; 
    } else if (std::strcmp("PE",argv[2]) == 0){
        format = 3; 
        if (argc < 6){
            std::cerr<<"For paired-end reads, please specify the output files as:"<<std::endl;
            std::cerr<<"output1.fastq output2.fastq"<<std::endl;
            return -1; 
        }
    } else {
        std::cerr<<"Invalid format, please specify either SE, PE, or I"<<std::endl; 
        return -1;
    }

    std::vector<std::string> filelist1;
	std::vector<std::string> filelist2;

    // open the correct file streams given the format
    std::ifstream datastream1;
    std::ifstream datastream2;

	// Create vector of ofstreams - one for each tau
	std::vector<std::ofstream> samplestreamvector1;
	std::vector<std::ofstream> samplestreamvector2;

    if (format != 3){

    	// add input read filenames to the filelist vector.
    	datastream1.open(argv[3]);
		std::string line;
    	std::string path = ""; // assume read files are in the same directory by default.
		do {
			int c = datastream1.peek();
			if (c == EOF) {
				if (datastream1.eof()){
					break;
				}
			}
			std::getline(datastream1, line);
			if (line.at(0) == '/' && line.rfind('.',line.length()) == std::string::npos) {path = line;}
			else if (line.rfind('.',line.length()) != std::string::npos) {filelist1.push_back(path+"/"+line);}
			else {std::cerr<<"Line \""<<line<<"\" not recognized. If it is a filename, it has to end with a file extension. If it is a path it has to start with a \"/\"."<<std::endl;}
		} while (datastream1);

		// add output file streams
		std::string baseoutputfilename(argv[4]);
		for(size_t i = 0; i < taus.size(); i++){
			std::string filename = baseoutputfilename;
			filename += "-";
			filename += std::to_string(taus[i]);
			filename += ".fastq"; // just realized this should depend on input type
			std::ofstream s(filename);
			samplestreamvector1.push_back(std::move(s));
		}
    } else {

		// add input read filenames to the filelist vector.
		datastream1.open(argv[3]);
		std::string line;
		std::string path = ""; // assume read files are in the same directory by default.
		do {
			int c = datastream1.peek();
			if (c == EOF) {
				if (datastream1.eof()){
					break;
				}
			}
			std::getline(datastream1, line);
			if (line.at(0) == '/' && line.rfind('.',line.length()) == std::string::npos) {
				path = line;
			} else if (line.rfind("_1",line.length()) != std::string::npos && line.rfind('.',line.length()) != std::string::npos) {
				filelist1.push_back(path+"/"+line);
			} else if (line.rfind("_2",line.length()) != std::string::npos && line.rfind('.',line.length()) != std::string::npos) {
				filelist2.push_back(path+"/"+line);
			} else {
				std::cerr<<"Line \""<<line<<"\" not recognized. If it is a filename, it has to end with a file extension and either \"_1\" or \"_2\" for paired end reads. If it is a path it has to start with a \"/\"."<<std::endl;
			}
		} while (datastream1);

		// add output file streams
		std::string baseoutputfilename1(argv[4]);
		std::string baseoutputfilename2(argv[5]);
		for(size_t i = 0; i < taus.size(); i++){
			std::string filename1 = baseoutputfilename1;
			std::string filename2 = baseoutputfilename2;
			filename1 += "-";
			filename2 += "-";
			filename1 += std::to_string(taus[i]);
			filename2 += std::to_string(taus[i]);
			filename1 += ".fastq"; // just realized this should depend on input type
			filename2 += ".fastq";
			std::ofstream s1(filename1);
			std::ofstream s2(filename2);
			samplestreamvector1.push_back(std::move(s1));
			samplestreamvector2.push_back(std::move(s2));
		}
    }

    // OPTIONAL ARGUMENTS
    int race_range = 10000;
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

    // Check if arguments are valid
    // Checking for taus is done when adding taus to the vector.
    if (race_range <= 0){ std::cerr<<"Invalid value for optional parameter --range"<<std::endl; return -1; }
    if (race_repetitions <= 0){ std::cerr<<"Invalid value for optional parameter --reps"<<std::endl; return -1; }
    if (hash_power <= 0){ std::cerr<<"Invalid value for optional parameter --hashes"<<std::endl; return -1; }
    if (kmer_k <= 0){ std::cerr<<"Invalid value for optional parameter --k"<<std::endl; return -1; }

    // done parsing information. Begin RACE algorithm: 

    // buffer for sequences and fasta/fastq chunks
    std::string sequence;
    std::string chunk1;
    std::string chunk2;

    // set up the hash function that will be used to hash input sequences
    SequenceMinHash hash = SequenceMinHash(race_repetitions*hash_power);
    int* raw_hashes = new int[race_repetitions*hash_power]; 
    int* rehashes = new int[race_repetitions];

    RACE sketch = RACE(race_repetitions,race_range); 

    for (int i = 0; i < filelist1.size(); ++i) {
		// determine file extension. Only checks for first file for paired-end reads.
		std::string filename1(filelist1[i]);
		std::string file_extension = "";
		size_t idx = filename1.rfind('.',filename1.length());
		if (idx != std::string::npos){
			file_extension = filename1.substr(idx+1, filename1.length() - idx);
		} else {
			std::cerr<<"Input file " << filename1 << " does not appear to have any file extension."<<std::endl;
			continue;
		}
		if (file_extension == "fq"){
			file_extension = "fastq";
		}
		if (file_extension != "fasta" && file_extension != "fastq"){
			std::cerr<<"Unknown file extension: "<<file_extension<<std::endl;
			std::cerr<<"Please specify either a file with the .fasta or .fastq extension."<<std::endl;
			continue;
		}

		// Reset and open appropriate datastreams according to format
		switch(format){
			case 1: // 1 = unpaired
				datastream1.close();
				datastream1.clear();
				datastream1.open(filelist1[i]);
				std::cout<<"Processing "<<filelist1[i]<<std::endl;
				break;
			case 2: // 2 = interleaved
				datastream1.close();
				datastream1.clear();
				datastream1.open(filelist1[i]);
				std::cout<<"Processing "<<filelist1[i]<<std::endl;
				break;
			case 3: // 3 = paired
				datastream1.close();
				datastream1.clear();
				datastream1.open(filelist1[i]);
				datastream2.close();
				datastream2.clear();
				datastream2.open(filelist2[i]);
				std::cout<<"Processing "<<filelist1[i]<<" and "<<filelist2[i]<<std::endl;
				break;
		}

		do{
			bool success = false;
			int c = datastream1.peek();
			if (c == EOF) {
				if (datastream1.eof()){
					break;
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
			double KDE = sketch.query_and_add(rehashes);
			// note: KDE is on a scale from [0,N] not the normalized interval [0,1]
			for(size_t i = 0; i < taus.size(); i++){
				double tau = taus[i];
				if (KDE < tau){
					switch(format){
						case 1: // 1 = unpaired
							samplestreamvector1[i] << chunk1;
							break;
						case 2: // 2 = interleaved
							samplestreamvector1[i] << chunk1;
							break;
						case 3: // 3 = paired
							samplestreamvector1[i] << chunk1;
							samplestreamvector2[i] << chunk2;
							break;
					}
				}
			}
		}
		while(datastream1);
    }
}
