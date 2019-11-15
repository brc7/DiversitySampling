#include "io.h"

/*
Copyright 2019, Benjamin Coleman, All rights reserved. 
Free for research use. For commercial use, contact 
Rice University Invention & Patent or the author

*/

bool SequenceFeatures(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat){

    int chunksize = 0; 
    char begin;
    int linesread = 0; 
    if (fastWhat == "fasta") {
        chunksize = 2;
        begin = '>';
    } else if (fastWhat == "fastq") {
        chunksize = 4;
        begin = '@';
    } else {
        std::cerr<<"Unsupported file type: "<<fastWhat<<std::endl; 
        return false; // unsupported file type
    }

    sequence = "";
    chunk = ""; 

    // read meta data
    std::string buffer;
    std::getline(in, buffer);
    linesread++;
    if (buffer.length() == 0) {
        std::cerr<<"Found empty line - probably reached the end of the input "<<fastWhat<<" file."<<std::endl; 
        std::cerr<<buffer<<std::endl; 
        return false; 
    }
    if (buffer.at(0) != begin){
        std::cerr<<"Error reading line of "<<fastWhat<<" file: Expected a line beginning with "<<begin<<" but instead found: "<<std::endl; 
        std::cerr<<buffer<<std::endl; 
        return false; 
    }

    chunk += buffer; 
    chunk += "\n"; 

    std::getline(in, sequence); 
    linesread++; 
    if (sequence.length() == 0){
        std::cerr<<"Error reading "<<fastWhat<<" file: Expected a sequence, but found empty line for read ID: "<<std::endl; 
        std::cerr<<chunk<<std::endl; 
        return false; 
    }
    chunk += sequence; 
    chunk += "\n"; 

    while((linesread < chunksize) && (in)){
        std::getline(in, buffer); 
        linesread++; 
        chunk += buffer; 
        chunk += "\n"; 
    }

    if (linesread != chunksize){
        std::cerr<<"Error reading "<<fastWhat<<" file: Expected "<<chunksize<<" lines of input, got "<<linesread<<" instead. The sequence data should still be valid. "<<std::endl; 
        return false; 
    }
    return true;
}

bool SequenceFeaturesSE(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat){
    return SequenceFeatures(in, sequence, chunk, fastWhat); 
}

bool SequenceFeaturesPE(std::istream& in1, std::istream& in2, std::string& sequence, std::string& chunk1, std::string& chunk2, std::string fastWhat){
    bool success1 = SequenceFeatures(in2, sequence, chunk2, fastWhat); 
    bool success2 = SequenceFeatures(in1, sequence, chunk1, fastWhat); 

    if (!success1){
        std::cerr<<"Error reading first "<<fastWhat<<" file"<<std::endl; 
        return false; 
    }
    if (!success2){
        std::cerr<<"Error reading second "<<fastWhat<<" file"<<std::endl; 
        return false; 
    }
    return (success1 && success2); 

}

bool SequenceFeaturesI(std::istream& in, std::string& sequence, std::string& chunk, std::string fastWhat){
    bool success = SequenceFeatures(in, sequence, chunk, fastWhat); 

    // now tack on another few lines for the second chunk of the interleaved file! 
    int chunksize = 0; 
    int linesread = 0; 
    if (fastWhat == "fasta") {
        chunksize = 2;
    } else if (fastWhat == "fastq") {
        chunksize = 4;
    }
    std::string buffer;
    while((linesread < chunksize) && (in)){
        std::getline(in, buffer); 
        linesread++; 
        chunk += buffer; 
        chunk += "\n"; 
    }
    return success; 
}


