# DiversitySampling


## Introduction
This repository implements a one-pass diversity sampling algorithm for DNA sequences. If you give this tool a stream of sequences, it automatically saves a diverse subset of them. Sequences are saved based on their distance from sequences that are already represented in the sample. It would be very slow to explicitly calculate the distances, so we use locality-sensitive hashing (LSH) to implement the fast algorithm described in our paper. 

This is interesting if you want to quickly estimate the biodiversity of a metagenomic sample, since different organisms often have sufficiently different sequences that our tool can ensure that each organism is represented in the sample. However, the tool should also work well for downsampling from other kinds of categories (chromosomes, genera, classes of bacteria, etc), provided the sequences are different enough. In a nutshell, our method attempts to save a set of sequences that maximizes genetic diversity. 

## Installation
The algorithm is implemented as a command line tool written in C++. It takes in a fasta or fastq file, downsamples the file using the RACE method, and outputs a fasta or fastq file of diverse sequences. To compile the tool, clone the repo and 
```
make binaries
```

The Makefile should produce build and bin directories and output the executable file samplerace to bin/. This should work fine on most Linux systems. If something goes wrong, it is probably because your C++ compiler does not support C++11 or OpenMP. In particular, on MacOS the g++ command aliases to an outdated version of clang that does not support the -fopenmp flag. You will need to change the CXX parameter in the Makefile to point to a compiler that implements OpenMP. Alternatively, you can remove the -fopenmp flag, but the hash computations will be slower. Windows does not include g++ by default, so you will want to install a compiler with OpenMP and C++11 support. 

## Algorithm and Hyperparameters
We use the RACE data structure, which is an efficient way to estimate kernel densities on streaming data. RACE is a small 2D array of integer counters indexed by a LSH function. These counters can tell whether we have seen data that is similar to a new sequence. The key idea is that we only store sequences if we haven't seen something similar before. This gives us a diverse sample. 

There are a couple of hyperparameters that you may want to change: 

- tau: This is a threshold for whether we should keep a new sequence or not. Increasing tau means that we will keep more sequences. Typical values for tau are between 0.1 and 100.0, depending on the size of the sample you wish to retain.
- range: This is the width of the RACE array. If there are many categories or organisms that you want to sample from, increasing range might help you get more diverse results. Increasing the range is essentially free, but keeping it below 5000 will lead to faster processing times. 
- reps: This is the depth of the RACE array. Increasing the reps will directly increase the time needed to process each input sequence, but you will be much less likely to accidentally discard a rare sequence. 
- hashes: This is the number of LSH functions we use for each row of the RACE array. Increasing this will directly increase the processing time but may also let you differentiate between sequences that are closer together in terms of edit distance. For instance, if you wanted to 
- k: This is the size of each k-mer that is fed to the LSH function (MinHash). Increasing k means that we can differentiate between more similar sequences. To differentiate between species in metagenomic studies, we found that k = 6 is a good choice. If you want to differentiate between mutations or organisms within the same species, try a larger value of k. 

A more in-depth explanation is available in our paper. Feel free to contact the authors with any questions. 

## How to run

Once you have the binaries compiled, you can run the algorithm on a test fastq file included with this repository by running 
```
bin/samplerace 0.1 data/SRR1056036.fastq data/output.fastq 
```
RACE will only require about 20 KB of RAM (not the > 10 GB needed by some other diversity sampling methods) and it can process about 10k sequences per second on a 2016 MacBook (with parallelism enabled). The full 
```
samplerace <tau> <input> <output> [--range race_range] [--reps race_reps] [--hashes n_minhashes] [-k kmer_size]
Positional arguments: 
tau: floating point RACE sampling threshold. Roughly determines how many samples you will store
input: path to input data file (.fastq or .fasta extension)
output: path to output sample file (same extension as input)
Optional arguments: 
[--range race_range]: (Optional, default 1000) Hash range for each ACE (B)
[--reps race_reps]: (Optional, default 4) Number of ACE repetitions (R)
[--hashes n_minhashes]: (Optional, default 4) Number of MinHashes for each ACE (n)
[--k kmer_size]: (Optional, default 6) Size of each MinHash k-mer (k)

Example usage:
samplerace 15.0 data/input.fastq data/output.fastq --range 100 --reps 5 --hashes 3 --k 5
samplerace 0.1 data/input.fasta data/output.fasta --range 500 --k 5

```

## Contact 
For questions about installing or using this software, fill out an issue on GitHub and we'll do our best to help. For questions about the RACE algorithm, contact Benjamin Coleman at Rice University. If you use this software, please cite our paper. 



