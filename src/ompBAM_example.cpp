#include "omp.h"

// [[Rcpp::export]]
int idxstats_pbam(
  std::string bam_file, int n_threads_to_use = 1, bool verbose = true
){

  // Ensure number of threads requested < number of system threads available
  unsigned int n_threads_to_really_use = use_threads_omp(n_threads_to_use);

  pbam_in inbam;
  inbam.openFile(bam_file, n_threads_to_really_use);
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrom_count = inbam.obtainChrs(s_chr_names, u32_chr_lens);

  // If obtainChrs returns zero or negative # chromosomes, BAM reading has failed
  if(chrom_count <= 0) return(-1); 
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(chrom_count);

  while(0 == inbam.fillReads()) {
    // OpenMP parallel FOR loop, each thread runs 1 loop simultaneously.
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_really_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      std::vector<uint32_t> read_counter(chrom_count);
      
      // Gets the first read from the thread read storage buffer
      pbam1_t read(inbam.supplyRead(i));
      // Keep looping while reads are valid
      while(read.validate()) {
        // Counts the read if it is mapped to a chromosome
        if(read.refID() >= 0 && read.refID() < chrom_count) {
          read_counter.at(read.refID())++;
        }
        
        // Gets the next read
        read = inbam.supplyRead(i);     
      }
      // Adds the counted reads to the final count
      // #pragma omp critical ensures only 1 thread at a time runs the following
      // block of code.
      #ifdef _OPENMP
      #pragma omp critical
      #endif
      for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
        total_reads.at(j) += read_counter.at(j);
      }
    }
    // At this stage, all threads would have read all their thread-specific reads
    // At the next call to pbam_in::fillReads(), if any reads were not read, it
    // will throw an error and fillReads() will return -1.
    // If we have finished reading the BAM file, fillReads() will return 1.
  }

  inbam.closeFile();

  // Prints out the count summary to console output
  if(verbose) {
    Rcout << bam_file << " summary:\n" << "Name\tLength\tNumber of reads\n";
    for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
      Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
        << total_reads.at(j) << '\n';
    }
  }
  return(0);
}
