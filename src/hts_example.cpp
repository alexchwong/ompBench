#include "omp.h"

unsigned int use_threads_hts(int n_threads) {
  return((unsigned int)n_threads);
}

// [[Rcpp::export]]
int idxstats_hts(
  std::string bam_path, int n_threads_to_use = 1, bool verbose = true
){

  // Ensure number of threads requested < number of system threads available
  unsigned int n_threads_to_really_use = use_threads_hts(n_threads_to_use);

  bam1_t *b = bam_init1();
  BGZF *fp = bgzf_open(bam_path.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);

  hts_tpool *p;
  const int queue_size = 64;

  if (n_threads_to_really_use > 1) {
      p = hts_tpool_init(n_threads_to_really_use);
      bgzf_thread_pool(fp, p, queue_size);
  }
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  
  if(header->n_targets <= 0) return(-1); 
  
  int chrom_count = header->n_targets;
  for (int i = 0; i < header->n_targets; ++i)
  {
      s_chr_names.push_back(header->target_name[i]);
      u32_chr_lens.push_back(header->target_len[i]);
  }
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(chrom_count);

  while (bam_read1(fp, b) >= 0)
  {
    if(b->core.tid >= 0 && b->core.tid < chrom_count) {
      total_reads.at(b->core.tid)++;
    }
  }

  bgzf_close(fp);

  // Prints out the count summary to console output
  if(verbose) {
    Rcout << bam_path << " summary:\n" << "Name\tLength\tNumber of reads\n";
    for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
      Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
        << total_reads.at(j) << '\n';
    }
  }
  return(0);
}

#define B_POOL_NUM 100000

// [[Rcpp::export]]
int idxstats_hts_omp(
  std::string bam_path, int n_threads_to_use = 1, bool verbose = true
){

  // Ensure number of threads requested < number of system threads available
  unsigned int n_threads_to_really_use = use_threads_omp(n_threads_to_use);

  BGZF *fp = bgzf_open(bam_path.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);
  
  hts_tpool *p;
  const int queue_size = 64;
  if (n_threads_to_really_use > 1) {
      p = hts_tpool_init(n_threads_to_really_use);
      bgzf_thread_pool(fp, p, queue_size);
  }
  
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  
  if(header->n_targets <= 0) return(-1); 
  
  int chrom_count = header->n_targets;
  for (int i = 0; i < header->n_targets; ++i)
  {
      s_chr_names.push_back(header->target_name[i]);
      u32_chr_lens.push_back(header->target_len[i]);
  }
  
  // Creates a data structure that stores per-chromosome read counts
  std::vector<uint32_t> total_reads(chrom_count);

  // Pre-allocate bam1_t pool
  std::vector<bam1_t *> bpool;
  for(int i = 0; i < B_POOL_NUM; i++) {
    bpool.push_back(bam_init1());
  }

  bool eofyet = false;
  while(!eofyet) {
    
    // Multi-threaded read a number of reads into memory
    int pool_size = 0;
    for(unsigned int i = 0; i < bpool.size(); i++) {
      int ret = bam_read1(fp, bpool.at(i));
      if(ret < 0) {
        break;
      } else {
        pool_size++;
      }
    }
      
    if(pool_size == 0) {
      break;
    }
    // read data divider
    std::vector<int> pool_starts;
    std::vector<int> pool_ends;
    int est_tp_size = 1 + (pool_size / n_threads_to_really_use);
    int curpos = 0;
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      if(curpos + est_tp_size > pool_size) {
        pool_starts.push_back(curpos);
        pool_ends.push_back(pool_size - 1);
        curpos = pool_size;
      } else {
        pool_starts.push_back(curpos);
        pool_ends.push_back(curpos + est_tp_size - 1);
        curpos += est_tp_size;
      }
    }
    for(unsigned int i = pool_starts.size(); i < n_threads_to_really_use; i++) {
      pool_starts.push_back(-1);
      pool_ends.push_back(-1);
    }
    
    // OpenMP parallel FOR loop, each thread runs 1 loop simultaneously.
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_really_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_really_use; i++) {
      if(pool_starts.at(i) >= 0) {
        std::vector<uint32_t> read_counter(chrom_count);
        
        // Gets the first read from the thread read storage buffer
        // Keep looping while reads are valid
        bam1_t *b;
        for(int j = pool_starts.at(i); j < pool_ends.at(i); j++) {
          // Counts the read if it is mapped to a chromosome
          
          b = bpool.at(j);
          if(b->core.tid >= 0 && b->core.tid < chrom_count) {
            total_reads.at(b->core.tid)++;
          }  
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
    }
  }
  
  for(unsigned int i = 0; i < bpool.size(); i++) {
    bam_destroy1(bpool.at(i));
  }
  bam_hdr_destroy(header);
  bgzf_close(fp);

  // Prints out the count summary to console output
  if(verbose) {
    Rcout << bam_path << " summary:\n" << "Name\tLength\tNumber of reads\n";
    for(unsigned int j = 0; j < (unsigned int)chrom_count; j++) {
      Rcout << s_chr_names.at(j) << '\t' << u32_chr_lens.at(j) << '\t'
        << total_reads.at(j) << '\n';
    }
  }
  return(0);
}