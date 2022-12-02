#include "Rcpp.h"
using namespace Rcpp;

// Required to print cout output generated by ompBAM
#define cout Rcpp::Rcout

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kseq.h"
#include "htslib/thread_pool.h"

unsigned int use_threads_hts(int n_threads = 1) {
  return((unsigned int)n_threads);
}

// [[Rcpp::export]]
int idxstats_hts(std::string bam_path, int n_threads_to_use = 1, bool verbose = true){

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
