#ifndef XAMDST_ENGINE_H
#define XAMDST_ENGINE_H

#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

#include "config.h"
#include "input.h"
#include "intervals.h"

typedef struct {
    uint64_t n_reads;
    uint64_t n_mapped;
    uint64_t n_pair_map;
    uint64_t n_pair_all;
    uint64_t n_pair_good;
    uint64_t n_sgltn;
    uint64_t n_read1;
    uint64_t n_read2;
    uint64_t n_dup;
    uint64_t n_rmdup1;
    uint64_t n_rmdup2;
    uint64_t n_diffchr;
    uint64_t n_pstrand;
    uint64_t n_mstrand;
    uint64_t n_qcfail;
    uint64_t n_data;
    uint64_t n_mdata;
    uint64_t n_qual;
    uint64_t n_tgt;
    uint64_t n_flk;
} read_stats_t;

typedef struct {
    uint64_t paired;
    uint64_t clipped;
    uint64_t unmatched;
    uint64_t ambiguous;
} fragment_stats_t;

typedef struct {
    uint64_t depth;
    uint64_t count;
} histogram_pair_t;

/* Sparse histogram backed by HTSlib's khash implementation. */
typedef struct depth_histogram depth_histogram_t;

depth_histogram_t *histogram_create(void);
void histogram_destroy(depth_histogram_t *histogram);
void histogram_add(depth_histogram_t *histogram, uint64_t depth, uint64_t count);
uint64_t histogram_total(const depth_histogram_t *histogram);
uint64_t histogram_sum(const depth_histogram_t *histogram);
/* Return a cached ascending view.  The view remains owned by the histogram
 * and is invalidated automatically if a later histogram_add() occurs. */
const histogram_pair_t *histogram_sorted(const depth_histogram_t *histogram,
                                         size_t *count);
uint64_t histogram_count_at_least(const depth_histogram_t *histogram, uint64_t threshold);
double histogram_median(const depth_histogram_t *histogram);

typedef struct {
    char *name;
    uint64_t length;
    uint64_t data;
    depth_histogram_t *depth;
    uint64_t coverage_data;
    depth_histogram_t *coverage_depth;
} chromosome_summary_t;

typedef struct {
    read_stats_t reads;
    fragment_stats_t fragments;
    /* target_data/flank_data are raw M/=/X bases in 3.1. */
    uint64_t target_data;
    uint64_t target_rmdup_data;
    uint64_t target_coverage_data;
    uint64_t flank_data;
    uint64_t flank_coverage_data;
    uint64_t target_bases;
    uint64_t flank_bases;
    uint64_t target_region_count;
    depth_histogram_t *target_depth;
    depth_histogram_t *target_rmdup_depth;
    depth_histogram_t *target_coverage_depth;
    depth_histogram_t *flank_depth;
    depth_histogram_t *flank_coverage_depth;
    depth_histogram_t *insert_sizes;
    depth_histogram_t *region_means;
    chromosome_summary_t *chromosomes;
    size_t nchromosomes;
} analysis_result_t;

void analysis_result_init(analysis_result_t *result, const interval_set_t *intervals);
void analysis_result_destroy(analysis_result_t *result);

typedef struct {
    void *opaque;
    int (*depth_row)(void *opaque, const char *chromosome, uint64_t position,
                     uint64_t raw, uint64_t rmdup, uint64_t coverage);
    int (*region_row)(void *opaque, const char *chromosome, const depth_region_t *region);
    int (*uncovered_row)(void *opaque, const char *chromosome, uint64_t start, uint64_t end);
} analysis_sink_t;

/* Run all input streams in coordinate order and populate result. */
int analysis_run(const xamdst_config_t *config, const input_set_t *inputs,
                 interval_set_t *intervals, analysis_result_t *result,
                 analysis_sink_t *sink, htsFile *bamout);

#endif
