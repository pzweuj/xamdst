#ifndef XAMDST_INTERVALS_H
#define XAMDST_INTERVALS_H

#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

#include "config.h"

typedef struct {
    uint64_t start; /* zero-based, half-open */
    uint64_t end;
    uint64_t length;
    int64_t *raw_diff;
    int64_t *rmdup_diff;
    int64_t *cov_diff;
    uint64_t *dirty_bits;
    size_t dirty_words;
    uint64_t raw_sum;
    uint64_t rmdup_sum;
    uint64_t cov_sum;
    uint64_t raw_covered;
    uint64_t rmdup_covered;
    uint64_t covered;
    double mean;
    double median;
    double cov_mean;
    double cov_median;
} depth_region_t;

typedef struct {
    depth_region_t *items;
    size_t count;
    size_t capacity;
} region_vec_t;

typedef struct {
    int tid;
    char *name;
    uint64_t length;
    region_vec_t target;
    region_vec_t flank;
    struct coverage_window *windows;
    size_t nwindows;
} chromosome_intervals_t;

typedef struct coverage_window {
    depth_region_t *region;
    int target;
} coverage_window_t;

typedef struct {
    chromosome_intervals_t *chromosomes;
    size_t nchromosomes;
    uint64_t target_bases;
    uint64_t flank_bases;
} interval_set_t;

int intervals_load(interval_set_t *set, const char *path, const sam_hdr_t *header,
                   int one_based, int flank);
void intervals_destroy(interval_set_t *set);

int region_ensure_buffers(depth_region_t *region, int summary_only);
void regions_release_buffers(region_vec_t *regions);

#endif
