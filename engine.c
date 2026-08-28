#include "engine.h"

#include "util.h"

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include <htslib/khash.h>

KHASH_MAP_INIT_INT64(xamdst_depth, uint64_t)

enum { HISTOGRAM_DENSE_MAX = 65535 };

struct depth_histogram {
    uint64_t *dense;
    size_t dense_capacity;
    khash_t(xamdst_depth) *overflow;
    uint64_t total;
    uint64_t sum;
    histogram_pair_t *sorted_cache;
    size_t sorted_count;
    int sorted_valid;
};

static int histogram_pair_compare(const void *lhs, const void *rhs)
{
    const histogram_pair_t *a = (const histogram_pair_t *)lhs;
    const histogram_pair_t *b = (const histogram_pair_t *)rhs;
    if (a->depth < b->depth) return -1;
    if (a->depth > b->depth) return 1;
    return 0;
}

depth_histogram_t *histogram_create(void)
{
    depth_histogram_t *histogram = xcalloc(1, sizeof(*histogram));
    histogram->overflow = kh_init(xamdst_depth);
    if (histogram->overflow == NULL) {
        xerror("failed to initialize depth histogram");
        exit(EXIT_FAILURE);
    }
    return histogram;
}

void histogram_destroy(depth_histogram_t *histogram)
{
    if (histogram == NULL)
        return;
    kh_destroy(xamdst_depth, histogram->overflow);
    free(histogram->dense);
    free(histogram->sorted_cache);
    free(histogram);
}

static void histogram_grow_dense(depth_histogram_t *histogram, size_t needed)
{
    if (needed <= histogram->dense_capacity)
        return;
    size_t next = histogram->dense_capacity == 0 ? 256 : histogram->dense_capacity;
    while (next < needed && next < HISTOGRAM_DENSE_MAX + (size_t)1) {
        if (next > (HISTOGRAM_DENSE_MAX + (size_t)1) / 2) {
            next = HISTOGRAM_DENSE_MAX + (size_t)1;
            break;
        }
        next *= 2;
    }
    if (next < needed)
        next = needed;
    if (next > HISTOGRAM_DENSE_MAX + (size_t)1)
        next = HISTOGRAM_DENSE_MAX + (size_t)1;
    histogram->dense = xreallocarray(histogram->dense, next, sizeof(*histogram->dense));
    if (next > histogram->dense_capacity)
        memset(histogram->dense + histogram->dense_capacity, 0,
               (next - histogram->dense_capacity) * sizeof(*histogram->dense));
    histogram->dense_capacity = next;
}

void histogram_add(depth_histogram_t *histogram, uint64_t depth, uint64_t count)
{
    if (count == 0)
        return;
    if (histogram->sorted_valid) {
        free(histogram->sorted_cache);
        histogram->sorted_cache = NULL;
        histogram->sorted_count = 0;
        histogram->sorted_valid = 0;
    }
    if (depth <= HISTOGRAM_DENSE_MAX) {
        histogram_grow_dense(histogram, (size_t)depth + 1);
        if (UINT64_MAX - histogram->dense[depth] < count) {
            xerror("depth histogram count overflow");
            exit(EXIT_FAILURE);
        }
        histogram->dense[depth] += count;
    } else {
        int ret = 0;
        khiter_t key = kh_put(xamdst_depth, histogram->overflow, depth, &ret);
        if (ret < 0 || key == kh_end(histogram->overflow)) {
            xerror("failed to grow depth histogram");
            exit(EXIT_FAILURE);
        }
        if (ret == 1)
            kh_val(histogram->overflow, key) = 0;
        if (UINT64_MAX - kh_val(histogram->overflow, key) < count) {
            xerror("depth histogram count overflow");
            exit(EXIT_FAILURE);
        }
        kh_val(histogram->overflow, key) += count;
    }
    if (UINT64_MAX - histogram->total < count) {
        xerror("depth histogram total overflow");
        exit(EXIT_FAILURE);
    }
    if (depth != 0 && count > UINT64_MAX / depth) {
        xerror("histogram sum overflow");
        exit(EXIT_FAILURE);
    }
    uint64_t term = depth * count;
    if (UINT64_MAX - histogram->sum < term) {
        xerror("histogram sum overflow");
        exit(EXIT_FAILURE);
    }
    histogram->total += count;
    histogram->sum += term;
}

uint64_t histogram_total(const depth_histogram_t *histogram)
{
    return histogram != NULL ? histogram->total : 0;
}

uint64_t histogram_sum(const depth_histogram_t *histogram)
{
    return histogram != NULL ? histogram->sum : 0;
}

const histogram_pair_t *histogram_sorted(const depth_histogram_t *histogram,
                                         size_t *count_out)
{
    depth_histogram_t *mutable_histogram = (depth_histogram_t *)histogram;
    if (!mutable_histogram->sorted_valid) {
        size_t count = (size_t)kh_size(histogram->overflow);
        for (size_t i = 0; i < histogram->dense_capacity; ++i)
            if (histogram->dense[i] != 0) {
                if (count == SIZE_MAX) {
                    xerror("depth histogram has too many keys");
                    exit(EXIT_FAILURE);
                }
                ++count;
            }
        mutable_histogram->sorted_cache = xcalloc(count == 0 ? 1 : count,
                                                  sizeof(*mutable_histogram->sorted_cache));
        mutable_histogram->sorted_count = count;
        size_t index = 0;
        for (size_t i = 0; i < histogram->dense_capacity; ++i) {
            if (histogram->dense[i] == 0)
                continue;
            mutable_histogram->sorted_cache[index++] =
                (histogram_pair_t){(uint64_t)i, histogram->dense[i]};
        }
        for (khiter_t key = kh_begin(histogram->overflow);
             key != kh_end(histogram->overflow); ++key) {
            if (!kh_exist(histogram->overflow, key))
                continue;
            mutable_histogram->sorted_cache[index++] = (histogram_pair_t){
                (uint64_t)kh_key(histogram->overflow, key), kh_val(histogram->overflow, key)};
        }
        qsort(mutable_histogram->sorted_cache, index,
              sizeof(*mutable_histogram->sorted_cache), histogram_pair_compare);
        mutable_histogram->sorted_count = index;
        mutable_histogram->sorted_valid = 1;
    }
    if (count_out != NULL)
        *count_out = histogram->sorted_count;
    return histogram->sorted_cache;
}

uint64_t histogram_count_at_least(const depth_histogram_t *histogram, uint64_t threshold)
{
    uint64_t result = 0;
    size_t first = threshold > HISTOGRAM_DENSE_MAX ? histogram->dense_capacity : (size_t)threshold;
    if (first < histogram->dense_capacity) {
        for (size_t i = first; i < histogram->dense_capacity; ++i) {
            if (histogram->dense[i] == 0)
                continue;
            if (UINT64_MAX - result < histogram->dense[i]) {
                xerror("histogram count overflow");
                exit(EXIT_FAILURE);
            }
            result += histogram->dense[i];
        }
    }
    for (khiter_t key = kh_begin(histogram->overflow);
         key != kh_end(histogram->overflow); ++key) {
        if (!kh_exist(histogram->overflow, key) ||
            (uint64_t)kh_key(histogram->overflow, key) < threshold)
            continue;
        if (UINT64_MAX - result < kh_val(histogram->overflow, key)) {
            xerror("histogram count overflow");
            exit(EXIT_FAILURE);
        }
        result += kh_val(histogram->overflow, key);
    }
    return result;
}

double histogram_median(const depth_histogram_t *histogram)
{
    uint64_t total = histogram_total(histogram);
    if (total == 0)
        return 0.0;
    size_t count = 0;
    const histogram_pair_t *pairs = histogram_sorted(histogram, &count);
    /* Avoid total + 1 wrapping for a saturated 64-bit histogram. */
    uint64_t lower = total / 2 + (total & 1);
    uint64_t upper = total / 2 + 1;
    uint64_t cumulative = 0;
    uint64_t lower_value = 0;
    uint64_t upper_value = 0;
    int found_lower = 0;
    for (size_t i = 0; i < count; ++i) {
        if (UINT64_MAX - cumulative < pairs[i].count) {
            xerror("histogram cumulative count overflow");
            exit(EXIT_FAILURE);
        }
        cumulative += pairs[i].count;
        if (!found_lower && cumulative >= lower) {
            lower_value = pairs[i].depth;
            found_lower = 1;
        }
        if (cumulative >= upper) {
            upper_value = pairs[i].depth;
            break;
        }
    }
    if (total & 1)
        return (double)lower_value;
    return ((double)lower_value + (double)upper_value) / 2.0;
}

void analysis_result_init(analysis_result_t *result, const interval_set_t *intervals)
{
    memset(result, 0, sizeof(*result));
    result->target_bases = intervals->target_bases;
    result->flank_bases = intervals->flank_bases;
    result->target_region_count = 0;
    result->target_depth = histogram_create();
    result->target_rmdup_depth = histogram_create();
    result->target_coverage_depth = histogram_create();
    result->flank_depth = histogram_create();
    result->flank_coverage_depth = histogram_create();
    result->insert_sizes = histogram_create();
    result->region_means = histogram_create();
    result->nchromosomes = intervals->nchromosomes;
    result->chromosomes = xcalloc(result->nchromosomes == 0 ? 1 : result->nchromosomes,
                                  sizeof(*result->chromosomes));
    for (size_t i = 0; i < result->nchromosomes; ++i) {
        const chromosome_intervals_t *source = &intervals->chromosomes[i];
        chromosome_summary_t *summary = &result->chromosomes[i];
        summary->name = xstrdup(source->name);
        summary->depth = histogram_create();
        summary->coverage_depth = histogram_create();
        for (size_t j = 0; j < source->target.count; ++j) {
            uint64_t region_length = source->target.items[j].length;
            if (UINT64_MAX - summary->length < region_length) {
                xerror("chromosome target length overflow");
                exit(EXIT_FAILURE);
            }
            summary->length += region_length;
        }
        if (source->target.count > UINT64_MAX - result->target_region_count) {
            xerror("target region count overflow");
            exit(EXIT_FAILURE);
        }
        result->target_region_count += source->target.count;
    }
}

void analysis_result_destroy(analysis_result_t *result)
{
    if (result == NULL)
        return;
    histogram_destroy(result->target_depth);
    histogram_destroy(result->target_rmdup_depth);
    histogram_destroy(result->target_coverage_depth);
    histogram_destroy(result->flank_depth);
    histogram_destroy(result->flank_coverage_depth);
    histogram_destroy(result->insert_sizes);
    histogram_destroy(result->region_means);
    for (size_t i = 0; i < result->nchromosomes; ++i) {
        free(result->chromosomes[i].name);
        histogram_destroy(result->chromosomes[i].depth);
        histogram_destroy(result->chromosomes[i].coverage_depth);
    }
    free(result->chromosomes);
    memset(result, 0, sizeof(*result));
}

typedef struct {
    htsFile *file;
    bam1_t *record;
    int has_record;
    int saw_record;
    int last_unmapped;
    int last_tid;
    hts_pos_t last_pos;
    size_t source_index;
} stream_cursor_t;

static int read_next(stream_cursor_t *cursor, sam_hdr_t *header)
{
    int ret = sam_read1(cursor->file, header, cursor->record);
    if (ret == -1) {
        cursor->has_record = 0;
        return 0;
    }
    if (ret < -1) {
        xerror("input %zu is truncated or malformed (HTSlib code %d)",
               cursor->source_index + 1, ret);
        return -1;
    }
    int tid = cursor->record->core.tid;
    hts_pos_t pos = cursor->record->core.pos;
    int unmapped = tid < 0 || (cursor->record->core.flag & BAM_FUNMAP) != 0;
    /* Coordinate-sorted streams may contain a terminal unmapped block.  Keep
     * the rule independent of whether an unmapped SAM record happens to carry
     * a reference id, so a mapped record can never follow that block. */
    int out_of_order = cursor->saw_record &&
        ((!unmapped && cursor->last_unmapped) ||
         (!unmapped && !cursor->last_unmapped &&
          (tid < cursor->last_tid || (tid == cursor->last_tid && pos < cursor->last_pos))));
    if (out_of_order) {
        xerror("input %zu is not coordinate sorted near tid %d position %lld",
               cursor->source_index + 1, tid, (long long)pos);
        return -1;
    }
    cursor->saw_record = 1;
    cursor->last_unmapped = unmapped;
    cursor->last_tid = tid;
    cursor->last_pos = pos;
    cursor->has_record = 1;
    return 0;
}

static int cursor_before(const stream_cursor_t *a, const stream_cursor_t *b)
{
    const bam1_core_t *ca = &a->record->core;
    const bam1_core_t *cb = &b->record->core;
    int a_unmapped = ca->tid < 0 || (ca->flag & BAM_FUNMAP) != 0;
    int b_unmapped = cb->tid < 0 || (cb->flag & BAM_FUNMAP) != 0;
    if (a_unmapped != b_unmapped)
        return !a_unmapped;
    if (ca->tid != cb->tid)
        return ca->tid < cb->tid;
    if (ca->pos != cb->pos)
        return ca->pos < cb->pos;
    return a->source_index < b->source_index;
}

static void heap_swap(size_t *heap, size_t a, size_t b)
{
    size_t temp = heap[a];
    heap[a] = heap[b];
    heap[b] = temp;
}

static void heap_push(size_t *heap, size_t *size, size_t value,
                      const stream_cursor_t *cursors)
{
    size_t position = (*size)++;
    heap[position] = value;
    while (position > 0) {
        size_t parent = (position - 1) / 2;
        if (cursor_before(&cursors[heap[parent]], &cursors[heap[position]]))
            break;
        heap_swap(heap, parent, position);
        position = parent;
    }
}

static size_t heap_pop(size_t *heap, size_t *size, const stream_cursor_t *cursors)
{
    size_t result = heap[0];
    heap[0] = heap[--*size];
    size_t position = 0;
    for (;;) {
        size_t left = position * 2 + 1;
        if (left >= *size)
            break;
        size_t right = left + 1;
        size_t smallest = left;
        if (right < *size && cursor_before(&cursors[heap[right]], &cursors[heap[left]]))
            smallest = right;
        if (cursor_before(&cursors[heap[position]], &cursors[heap[smallest]]))
            break;
        heap_swap(heap, position, smallest);
        position = smallest;
    }
    return result;
}

static void add_diff(int64_t *array, size_t index, int64_t amount)
{
    if ((amount > 0 && array[index] > INT64_MAX - amount) ||
        (amount < 0 && array[index] < INT64_MIN - amount)) {
        xerror("depth difference overflow");
        exit(EXIT_FAILURE);
    }
    array[index] += amount;
}

static int add_level(int64_t *level, int64_t delta)
{
    if ((delta > 0 && *level > INT64_MAX - delta) ||
        (delta < 0 && *level < INT64_MIN - delta))
        return -1;
    *level += delta;
    return 0;
}

static int increment_counter(uint64_t *counter, const char *name)
{
    if (*counter == UINT64_MAX) {
        xerror("%s counter overflow", name);
        return -1;
    }
    ++*counter;
    return 0;
}

static int add_counter(uint64_t *counter, uint64_t amount, const char *name)
{
    if (amount > UINT64_MAX - *counter) {
        xerror("%s counter overflow", name);
        return -1;
    }
    *counter += amount;
    return 0;
}

/* Add a signed coverage event to one region.  The caller has already found
 * the intersecting region, so this helper contains no search and can be used
 * by the unified target/flank window walk below. */
static int update_region_delta(depth_region_t *region, uint64_t start, uint64_t end,
                               int raw_delta, int rmdup_delta, int cov_delta,
                               int summary_only)
{
    if (region == NULL || start >= end)
        return 0;
    uint64_t overlap_start = start > region->start ? start : region->start;
    uint64_t overlap_end = end < region->end ? end : region->end;
    if (overlap_start >= overlap_end)
        return 0;
    size_t left = (size_t)(overlap_start - region->start);
    size_t right = (size_t)(overlap_end - region->start);
    if (region_ensure_buffers(region, summary_only) != 0)
        return -1;
    if (raw_delta != 0) {
        add_diff(region->raw_diff, left, raw_delta);
        add_diff(region->raw_diff, right, -raw_delta);
    }
    if (rmdup_delta != 0) {
        add_diff(region->rmdup_diff, left, rmdup_delta);
        add_diff(region->rmdup_diff, right, -rmdup_delta);
    }
    if (cov_delta != 0) {
        add_diff(region->cov_diff, left, cov_delta);
        add_diff(region->cov_diff, right, -cov_delta);
    }
    if (summary_only && region->dirty_bits != NULL) {
        /* The bitmap is a boundary index, rather than a covered-span mask.
         * Recording just the two endpoints lets summary-only finalization
         * jump directly between constant-depth spans.  A bit may remain set
         * after opposing events cancel; that is harmless and only creates a
         * cheap zero-length span check. */
        region->dirty_bits[left / 64] |= UINT64_C(1) << (left % 64);
        region->dirty_bits[right / 64] |= UINT64_C(1) << (right % 64);
    }
    return 1;
}

/* Target and flank regions are disjoint and are stored in one sorted window
 * view.  Locate the first window once per CIGAR reference segment, then walk
 * both target and flank regions through the same loop.  This removes the two
 * independent binary searches that used to dominate high-depth panel runs. */
static int update_chromosome_delta(chromosome_intervals_t *chromosome,
                                   uint64_t start, uint64_t end,
                                   int raw_delta, int rmdup_delta, int cov_delta,
                                   int summary_only)
{
    if (chromosome == NULL || start >= end || chromosome->nwindows == 0)
        return 0;
    size_t low = 0;
    size_t high = chromosome->nwindows;
    while (low < high) {
        size_t middle = low + (high - low) / 2;
        if (chromosome->windows[middle].region->end <= start)
            low = middle + 1;
        else
            high = middle;
    }
    int hit = 0;
    for (size_t i = low; i < chromosome->nwindows; ++i) {
        coverage_window_t *window = &chromosome->windows[i];
        if (window->region->start >= end)
            break;
        int result = update_region_delta(window->region, start, end,
                                         raw_delta, rmdup_delta, cov_delta,
                                         summary_only);
        if (result < 0)
            return -1;
        if (result > 0)
            hit |= window->target ? 1 : 2;
    }
    return hit;
}

typedef struct {
    uint64_t start;
    uint64_t end;
    unsigned kind; /* bit 1: M/=/X (raw), bit 2: M/=/X/D (coverage) */
} fragment_segment_t;

typedef struct {
    size_t input_index;
    int tid;
    hts_pos_t pos;
    hts_pos_t mate_pos;
    uint16_t flag;
    int clean;
    char *qname;
    char *rg;
    fragment_segment_t *segments;
    size_t nsegments;
} fragment_pending_t;

typedef struct {
    const xamdst_config_t *config;
    interval_set_t *intervals;
    analysis_result_t *result;
    analysis_sink_t *sink;
    fragment_pending_t *pending;
    size_t npending;
    size_t pending_capacity;
} run_context_t;

static int checked_sum_product(uint64_t *sum, uint64_t value, uint64_t count,
                               const char *what)
{
    if (value != 0 && count > UINT64_MAX / value) {
        xerror("%s overflow", what);
        return -1;
    }
    uint64_t product = value * count;
    if (product > UINT64_MAX - *sum) {
        xerror("%s overflow", what);
        return -1;
    }
    *sum += product;
    return 0;
}

/* Return the first recorded difference boundary at or after ``offset``.
 * ``region->dirty_bits`` has one bit for every possible diff-array slot,
 * including the sentinel at region->length.  The finalization loop uses this
 * to avoid inspecting every base of a long, sparse region. */
static uint64_t next_summary_boundary(const depth_region_t *region, uint64_t offset)
{
    if (offset >= region->length)
        return region->length;
    size_t word = (size_t)(offset / 64);
    unsigned bit = (unsigned)(offset % 64);
    uint64_t value = region->dirty_bits[word] & (UINT64_MAX << bit);
    for (;;) {
        if (value != 0) {
#if defined(__GNUC__) || defined(__clang__)
            unsigned delta = (unsigned)__builtin_ctzll(value);
#else
            unsigned delta = 0;
            while (((value >> delta) & UINT64_C(1)) == 0)
                ++delta;
#endif
            uint64_t result = (uint64_t)word * 64U + delta;
            return result <= region->length ? result : region->length;
        }
        ++word;
        if (word >= region->dirty_words)
            return region->length;
        value = region->dirty_bits[word];
    }
}

/* Summary-only output never needs individual depth rows.  Walk only the
 * constant-depth spans delimited by difference events, which is substantially
 * cheaper for sparse high-depth panels while preserving all histogram and
 * region statistics. */
static int finalize_region_summary(run_context_t *context,
                                   chromosome_intervals_t *chromosome,
                                   depth_region_t *region, int target)
{
    depth_histogram_t *local_raw = target ? histogram_create() : NULL;
    depth_histogram_t *local_cov = target ? histogram_create() : NULL;
    int64_t raw_level = 0, rmdup_level = 0, cov_level = 0;
    uint64_t uncover_start = 0;
    int in_uncover = 0;
    uint64_t offset = 0;
    while (offset < region->length) {
        if (add_level(&raw_level, region->raw_diff[offset]) ||
            add_level(&rmdup_level, region->rmdup_diff[offset]) ||
            add_level(&cov_level, region->cov_diff[offset]) ||
            raw_level < 0 || rmdup_level < 0 || cov_level < 0) {
            xerror("internal negative/overflow depth while finalizing %s:%llu-%llu",
                   chromosome->name, (unsigned long long)region->start,
                   (unsigned long long)region->end);
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        uint64_t next = next_summary_boundary(region, offset + 1);
        if (next < offset + 1)
            next = offset + 1;
        uint64_t span = next - offset;
        uint64_t raw = (uint64_t)raw_level;
        uint64_t rmdup = (uint64_t)rmdup_level;
        uint64_t coverage = (uint64_t)cov_level;
        if (checked_sum_product(&region->raw_sum, raw, span, "raw depth sum") ||
            checked_sum_product(&region->rmdup_sum, rmdup, span, "rmdup depth sum") ||
            checked_sum_product(&region->cov_sum, coverage, span, "coverage depth sum")) {
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        if (target) {
            histogram_add(context->result->target_depth, raw, span);
            histogram_add(context->result->target_rmdup_depth, rmdup, span);
            histogram_add(context->result->target_coverage_depth, coverage, span);
            histogram_add(context->result->chromosomes[chromosome->tid].depth, raw, span);
            histogram_add(context->result->chromosomes[chromosome->tid].coverage_depth,
                          coverage, span);
            histogram_add(local_raw, raw, span);
            histogram_add(local_cov, coverage, span);
        } else {
            histogram_add(context->result->flank_depth, raw, span);
            histogram_add(context->result->flank_coverage_depth, coverage, span);
        }
        if (raw > 0) {
            if (span > UINT64_MAX - region->raw_covered) {
                xerror("raw covered count overflow");
                histogram_destroy(local_raw);
                histogram_destroy(local_cov);
                return -1;
            }
            region->raw_covered += span;
        }
        if (rmdup > 0) {
            if (span > UINT64_MAX - region->rmdup_covered) {
                xerror("rmdup covered count overflow");
                histogram_destroy(local_raw);
                histogram_destroy(local_cov);
                return -1;
            }
            region->rmdup_covered += span;
        }
        if (coverage > 0) {
            if (span > UINT64_MAX - region->covered) {
                xerror("coverage count overflow");
                histogram_destroy(local_raw);
                histogram_destroy(local_cov);
                return -1;
            }
            region->covered += span;
        }
        int uncover = target && raw < (uint64_t)context->config->uncover;
        if (uncover && !in_uncover) {
            uncover_start = region->start + offset;
            in_uncover = 1;
        } else if (!uncover && in_uncover) {
            if (context->sink->uncovered_row(context->sink->opaque, chromosome->name,
                                             uncover_start, region->start + offset)) {
                histogram_destroy(local_raw);
                histogram_destroy(local_cov);
                return -1;
            }
            in_uncover = 0;
        }
        offset = next;
    }
    if (add_level(&raw_level, region->raw_diff[region->length]) ||
        add_level(&rmdup_level, region->rmdup_diff[region->length]) ||
        add_level(&cov_level, region->cov_diff[region->length]) ||
        raw_level != 0 || rmdup_level != 0 || cov_level != 0) {
        xerror("unbalanced depth differences while finalizing %s:%llu-%llu",
               chromosome->name, (unsigned long long)region->start,
               (unsigned long long)region->end);
        histogram_destroy(local_raw);
        histogram_destroy(local_cov);
        return -1;
    }
    if (in_uncover && context->sink->uncovered_row(context->sink->opaque, chromosome->name,
                                                   uncover_start, region->end)) {
        histogram_destroy(local_raw);
        histogram_destroy(local_cov);
        return -1;
    }
    if (region->length > 0) {
        region->mean = (double)region->raw_sum / (double)region->length;
        region->cov_mean = (double)region->cov_sum / (double)region->length;
        region->median = target ? histogram_median(local_raw) : 0.0;
        region->cov_median = target ? histogram_median(local_cov) : 0.0;
    }
    if (target) {
        chromosome_summary_t *summary = &context->result->chromosomes[chromosome->tid];
        if (UINT64_MAX - context->result->target_data < region->raw_sum ||
            UINT64_MAX - context->result->target_rmdup_data < region->rmdup_sum ||
            UINT64_MAX - context->result->target_coverage_data < region->cov_sum ||
            UINT64_MAX - summary->data < region->raw_sum ||
            UINT64_MAX - summary->coverage_data < region->cov_sum) {
            xerror("target depth total overflow");
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        context->result->target_data += region->raw_sum;
        context->result->target_rmdup_data += region->rmdup_sum;
        context->result->target_coverage_data += region->cov_sum;
        summary->data += region->raw_sum;
        summary->coverage_data += region->cov_sum;
        histogram_add(context->result->region_means, (uint64_t)region->mean, 1);
        if (context->sink->region_row(context->sink->opaque, chromosome->name, region)) {
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
    } else {
        if (UINT64_MAX - context->result->flank_data < region->raw_sum ||
            UINT64_MAX - context->result->flank_coverage_data < region->cov_sum) {
            xerror("flank depth total overflow");
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        context->result->flank_data += region->raw_sum;
        context->result->flank_coverage_data += region->cov_sum;
    }
    histogram_destroy(local_raw);
    histogram_destroy(local_cov);
    return 0;
}

static int finalize_region(run_context_t *context, chromosome_intervals_t *chromosome,
                           depth_region_t *region, int target)
{
    uint64_t length = region->length;
    region->raw_sum = region->rmdup_sum = region->cov_sum = 0;
    region->raw_covered = region->rmdup_covered = region->covered = 0;
    region->mean = region->median = region->cov_mean = region->cov_median = 0.0;

    /* A lazily allocated region has never intersected a CIGAR segment, so
     * every base is exactly zero.  Account for the whole span in one
     * histogram update; target depth rows still have to be emitted because
     * they are part of the public per-position report. */
    if (region->raw_diff == NULL) {
        if (target) {
            histogram_add(context->result->target_depth, 0, length);
            histogram_add(context->result->target_rmdup_depth, 0, length);
            histogram_add(context->result->target_coverage_depth, 0, length);
            histogram_add(context->result->chromosomes[chromosome->tid].depth, 0, length);
            histogram_add(context->result->chromosomes[chromosome->tid].coverage_depth, 0, length);
            histogram_add(context->result->region_means, 0, 1);
            if (context->sink->uncovered_row(context->sink->opaque, chromosome->name,
                                             region->start, region->end))
                return -1;
            if (!context->config->summary_only)
                for (uint64_t offset = 0; offset < length; ++offset)
                    if (context->sink->depth_row(context->sink->opaque, chromosome->name,
                                                 region->start + offset + 1, 0, 0, 0))
                        return -1;
            if (context->sink->region_row(context->sink->opaque, chromosome->name, region))
                return -1;
        } else {
            histogram_add(context->result->flank_depth, 0, length);
            histogram_add(context->result->flank_coverage_depth, 0, length);
        }
        return 0;
    }

    if (context->config->summary_only)
        return finalize_region_summary(context, chromosome, region, target);

    depth_histogram_t *local_raw = target ? histogram_create() : NULL;
    depth_histogram_t *local_cov = target ? histogram_create() : NULL;
    int64_t raw_level = 0;
    int64_t rmdup_level = 0;
    int64_t cov_level = 0;
    uint64_t uncover_start = 0;
    int in_uncover = 0;

    for (uint64_t offset = 0; offset < length; ++offset) {
        if (region->raw_diff != NULL) {
            if (add_level(&raw_level, region->raw_diff[offset]) ||
                add_level(&rmdup_level, region->rmdup_diff[offset]) ||
                add_level(&cov_level, region->cov_diff[offset])) {
                xerror("depth level overflow while finalizing %s:%llu-%llu",
                       chromosome->name, (unsigned long long)region->start,
                       (unsigned long long)region->end);
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
            }
        }
        if (raw_level < 0 || rmdup_level < 0 || cov_level < 0) {
            xerror("internal negative depth while finalizing %s:%llu-%llu",
                   chromosome->name, (unsigned long long)region->start,
                   (unsigned long long)region->end);
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        uint64_t raw = (uint64_t)raw_level;
        uint64_t rmdup = (uint64_t)rmdup_level;
        uint64_t coverage = (uint64_t)cov_level;
        if (UINT64_MAX - region->raw_sum < raw ||
            UINT64_MAX - region->rmdup_sum < rmdup ||
            UINT64_MAX - region->cov_sum < coverage) {
            xerror("depth sum overflow while finalizing %s:%llu-%llu",
                   chromosome->name, (unsigned long long)region->start,
                   (unsigned long long)region->end);
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        region->raw_sum += raw;
        region->rmdup_sum += rmdup;
        region->cov_sum += coverage;
        if (local_raw != NULL) {
            histogram_add(local_raw, raw, 1);
            histogram_add(local_cov, coverage, 1);
        }
        if (target) {
            histogram_add(context->result->target_depth, raw, 1);
            histogram_add(context->result->target_rmdup_depth, rmdup, 1);
            histogram_add(context->result->target_coverage_depth, coverage, 1);
            histogram_add(context->result->chromosomes[chromosome->tid].depth, raw, 1);
            histogram_add(context->result->chromosomes[chromosome->tid].coverage_depth,
                          coverage, 1);
        } else {
            histogram_add(context->result->flank_depth, raw, 1);
            histogram_add(context->result->flank_coverage_depth, coverage, 1);
        }
        if (raw > 0) ++region->raw_covered;
        if (rmdup > 0) ++region->rmdup_covered;
        if (coverage > 0) ++region->covered;

        /* uncover.bed follows the public raw-depth series in 3.1. */
        int uncover = target && raw < (uint64_t)context->config->uncover;
        if (uncover && !in_uncover) {
            uncover_start = region->start + offset;
            in_uncover = 1;
        } else if (!uncover && in_uncover) {
            if (context->sink->uncovered_row(context->sink->opaque, chromosome->name,
                                             uncover_start, region->start + offset)) {
                histogram_destroy(local_raw);
                histogram_destroy(local_cov);
                return -1;
            }
            in_uncover = 0;
        }
        if (target && !context->config->summary_only &&
            context->sink->depth_row(context->sink->opaque, chromosome->name,
                                     region->start + offset + 1, raw, rmdup, coverage)) {
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
    }
    if (add_level(&raw_level, region->raw_diff[length]) ||
        add_level(&rmdup_level, region->rmdup_diff[length]) ||
        add_level(&cov_level, region->cov_diff[length]) ||
        raw_level != 0 || rmdup_level != 0 || cov_level != 0) {
        xerror("unbalanced depth differences while finalizing %s:%llu-%llu",
               chromosome->name, (unsigned long long)region->start,
               (unsigned long long)region->end);
        histogram_destroy(local_raw);
        histogram_destroy(local_cov);
        return -1;
    }
    if (in_uncover && context->sink->uncovered_row(context->sink->opaque, chromosome->name,
                                                   uncover_start, region->end)) {
        histogram_destroy(local_raw);
        histogram_destroy(local_cov);
        return -1;
    }
    if (length > 0) {
        region->mean = (double)region->raw_sum / (double)length;
        region->cov_mean = (double)region->cov_sum / (double)length;
        region->median = target ? histogram_median(local_raw) : 0.0;
        region->cov_median = target ? histogram_median(local_cov) : 0.0;
    }
    if (target) {
        if (UINT64_MAX - context->result->target_data < region->raw_sum ||
            UINT64_MAX - context->result->target_rmdup_data < region->rmdup_sum ||
            UINT64_MAX - context->result->target_coverage_data < region->cov_sum ||
            UINT64_MAX - context->result->chromosomes[chromosome->tid].data < region->raw_sum ||
            UINT64_MAX - context->result->chromosomes[chromosome->tid].coverage_data < region->cov_sum) {
            xerror("target depth total overflow");
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        context->result->target_data += region->raw_sum;
        context->result->target_rmdup_data += region->rmdup_sum;
        context->result->target_coverage_data += region->cov_sum;
        context->result->chromosomes[chromosome->tid].data += region->raw_sum;
        context->result->chromosomes[chromosome->tid].coverage_data += region->cov_sum;
        histogram_add(context->result->region_means, (uint64_t)region->mean, 1);
        if (context->sink->region_row(context->sink->opaque, chromosome->name, region)) {
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
    } else {
        if (UINT64_MAX - context->result->flank_data < region->raw_sum ||
            UINT64_MAX - context->result->flank_coverage_data < region->cov_sum) {
            xerror("flank depth total overflow");
            histogram_destroy(local_raw);
            histogram_destroy(local_cov);
            return -1;
        }
        context->result->flank_data += region->raw_sum;
        context->result->flank_coverage_data += region->cov_sum;
    }
    histogram_destroy(local_raw);
    histogram_destroy(local_cov);
    return 0;
}

static int finish_chromosome(run_context_t *context, size_t index)
{
    chromosome_intervals_t *chromosome = &context->intervals->chromosomes[index];
    for (size_t i = 0; i < chromosome->target.count; ++i)
        if (finalize_region(context, chromosome, &chromosome->target.items[i], 1))
            return -1;
    for (size_t i = 0; i < chromosome->flank.count; ++i)
        if (finalize_region(context, chromosome, &chromosome->flank.items[i], 0))
            return -1;
    regions_release_buffers(&chromosome->target);
    regions_release_buffers(&chromosome->flank);
    return 0;
}

static void fragment_pending_destroy(fragment_pending_t *pending)
{
    if (pending == NULL)
        return;
    free(pending->qname);
    free(pending->rg);
    free(pending->segments);
    memset(pending, 0, sizeof(*pending));
}

static void fragment_remove(run_context_t *context, size_t index, int unmatched)
{
    if (index >= context->npending)
        return;
    if (unmatched)
        (void)increment_counter(&context->result->fragments.unmatched, "unmatched fragment");
    fragment_pending_destroy(&context->pending[index]);
    if (index + 1 < context->npending)
        context->pending[index] = context->pending[context->npending - 1];
    --context->npending;
}

static void fragment_clear(run_context_t *context)
{
    while (context->npending > 0)
        fragment_remove(context, context->npending - 1, 1);
}

static int fragment_add_segment(fragment_segment_t **segments, size_t *count,
                                size_t *capacity, uint64_t start, uint64_t end,
                                unsigned kind)
{
    if (start >= end)
        return 0;
    if (*count == *capacity) {
        size_t next = *capacity == 0 ? 8 : *capacity * 2;
        if (next < *capacity) {
            xerror("too many CIGAR segments in fragment");
            return -1;
        }
        *segments = xreallocarray(*segments, next, sizeof(**segments));
        *capacity = next;
    }
    (*segments)[(*count)++] = (fragment_segment_t){start, end, kind};
    return 0;
}

static const char *record_rg(const bam1_t *record)
{
    uint8_t *value = bam_aux_get(record, "RG");
    return value != NULL ? bam_aux2Z(value) : NULL;
}

static int validate_qname(const bam1_t *record, const char *chromosome)
{
    if (record->core.l_qname == 0 || record->data == NULL || record->l_data < 0 ||
        (size_t)record->core.l_qname > (size_t)record->l_data ||
        record->data[record->core.l_qname - 1] != '\0') {
        xerror("mapped record on %s has an invalid query name", chromosome);
        return -1;
    }
    return 0;
}

static int fragment_match(const fragment_pending_t *pending, size_t input_index,
                          const bam1_core_t *core, const char *qname, const char *rg)
{
    if (pending->input_index != input_index || pending->tid != core->tid ||
        pending->mate_pos < 0 || core->pos != pending->mate_pos ||
        core->mtid != pending->tid || core->mpos != pending->pos ||
        strcmp(pending->qname, qname) != 0)
        return 0;
    if ((pending->rg == NULL) != (rg == NULL) ||
        (pending->rg != NULL && strcmp(pending->rg, rg) != 0))
        return 0;
    int p_read1 = (pending->flag & BAM_FREAD1) != 0;
    int p_read2 = (pending->flag & BAM_FREAD2) != 0;
    int c_read1 = (core->flag & BAM_FREAD1) != 0;
    int c_read2 = (core->flag & BAM_FREAD2) != 0;
    return (p_read1 && c_read2) || (p_read2 && c_read1);
}

static int fragment_subtract(run_context_t *context, const fragment_pending_t *pending,
                             const fragment_segment_t *current, size_t ncurrent,
                             int clean)
{
    int clipped = 0;
    for (size_t i = 0; i < pending->nsegments; ++i) {
        const fragment_segment_t *left = &pending->segments[i];
        for (size_t j = 0; j < ncurrent; ++j) {
            const fragment_segment_t *right = &current[j];
            uint64_t start = left->start > right->start ? left->start : right->start;
            uint64_t end = left->end < right->end ? left->end : right->end;
            if (start >= end)
                continue;
            int raw = (left->kind & 1U) && (right->kind & 1U);
            int coverage = (left->kind & 2U) && (right->kind & 2U);
            if (!raw && !coverage)
                continue;
            if (update_chromosome_delta(&context->intervals->chromosomes[pending->tid],
                                        start, end, raw ? -1 : 0,
                                        raw && pending->clean && clean ? -1 : 0,
                                        coverage ? -1 : 0,
                                        context->config->summary_only) < 0)
                return -1;
            clipped = 1;
        }
    }
    if (clipped) {
        if (increment_counter(&context->result->fragments.clipped, "clipped fragment overlap"))
            return -1;
    }
    return 0;
}

static int fragment_maybe_store(run_context_t *context, size_t input_index,
                                const bam1_core_t *core, const char *qname,
                                const char *rg, int clean,
                                fragment_segment_t *segments, size_t nsegments)
{
    if (!context->config->fragment_mode || !(core->flag & BAM_FPAIRED) ||
        !(core->flag & BAM_FPROPER_PAIR) || (core->flag & BAM_FMUNMAP) ||
        core->mtid != core->tid || core->tid < 0 || core->pos < 0 || core->mpos < 0 ||
        !(((core->flag & BAM_FREAD1) != 0) ^ ((core->flag & BAM_FREAD2) != 0)))
        return 0;
    if (qname == NULL) {
        xerror("eligible fragment record has no query name");
        return -1;
    }

    /* Try an exact reciprocal match before expiring entries.  This keeps the
     * normal fast path independent of the cache cleanup policy. */
    for (size_t i = 0; i < context->npending; ++i) {
        if (!fragment_match(&context->pending[i], input_index, core, qname, rg))
            continue;
        if (fragment_subtract(context, &context->pending[i], segments, nsegments, clean) < 0)
            return -1;
        fragment_remove(context, i, 0);
        if (increment_counter(&context->result->fragments.paired, "paired fragment"))
            return -1;
        return 0;
    }
    /* Same logical key with a conflicting mate coordinate is deliberately not
     * clipped.  Keep a diagnostic counter while retaining both records so a
     * later exact reciprocal pair can still be recognized.  This check must
     * precede expiry: a conflicting record can arrive after the old expected
     * mate position and should still be visible as ambiguous, not only as an
     * unmatched cache entry. */
    for (size_t i = 0; i < context->npending; ++i) {
        const fragment_pending_t *old = &context->pending[i];
        if (old->input_index == input_index && old->tid == core->tid &&
            strcmp(old->qname, qname) == 0 &&
            ((old->rg == NULL && rg == NULL) ||
             (old->rg != NULL && rg != NULL && strcmp(old->rg, rg) == 0))) {
            if (increment_counter(&context->result->fragments.ambiguous,
                                  "ambiguous fragment"))
                return -1;
            break;
        }
    }
    /* A coordinate-sorted stream cannot deliver a mate after its expected
     * coordinate.  Retire such entries promptly to keep the cache bounded. */
    for (size_t i = 0; i < context->npending;) {
        fragment_pending_t *old = &context->pending[i];
        if (old->input_index == input_index && old->tid == core->tid &&
            old->mate_pos >= 0 && core->pos > old->mate_pos)
            fragment_remove(context, i, 1);
        else
            ++i;
    }
    if (context->npending == context->pending_capacity) {
        size_t next = context->pending_capacity == 0 ? 128 : context->pending_capacity * 2;
        if (next < context->pending_capacity) {
            xerror("fragment pending cache overflow");
            return -1;
        }
        context->pending = xreallocarray(context->pending, next, sizeof(*context->pending));
        context->pending_capacity = next;
    }
    fragment_pending_t *entry = &context->pending[context->npending++];
    memset(entry, 0, sizeof(*entry));
    entry->input_index = input_index;
    entry->tid = core->tid;
    entry->pos = core->pos;
    entry->mate_pos = core->mpos;
    entry->flag = core->flag;
    entry->clean = clean;
    entry->qname = xstrdup(qname);
    entry->rg = rg != NULL ? xstrdup(rg) : NULL;
    entry->segments = segments;
    entry->nsegments = nsegments;
    return 1; /* ownership of segments moved to the pending cache */
}

static int process_alignment(run_context_t *context, sam_hdr_t *header,
                             const bam1_t *record, size_t input_index, htsFile *bamout)
{
    bam1_core_t const *core = &record->core;
    if (core->l_qseq < 0) {
        xerror("input contains a negative query length");
        return -1;
    }
    if (core->tid < 0 && !(core->flag & BAM_FUNMAP)) {
        xerror("mapped record has no reference id");
        return -1;
    }
    if (record->l_data < 0) {
        xerror("input contains a negative BAM data length");
        return -1;
    }
    read_stats_t *stats = &context->result->reads;
    if (increment_counter(&stats->n_reads, "read") ||
        add_counter(&stats->n_data, (uint64_t)core->l_qseq, "raw data"))
        return -1;
    int qc_fail = (core->flag & BAM_FQCFAIL) != 0;
    int duplicate = (core->flag & BAM_FDUP) != 0;
    int mapped = core->tid >= 0 && !(core->flag & BAM_FUNMAP);
    if (qc_fail && increment_counter(&stats->n_qcfail, "QC-fail read"))
        return -1;
    if (core->flag & BAM_FPAIRED) {
        if (increment_counter(&stats->n_pair_all, "paired-read") ||
            ((core->flag & BAM_FPROPER_PAIR) && mapped && !(core->flag & BAM_FMUNMAP) &&
             increment_counter(&stats->n_pair_good, "proper-pair")) ||
            ((core->flag & BAM_FREAD1) && increment_counter(&stats->n_read1, "read1")) ||
            ((core->flag & BAM_FREAD2) && increment_counter(&stats->n_read2, "read2")))
            return -1;
        if (mapped && (core->flag & BAM_FMUNMAP) &&
            increment_counter(&stats->n_sgltn, "singleton"))
            return -1;
        if (mapped && !(core->flag & BAM_FMUNMAP)) {
            if (increment_counter(&stats->n_pair_map, "mapped pair") ||
                (core->mtid != core->tid && increment_counter(&stats->n_diffchr, "different-reference pair")))
                return -1;
        }
    }
    if (!mapped)
        return 0;
    if ((size_t)core->tid >= context->intervals->nchromosomes) {
        xerror("input contains reference id %d outside its header", core->tid);
        return -1;
    }
    if (core->pos < 0) {
        xerror("mapped record on %s has a negative position",
               context->intervals->chromosomes[core->tid].name);
        return -1;
    }
    chromosome_intervals_t *chromosome = &context->intervals->chromosomes[core->tid];
    if ((uint64_t)core->pos > chromosome->length) {
        xerror("mapped record on %s starts beyond chromosome length",
               chromosome->name);
        return -1;
    }
    if (core->n_cigar == 0 || record->data == NULL) {
        xerror("mapped record on %s has no CIGAR", chromosome->name);
        return -1;
    }
    size_t cigar_bytes;
    if (size_mul((size_t)core->n_cigar, sizeof(uint32_t), &cigar_bytes) != 0 ||
        (size_t)core->l_qname > (size_t)record->l_data ||
        cigar_bytes > (size_t)record->l_data - (size_t)core->l_qname) {
        xerror("mapped record on %s has a truncated CIGAR", chromosome->name);
        return -1;
    }
    if (validate_qname(record, chromosome->name) != 0)
        return -1;
    if (increment_counter(&stats->n_mapped, "mapped read") ||
        add_counter(&stats->n_mdata, (uint64_t)core->l_qseq, "mapped data") ||
        ((core->flag & BAM_FREVERSE) && increment_counter(&stats->n_mstrand, "reverse-strand read")) ||
        (!(core->flag & BAM_FREVERSE) && increment_counter(&stats->n_pstrand, "forward-strand read")) ||
        (duplicate && increment_counter(&stats->n_dup, "duplicate read")))
        return -1;
    int high_mapq = core->qual >= context->config->mapq;
    if (high_mapq && increment_counter(&stats->n_qual, "mapQ-qualified read"))
        return -1;
    int clean = !qc_fail && !duplicate && high_mapq;
    if (clean) {
        if ((core->flag & BAM_FREAD1 && increment_counter(&stats->n_rmdup1, "rmdup read1")) ||
            (core->flag & BAM_FREAD2 && increment_counter(&stats->n_rmdup2, "rmdup read2")))
            return -1;
    }
    if (core->isize > 0 && core->isize < context->config->isize)
        histogram_add(context->result->insert_sizes, (uint64_t)core->isize, 1);

    uint32_t *cigar = bam_get_cigar(record);
    fragment_segment_t *segments = NULL;
    size_t nsegments = 0;
    size_t segment_capacity = 0;
    uint64_t reference_position = (uint64_t)core->pos;
    uint64_t query_position = 0;
    int target_hit = 0;
    int flank_hit = 0;
    for (uint32_t i = 0; i < core->n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        uint32_t length = bam_cigar_oplen(cigar[i]);
        if (length == 0) {
            xerror("mapped record on %s contains a zero-length CIGAR operation",
                   chromosome->name);
            free(segments);
            return -1;
        }
        int consumes_query = (op == BAM_CMATCH || op == BAM_CINS ||
                              op == BAM_CSOFT_CLIP || op == BAM_CEQUAL ||
                              op == BAM_CDIFF);
        if (consumes_query) {
            if ((uint64_t)length > UINT64_MAX - query_position ||
                query_position + (uint64_t)length > (uint64_t)core->l_qseq) {
                xerror("CIGAR query length exceeds sequence length on %s",
                       chromosome->name);
                free(segments);
                return -1;
            }
            query_position += (uint64_t)length;
        }
        int consumes_reference = (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ||
                                  op == BAM_CDEL || op == BAM_CREF_SKIP);
        if (consumes_reference) {
            if ((uint64_t)length > UINT64_MAX - reference_position) {
                xerror("CIGAR reference position overflow");
                free(segments);
                return -1;
            }
            uint64_t end = reference_position + (uint64_t)length;
            if (end > chromosome->length) {
                xerror("CIGAR of record on %s extends beyond chromosome length",
                       chromosome->name);
                free(segments);
                return -1;
            }
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                if (context->config->fragment_mode &&
                    fragment_add_segment(&segments, &nsegments, &segment_capacity,
                                         reference_position, end, 3U) != 0) {
                    free(segments);
                    return -1;
                }
                int hit_result = update_chromosome_delta(chromosome, reference_position, end,
                                                         1, clean ? 1 : 0, 1,
                                                         context->config->summary_only);
                if (hit_result < 0) {
                    free(segments);
                    return -1;
                }
                /* The unified walk returns a hit for either target or flank;
                 * retain the individual flags for read-level diagnostics. */
                target_hit |= (hit_result & 1) != 0;
                flank_hit |= (hit_result & 2) != 0;
            } else if (op == BAM_CDEL) {
                if (context->config->fragment_mode &&
                    fragment_add_segment(&segments, &nsegments, &segment_capacity,
                                         reference_position, end, 2U) != 0) {
                    free(segments);
                    return -1;
                }
                int hit_result = update_chromosome_delta(chromosome, reference_position, end,
                                                         0, 0, 1,
                                                         context->config->summary_only);
                if (hit_result < 0) {
                    free(segments);
                    return -1;
                }
                target_hit |= (hit_result & 1) != 0;
                flank_hit |= (hit_result & 2) != 0;
            }
            reference_position = end;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CPAD) {
            /* These operations do not advance the reference position. */
        } else {
            xerror("unsupported CIGAR operation %d", op);
            free(segments);
            return -1;
        }
    }
    if (query_position != (uint64_t)core->l_qseq) {
        xerror("CIGAR query length does not match sequence length on %s",
               chromosome->name);
        free(segments);
        return -1;
    }
    if (target_hit && increment_counter(&stats->n_tgt, "target read")) {
        free(segments);
        return -1;
    }
    if (flank_hit && increment_counter(&stats->n_flk, "flank read")) {
        free(segments);
        return -1;
    }
    if (target_hit && bamout != NULL && sam_write1(bamout, header, record) < 0) {
        xerror("failed to write target BAM record");
        free(segments);
        return -1;
    }
    const char *qname = NULL;
    const char *rg = NULL;
    if (context->config->fragment_mode) {
        qname = bam_get_qname(record);
        rg = record_rg(record);
    }
    int fragment_status = fragment_maybe_store(context, input_index, core, qname, rg,
                                               clean, segments, nsegments);
    if (fragment_status < 0) {
        free(segments);
        return -1;
    }
    if (fragment_status == 0)
        free(segments);
    return 0;
}

typedef struct {
    depth_region_t *region;
    uint64_t start;
    uint64_t end;
    unsigned kind; /* bit 1: M/=/X (raw), bit 2: M/=/X/D (coverage) */
} alignment_delta_t;

typedef struct {
    bam1_t *record;
    size_t input_index;
    int status;
    int mapped;
    int clean;
    int target_hit;
    int flank_hit;
    alignment_delta_t *deltas;
    size_t ndeltas;
    size_t delta_capacity;
    fragment_segment_t *segments;
    size_t nsegments;
} alignment_event_t;

static int alignment_add_delta(alignment_event_t *event, depth_region_t *region,
                               uint64_t start, uint64_t end, unsigned kind)
{
    if (region == NULL || start >= end)
        return 0;
    if (event->ndeltas == event->delta_capacity) {
        size_t next = event->delta_capacity == 0 ? 8 : event->delta_capacity * 2;
        if (next < event->delta_capacity) {
            xerror("too many coverage events in one alignment");
            return -1;
        }
        event->deltas = xreallocarray(event->deltas, next, sizeof(*event->deltas));
        event->delta_capacity = next;
    }
    event->deltas[event->ndeltas++] = (alignment_delta_t){region, start, end, kind};
    return 0;
}

/* Expand one CIGAR reference segment into direct region events.  The worker
 * performs the binary search and window walk once; the reducer can then apply
 * each event without repeating interval lookup or touching unrelated
 * windows. */
static int prepare_alignment_deltas(const chromosome_intervals_t *chromosome,
                                   uint64_t start, uint64_t end, unsigned kind,
                                   alignment_event_t *event)
{
    if (chromosome == NULL || start >= end || chromosome->nwindows == 0)
        return 0;
    size_t low = 0, high = chromosome->nwindows;
    while (low < high) {
        size_t middle = low + (high - low) / 2;
        if (chromosome->windows[middle].region->end <= start)
            low = middle + 1;
        else
            high = middle;
    }
    for (size_t i = low; i < chromosome->nwindows; ++i) {
        const coverage_window_t *window = &chromosome->windows[i];
        if (window->region->start >= end)
            break;
        uint64_t overlap_start = start > window->region->start ?
                                 start : window->region->start;
        uint64_t overlap_end = end < window->region->end ? end : window->region->end;
        if (overlap_start >= overlap_end)
            continue;
        if (alignment_add_delta(event, window->region, overlap_start, overlap_end,
                                kind) != 0)
            return -1;
        if (window->target)
            event->target_hit = 1;
        else
            event->flank_hit = 1;
    }
    return 0;
}

/* Worker-side CIGAR validation and segmentation.  This function only reads
 * immutable header/interval/config state and writes its private event. */
static int prepare_alignment_event(const run_context_t *context,
                                   alignment_event_t *event)
{
    const bam1_t *record = event->record;
    const bam1_core_t *core = &record->core;
    event->mapped = core->tid >= 0 && !(core->flag & BAM_FUNMAP);
    event->clean = !(core->flag & (BAM_FQCFAIL | BAM_FDUP)) &&
                   core->qual >= context->config->mapq;
    if (core->l_qseq < 0 || record->l_data < 0) {
        xerror("input contains a negative BAM/query length");
        return -1;
    }
    if (core->tid < 0 && !(core->flag & BAM_FUNMAP)) {
        xerror("mapped record has no reference id");
        return -1;
    }
    if (!event->mapped)
        return 0;
    if ((size_t)core->tid >= context->intervals->nchromosomes) {
        xerror("input contains a reference id outside its header");
        return -1;
    }
    if (core->pos < 0) {
        xerror("mapped record has a negative position");
        return -1;
    }
    chromosome_intervals_t *chromosome = &context->intervals->chromosomes[core->tid];
    if ((uint64_t)core->pos > chromosome->length || core->n_cigar == 0 ||
        record->data == NULL) {
        xerror("mapped record on %s has an invalid CIGAR/position", chromosome->name);
        return -1;
    }
    size_t cigar_bytes;
    if (size_mul((size_t)core->n_cigar, sizeof(uint32_t), &cigar_bytes) != 0 ||
        (size_t)core->l_qname > (size_t)record->l_data ||
        cigar_bytes > (size_t)record->l_data - (size_t)core->l_qname) {
        xerror("mapped record on %s has a truncated CIGAR", chromosome->name);
        return -1;
    }
    if (validate_qname(record, chromosome->name) != 0)
        return -1;
    uint32_t *cigar = bam_get_cigar(record);
    uint64_t reference_position = (uint64_t)core->pos;
    uint64_t query_position = 0;
    size_t capacity = 0;
    for (uint32_t i = 0; i < core->n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        uint32_t length = bam_cigar_oplen(cigar[i]);
        if (length == 0) {
            xerror("mapped record on %s contains a zero-length CIGAR operation",
                   chromosome->name);
            return -1;
        }
        int consumes_query = (op == BAM_CMATCH || op == BAM_CINS ||
                              op == BAM_CSOFT_CLIP || op == BAM_CEQUAL ||
                              op == BAM_CDIFF);
        if (consumes_query) {
            if ((uint64_t)length > UINT64_MAX - query_position ||
                query_position + (uint64_t)length > (uint64_t)core->l_qseq) {
                xerror("CIGAR query length exceeds sequence length on %s",
                       chromosome->name);
                return -1;
            }
            query_position += (uint64_t)length;
        }
        int consumes_reference = (op == BAM_CMATCH || op == BAM_CEQUAL ||
                                  op == BAM_CDIFF || op == BAM_CDEL ||
                                  op == BAM_CREF_SKIP);
        if (!consumes_reference) {
            if (op != BAM_CINS && op != BAM_CSOFT_CLIP && op != BAM_CHARD_CLIP &&
                op != BAM_CPAD) {
                xerror("unsupported CIGAR operation %d", op);
                return -1;
            }
            continue;
        }
        if ((uint64_t)length > UINT64_MAX - reference_position) {
            xerror("CIGAR reference position overflow");
            return -1;
        }
        uint64_t end = reference_position + (uint64_t)length;
        if (end > chromosome->length) {
            xerror("CIGAR of record on %s extends beyond chromosome length",
                   chromosome->name);
            return -1;
        }
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
            unsigned kind = (op == BAM_CDEL) ? 2U : 3U;
            if (context->config->fragment_mode &&
                fragment_add_segment(&event->segments, &event->nsegments, &capacity,
                                     reference_position, end, kind) != 0)
                return -1;
            if (prepare_alignment_deltas(chromosome, reference_position, end, kind,
                                         event) != 0)
                return -1;
        }
        reference_position = end;
    }
    if (query_position != (uint64_t)core->l_qseq) {
        xerror("CIGAR query length does not match sequence length on %s",
               chromosome->name);
        return -1;
    }
    return 0;
}

static int apply_alignment_event(run_context_t *context, sam_hdr_t *header,
                                 alignment_event_t *event, htsFile *bamout)
{
    const bam1_t *record = event->record;
    const bam1_core_t *core = &record->core;
    read_stats_t *stats = &context->result->reads;
    if (increment_counter(&stats->n_reads, "read") ||
        add_counter(&stats->n_data, (uint64_t)core->l_qseq, "raw data"))
        return -1;
    int qc_fail = (core->flag & BAM_FQCFAIL) != 0;
    int duplicate = (core->flag & BAM_FDUP) != 0;
    if (qc_fail && increment_counter(&stats->n_qcfail, "QC-fail read"))
        return -1;
    if (core->flag & BAM_FPAIRED) {
        int mapped = event->mapped;
        if (increment_counter(&stats->n_pair_all, "paired-read") ||
            ((core->flag & BAM_FPROPER_PAIR) && mapped && !(core->flag & BAM_FMUNMAP) &&
             increment_counter(&stats->n_pair_good, "proper-pair")) ||
            ((core->flag & BAM_FREAD1) && increment_counter(&stats->n_read1, "read1")) ||
            ((core->flag & BAM_FREAD2) && increment_counter(&stats->n_read2, "read2")))
            return -1;
        if (mapped && (core->flag & BAM_FMUNMAP) &&
            increment_counter(&stats->n_sgltn, "singleton"))
            return -1;
        if (mapped && !(core->flag & BAM_FMUNMAP) &&
            (increment_counter(&stats->n_pair_map, "mapped pair") ||
             (core->mtid != core->tid &&
              increment_counter(&stats->n_diffchr, "different-reference pair"))))
            return -1;
    }
    if (!event->mapped)
        return 0;
    if (increment_counter(&stats->n_mapped, "mapped read") ||
        add_counter(&stats->n_mdata, (uint64_t)core->l_qseq, "mapped data") ||
        ((core->flag & BAM_FREVERSE) && increment_counter(&stats->n_mstrand, "reverse-strand read")) ||
        (!(core->flag & BAM_FREVERSE) && increment_counter(&stats->n_pstrand, "forward-strand read")) ||
        (duplicate && increment_counter(&stats->n_dup, "duplicate read")))
        return -1;
    int high_mapq = core->qual >= context->config->mapq;
    if (high_mapq && increment_counter(&stats->n_qual, "mapQ-qualified read"))
        return -1;
    int clean = !qc_fail && !duplicate && high_mapq;
    if (clean && (((core->flag & BAM_FREAD1) &&
                   increment_counter(&stats->n_rmdup1, "rmdup read1")) ||
                  ((core->flag & BAM_FREAD2) &&
                   increment_counter(&stats->n_rmdup2, "rmdup read2"))))
        return -1;
    if (core->isize > 0 && core->isize < context->config->isize)
        histogram_add(context->result->insert_sizes, (uint64_t)core->isize, 1);

    for (size_t i = 0; i < event->ndeltas; ++i) {
        const alignment_delta_t *delta = &event->deltas[i];
        int match = (delta->kind & 1U) != 0;
        if (update_region_delta(delta->region, delta->start, delta->end,
                                match ? 1 : 0,
                                match && clean ? 1 : 0,
                                (delta->kind & 2U) ? 1 : 0,
                                context->config->summary_only) < 0)
            return -1;
    }
    if (event->target_hit && increment_counter(&stats->n_tgt, "target read"))
        return -1;
    if (event->flank_hit && increment_counter(&stats->n_flk, "flank read"))
        return -1;
    if (event->target_hit && bamout != NULL && sam_write1(bamout, header, record) < 0) {
        xerror("failed to write target BAM record");
        return -1;
    }
    const char *qname = NULL;
    const char *rg = NULL;
    if (context->config->fragment_mode) {
        qname = bam_get_qname(record);
        rg = record_rg(record);
    }
    int fragment_status = fragment_maybe_store(context, event->input_index, core,
                                               qname, rg, clean,
                                               event->segments, event->nsegments);
    if (fragment_status < 0)
        return -1;
    if (fragment_status == 1) {
        event->segments = NULL;
        event->nsegments = 0;
    }
    return 0;
}

typedef struct {
    pthread_t *threads;
    size_t nthreads;
    pthread_mutex_t mutex;
    pthread_cond_t work_ready;
    pthread_cond_t work_done;
    alignment_event_t *events;
    size_t event_count;
    size_t next_event;
    size_t completed;
    int stopping;
    const run_context_t *context;
} compute_pool_t;

static void *compute_worker(void *opaque)
{
    compute_pool_t *pool = (compute_pool_t *)opaque;
    for (;;) {
        pthread_mutex_lock(&pool->mutex);
        while (!pool->stopping && pool->next_event >= pool->event_count)
            pthread_cond_wait(&pool->work_ready, &pool->mutex);
        if (pool->stopping) {
            pthread_mutex_unlock(&pool->mutex);
            return NULL;
        }
        size_t index = pool->next_event++;
        alignment_event_t *event = &pool->events[index];
        pthread_mutex_unlock(&pool->mutex);
        event->status = prepare_alignment_event(pool->context, event);
        pthread_mutex_lock(&pool->mutex);
        ++pool->completed;
        if (pool->completed == pool->event_count)
            pthread_cond_signal(&pool->work_done);
        pthread_mutex_unlock(&pool->mutex);
    }
}

static int compute_pool_init(compute_pool_t *pool, size_t nthreads,
                             const run_context_t *context)
{
    memset(pool, 0, sizeof(*pool));
    pool->context = context;
    pool->nthreads = nthreads;
    pool->threads = xcalloc(nthreads, sizeof(*pool->threads));
    int mutex_initialized = pthread_mutex_init(&pool->mutex, NULL) == 0;
    int ready_initialized = mutex_initialized &&
                            pthread_cond_init(&pool->work_ready, NULL) == 0;
    int done_initialized = ready_initialized &&
                           pthread_cond_init(&pool->work_done, NULL) == 0;
    if (!mutex_initialized || !ready_initialized || !done_initialized) {
        xerror("failed to initialize compute worker synchronization");
        if (done_initialized) pthread_cond_destroy(&pool->work_done);
        if (ready_initialized) pthread_cond_destroy(&pool->work_ready);
        if (mutex_initialized) pthread_mutex_destroy(&pool->mutex);
        free(pool->threads);
        memset(pool, 0, sizeof(*pool));
        return -1;
    }
    for (size_t i = 0; i < nthreads; ++i) {
        if (pthread_create(&pool->threads[i], NULL, compute_worker, pool) != 0) {
            xerror("failed to create compute worker");
            pthread_mutex_lock(&pool->mutex);
            pool->stopping = 1;
            pthread_cond_broadcast(&pool->work_ready);
            pthread_mutex_unlock(&pool->mutex);
            for (size_t j = 0; j < i; ++j)
                pthread_join(pool->threads[j], NULL);
            pthread_cond_destroy(&pool->work_ready);
            pthread_cond_destroy(&pool->work_done);
            pthread_mutex_destroy(&pool->mutex);
            free(pool->threads);
            memset(pool, 0, sizeof(*pool));
            return -1;
        }
    }
    return 0;
}

static int compute_pool_run(compute_pool_t *pool, alignment_event_t *events,
                            size_t count)
{
    pthread_mutex_lock(&pool->mutex);
    pool->events = events;
    pool->event_count = count;
    pool->next_event = 0;
    pool->completed = 0;
    pthread_cond_broadcast(&pool->work_ready);
    while (pool->completed < count)
        pthread_cond_wait(&pool->work_done, &pool->mutex);
    pthread_mutex_unlock(&pool->mutex);
    for (size_t i = 0; i < count; ++i)
        if (events[i].status != 0)
            return -1;
    return 0;
}

static void compute_pool_destroy(compute_pool_t *pool)
{
    if (pool == NULL || pool->threads == NULL)
        return;
    pthread_mutex_lock(&pool->mutex);
    pool->stopping = 1;
    pthread_cond_broadcast(&pool->work_ready);
    pthread_mutex_unlock(&pool->mutex);
    for (size_t i = 0; i < pool->nthreads; ++i)
        pthread_join(pool->threads[i], NULL);
    pthread_cond_destroy(&pool->work_ready);
    pthread_cond_destroy(&pool->work_done);
    pthread_mutex_destroy(&pool->mutex);
    free(pool->threads);
    memset(pool, 0, sizeof(*pool));
}

static void alignment_event_destroy(alignment_event_t *event)
{
    if (event == NULL)
        return;
    bam_destroy1(event->record);
    free(event->deltas);
    free(event->segments);
    memset(event, 0, sizeof(*event));
}

static int advance_cursor(stream_cursor_t *cursor, sam_hdr_t *header,
                          size_t *heap, size_t *heap_size,
                          const stream_cursor_t *cursors)
{
    if (read_next(cursor, header) != 0)
        return -1;
    if (cursor->has_record)
        heap_push(heap, heap_size, cursor->source_index, cursors);
    return 0;
}

static int analysis_run_parallel(const xamdst_config_t *config,
                                 const input_set_t *inputs,
                                 interval_set_t *intervals,
                                 analysis_result_t *result,
                                 analysis_sink_t *sink, htsFile *bamout)
{
    run_context_t context = {
        .config = config,
        .intervals = intervals,
        .result = result,
        .sink = sink,
    };
    stream_cursor_t *cursors = xcalloc(inputs->count, sizeof(*cursors));
    size_t *heap = xcalloc(inputs->count == 0 ? 1 : inputs->count, sizeof(*heap));
    size_t heap_size = 0;
    int error = 0;
    for (size_t i = 0; i < inputs->count; ++i) {
        cursors[i].file = inputs->files[i];
        cursors[i].source_index = i;
        cursors[i].record = bam_init1();
        if (cursors[i].record == NULL || read_next(&cursors[i], inputs->header) != 0) {
            xerror("failed to initialize input cursor %zu", i + 1);
            error = 1;
            break;
        }
        if (cursors[i].has_record)
            heap_push(heap, &heap_size, i, cursors);
    }

    size_t batch_size = (size_t)config->compute_threads * 256U;
    if (batch_size < 1024U)
        batch_size = 1024U;
    if (batch_size > 8192U)
        batch_size = 8192U;
    alignment_event_t *events = xcalloc(batch_size, sizeof(*events));
    compute_pool_t pool;
    memset(&pool, 0, sizeof(pool));
    if (!error && compute_pool_init(&pool, (size_t)config->compute_threads, &context) != 0)
        error = 1;
    int active_tid = -1;
    while (!error && heap_size > 0) {
        size_t cursor_index = heap_pop(heap, &heap_size, cursors);
        stream_cursor_t *cursor = &cursors[cursor_index];
        bam1_t *record = cursor->record;
        int tid = record->core.tid;
        int unmapped = tid < 0 || (record->core.flag & BAM_FUNMAP) != 0;
        if (tid < -1 || (tid >= 0 && (size_t)tid >= intervals->nchromosomes)) {
            xerror("input contains an invalid reference id %d", tid);
            error = 1;
            break;
        }
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            if (advance_cursor(cursor, inputs->header, heap, &heap_size, cursors) != 0)
                error = 1;
            continue;
        }
        if (unmapped) {
            if (process_alignment(&context, inputs->header, record, cursor_index, bamout) != 0 ||
                advance_cursor(cursor, inputs->header, heap, &heap_size, cursors) != 0)
                error = 1;
            continue;
        }
        if (active_tid < 0) {
            for (int missing = 0; missing < tid && !error; ++missing)
                if (finish_chromosome(&context, (size_t)missing) != 0)
                    error = 1;
            if (!error)
                active_tid = tid;
        } else if (tid != active_tid) {
            fragment_clear(&context);
            if (finish_chromosome(&context, (size_t)active_tid) != 0)
                error = 1;
            for (int missing = active_tid + 1; missing < tid && !error; ++missing)
                if (finish_chromosome(&context, (size_t)missing) != 0)
                    error = 1;
            if (!error)
                active_tid = tid;
        }
        if (error)
            break;

        size_t count = 0;
        for (;;) {
            /* Transfer the populated bam1_t to the worker event and give the
             * cursor a fresh container for subsequent streaming reads.  This
             * avoids copying the BAM payload for every batch item. */
            events[count].record = record;
            events[count].input_index = cursor_index;
            ++count;
            cursor->record = bam_init1();
            if (cursor->record == NULL) {
                xerror("failed to allocate BAM cursor record for compute worker");
                error = 1;
                break;
            }
            if (advance_cursor(cursor, inputs->header, heap, &heap_size, cursors) != 0) {
                error = 1;
                break;
            }
            if (count == batch_size || heap_size == 0)
                break;
            stream_cursor_t *next_cursor = &cursors[heap[0]];
            bam1_t *next_record = next_cursor->record;
            int next_tid = next_record->core.tid;
            if ((next_record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) ||
                next_tid < 0 || (next_record->core.flag & BAM_FUNMAP) ||
                next_tid != active_tid)
                break;
            cursor_index = heap_pop(heap, &heap_size, cursors);
            cursor = &cursors[cursor_index];
            record = cursor->record;
        }
        if (!error && compute_pool_run(&pool, events, count) != 0)
            error = 1;
        for (size_t i = 0; i < count; ++i) {
            if (!error && events[i].status == 0 &&
                apply_alignment_event(&context, inputs->header, &events[i], bamout) != 0)
                error = 1;
            alignment_event_destroy(&events[i]);
        }
    }
    if (!error) {
        fragment_clear(&context);
        if (active_tid >= 0) {
            if (finish_chromosome(&context, (size_t)active_tid) != 0)
                error = 1;
            for (int missing = active_tid + 1; missing < (int)intervals->nchromosomes && !error; ++missing)
                if (finish_chromosome(&context, (size_t)missing) != 0)
                    error = 1;
        } else {
            for (size_t missing = 0; missing < intervals->nchromosomes && !error; ++missing)
                if (finish_chromosome(&context, missing) != 0)
                    error = 1;
        }
    }
    for (size_t i = 0; i < batch_size; ++i)
        alignment_event_destroy(&events[i]);
    fragment_clear(&context);
    compute_pool_destroy(&pool);
    for (size_t i = 0; i < inputs->count; ++i)
        bam_destroy1(cursors[i].record);
    free(events);
    free(heap);
    free(cursors);
    free(context.pending);
    return error ? -1 : 0;
}

int analysis_run(const xamdst_config_t *config, const input_set_t *inputs,
                 interval_set_t *intervals, analysis_result_t *result,
                 analysis_sink_t *sink, htsFile *bamout)
{
    if (sink == NULL || sink->depth_row == NULL || sink->region_row == NULL ||
        sink->uncovered_row == NULL) {
        xerror("analysis output sink is incomplete");
        return -1;
    }
    if (config->compute_threads > 1)
        return analysis_run_parallel(config, inputs, intervals, result, sink, bamout);
    run_context_t context = {
        .config = config,
        .intervals = intervals,
        .result = result,
        .sink = sink,
    };
    stream_cursor_t *cursors = xcalloc(inputs->count, sizeof(*cursors));
    size_t *heap = xcalloc(inputs->count == 0 ? 1 : inputs->count, sizeof(*heap));
    size_t heap_size = 0;
    int error = 0;
    for (size_t i = 0; i < inputs->count; ++i) {
        cursors[i].file = inputs->files[i];
        cursors[i].source_index = i;
        cursors[i].record = bam_init1();
        if (cursors[i].record == NULL) {
            xerror("failed to allocate BAM record");
            error = 1;
            break;
        }
        if (read_next(&cursors[i], inputs->header)) {
            error = 1;
            break;
        }
        if (cursors[i].has_record)
            heap_push(heap, &heap_size, i, cursors);
    }

    int active_tid = -1;
    while (!error && heap_size > 0) {
        size_t cursor_index = heap_pop(heap, &heap_size, cursors);
        stream_cursor_t *cursor = &cursors[cursor_index];
        bam1_t *record = cursor->record;
        int tid = record->core.tid;
        if (tid < -1 || (tid >= 0 && (size_t)tid >= intervals->nchromosomes)) {
            xerror("input contains an invalid reference id %d", tid);
            error = 1;
            break;
        }
        if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            if (read_next(cursor, inputs->header))
                error = 1;
            if (!error && cursor->has_record)
                heap_push(heap, &heap_size, cursor_index, cursors);
            continue;
        }
        if (!(record->core.flag & BAM_FUNMAP) && tid >= 0) {
            if (active_tid < 0) {
                for (int missing = 0; missing < tid; ++missing)
                    if (finish_chromosome(&context, (size_t)missing)) { error = 1; break; }
                if (!error) {
                    active_tid = tid;
                }
            } else if (tid != active_tid) {
                fragment_clear(&context);
                if (finish_chromosome(&context, (size_t)active_tid)) {
                    error = 1;
                } else {
                    for (int missing = active_tid + 1; missing < tid; ++missing)
                        if (finish_chromosome(&context, (size_t)missing)) { error = 1; break; }
                    if (!error) {
                        active_tid = tid;
                    }
                }
            }
            if (!error && process_alignment(&context, inputs->header, record,
                                            cursor_index, bamout))
                error = 1;
        } else if (!error) {
            /* Primary unmapped reads contribute to read totals but not coverage. */
            if (process_alignment(&context, inputs->header, record,
                                  cursor_index, bamout))
                error = 1;
        }
        if (!error && read_next(cursor, inputs->header))
            error = 1;
        if (!error && cursor->has_record)
            heap_push(heap, &heap_size, cursor_index, cursors);
    }
    if (!error) {
        fragment_clear(&context);
        if (active_tid >= 0) {
            if (finish_chromosome(&context, (size_t)active_tid))
                error = 1;
            for (int missing = active_tid + 1; !error && missing < (int)intervals->nchromosomes; ++missing)
                if (finish_chromosome(&context, (size_t)missing)) error = 1;
        } else {
            for (size_t missing = 0; missing < intervals->nchromosomes; ++missing)
                if (finish_chromosome(&context, missing)) { error = 1; break; }
        }
    }
    fragment_clear(&context);
    for (size_t i = 0; i < inputs->count; ++i)
        bam_destroy1(cursors[i].record);
    free(heap);
    free(cursors);
    free(context.pending);
    return error ? -1 : 0;
}
