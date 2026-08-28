#include "intervals.h"

#include "util.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

typedef struct {
    uint64_t start;
    uint64_t end;
} raw_interval_t;

typedef struct {
    raw_interval_t *items;
    size_t count;
    size_t capacity;
} raw_vec_t;

static void raw_push(raw_vec_t *vec, uint64_t start, uint64_t end)
{
    if (vec->count == vec->capacity) {
        size_t next = vec->capacity == 0 ? 16 : vec->capacity * 2;
        if (next < vec->capacity) {
            xerror("too many BED intervals");
            exit(EXIT_FAILURE);
        }
        vec->items = xreallocarray(vec->items, next, sizeof(*vec->items));
        vec->capacity = next;
    }
    vec->items[vec->count++] = (raw_interval_t){start, end};
}

static void region_push(region_vec_t *vec, uint64_t start, uint64_t end)
{
    if (vec->count == vec->capacity) {
        size_t next = vec->capacity == 0 ? 16 : vec->capacity * 2;
        if (next < vec->capacity) {
            xerror("too many BED regions");
            exit(EXIT_FAILURE);
        }
        vec->items = xreallocarray(vec->items, next, sizeof(*vec->items));
        vec->capacity = next;
    }
    vec->items[vec->count++] = (depth_region_t){
        .start = start,
        .end = end,
        .length = end - start,
    };
}

static int raw_compare(const void *lhs, const void *rhs)
{
    const raw_interval_t *a = (const raw_interval_t *)lhs;
    const raw_interval_t *b = (const raw_interval_t *)rhs;
    if (a->start < b->start) return -1;
    if (a->start > b->start) return 1;
    if (a->end < b->end) return -1;
    if (a->end > b->end) return 1;
    return 0;
}

static void merge_raw(raw_vec_t *vec)
{
    if (vec->count == 0)
        return;
    qsort(vec->items, vec->count, sizeof(*vec->items), raw_compare);
    size_t out = 0;
    for (size_t i = 0; i < vec->count; ++i) {
        raw_interval_t current = vec->items[i];
        if (out == 0 || current.start > vec->items[out - 1].end) {
            vec->items[out++] = current;
        } else if (current.end > vec->items[out - 1].end) {
            vec->items[out - 1].end = current.end;
        }
    }
    vec->count = out;
}

static void raw_to_regions(const raw_vec_t *raw, region_vec_t *regions)
{
    for (size_t i = 0; i < raw->count; ++i)
        region_push(regions, raw->items[i].start, raw->items[i].end);
}

static int build_coverage_windows(chromosome_intervals_t *chrom)
{
    if (chrom->flank.count > SIZE_MAX - chrom->target.count) {
        xerror("too many coverage windows on chromosome '%s'", chrom->name);
        return -1;
    }
    size_t total = chrom->target.count + chrom->flank.count;
    if (total == 0)
        return 0;
    chrom->windows = xcalloc(total, sizeof(*chrom->windows));
    size_t target = 0;
    size_t flank = 0;
    while (target < chrom->target.count || flank < chrom->flank.count) {
        int take_target;
        if (flank == chrom->flank.count)
            take_target = 1;
        else if (target == chrom->target.count)
            take_target = 0;
        else
            take_target = chrom->target.items[target].start <= chrom->flank.items[flank].start;
        coverage_window_t *window = &chrom->windows[chrom->nwindows++];
        window->target = take_target;
        window->region = take_target ? &chrom->target.items[target++] : &chrom->flank.items[flank++];
    }
    return 0;
}

static int read_gz_line(gzFile file, char **buffer, size_t *capacity, size_t *length)
{
    *length = 0;
    for (;;) {
        int value = gzgetc(file);
        if (value == -1) {
            int error_number = Z_OK;
            (void)gzerror(file, &error_number);
            if (gzeof(file) && (error_number == Z_OK || error_number == Z_STREAM_END)) {
                if (*length == 0)
                    return 0;
                if (*capacity == 0) {
                    *buffer = xreallocarray(*buffer, 1, sizeof(**buffer));
                    *capacity = 1;
                }
                if (*length > 0 && (*buffer)[*length - 1] == '\r')
                    --*length;
                (*buffer)[*length] = '\0';
                return 1;
            }
            xerror("failed while reading BED file (zlib error %d)", error_number);
            return -1;
        }
        if (value == '\n')
            break;
        if (value == '\0') {
            xerror("BED file contains a NUL byte");
            return -1;
        }
        if (*length == SIZE_MAX || *length + 1 >= *capacity) {
            if (*length == SIZE_MAX) {
                xerror("BED line is too long");
                return -1;
            }
            size_t next = *capacity == 0 ? 4096 : *capacity * 2;
            if (next < *capacity) {
                xerror("BED line is too long");
                return -1;
            }
            *buffer = xreallocarray(*buffer, next, sizeof(**buffer));
            *capacity = next;
        }
        (*buffer)[(*length)++] = (char)value;
    }
    if (*length > 0 && (*buffer)[*length - 1] == '\r')
        --*length;
    if (*capacity == 0) {
        *buffer = xreallocarray(*buffer, 1, sizeof(**buffer));
        *capacity = 1;
    }
    (*buffer)[*length] = '\0';
    return 1;
}

static int parse_u64_token(const char *token, uint64_t *value)
{
    /* strtoull() accepts a leading minus and may otherwise turn it into a
     * large unsigned value, so reject signs explicitly before conversion. */
    if (token == NULL || *token == '\0' || token[0] == '-' || token[0] == '+')
        return -1;
    errno = 0;
    char *end = NULL;
    unsigned long long parsed = strtoull(token, &end, 10);
    if (errno == ERANGE || end == token || *end != '\0')
        return -1;
    *value = (uint64_t)parsed;
    return 0;
}

static size_t chromosome_index(const interval_set_t *set, const char *name)
{
    for (size_t i = 0; i < set->nchromosomes; ++i)
        if (strcmp(set->chromosomes[i].name, name) == 0)
            return i;
    return SIZE_MAX;
}

static void build_flank_for_chromosome(chromosome_intervals_t *chrom, int flank)
{
    raw_vec_t expanded = {0};
    raw_vec_t pieces = {0};
    uint64_t flank_size = (uint64_t)flank;
    for (size_t i = 0; i < chrom->target.count; ++i) {
        depth_region_t *target = &chrom->target.items[i];
        uint64_t start = target->start > flank_size ? target->start - flank_size : 0;
        uint64_t end = target->end;
        if (flank_size > chrom->length - end)
            end = chrom->length;
        else
            end += flank_size;
        if (start < end)
            raw_push(&expanded, start, end);
    }

    /* Subtract the complete merged target union from each expanded window. */
    for (size_t i = 0; i < expanded.count; ++i) {
        uint64_t cursor = expanded.items[i].start;
        uint64_t end = expanded.items[i].end;
        /* Expanded windows can overlap when flank is large.  Locate the
         * first target independently for each window; carrying a monotonic
         * cursor across windows would skip targets covered by a later
         * expansion and incorrectly put them back into the flank union. */
        size_t first_target = 0;
        size_t target_high = chrom->target.count;
        while (first_target < target_high) {
            size_t middle = first_target + (target_high - first_target) / 2;
            if (chrom->target.items[middle].end <= cursor)
                first_target = middle + 1;
            else
                target_high = middle;
        }
        for (size_t j = first_target; j < chrom->target.count && cursor < end; ++j) {
            const depth_region_t *target = &chrom->target.items[j];
            if (target->end <= cursor)
                continue;
            if (target->start >= end)
                break;
            if (target->start > cursor)
                raw_push(&pieces, cursor, target->start < end ? target->start : end);
            if (target->end > cursor)
                cursor = target->end;
        }
        if (cursor < end)
            raw_push(&pieces, cursor, end);
    }
    merge_raw(&pieces);
    raw_to_regions(&pieces, &chrom->flank);
    free(expanded.items);
    free(pieces.items);
}

int intervals_load(interval_set_t *set, const char *path, const sam_hdr_t *header,
                   int one_based, int flank)
{
    memset(set, 0, sizeof(*set));
    if (flank < 0) {
        xerror("flank cannot be negative");
        return -1;
    }
    int nref = sam_hdr_nref(header);
    if (nref < 0) {
        xerror("invalid SAM header: cannot enumerate references");
        return -1;
    }
    set->nchromosomes = (size_t)nref;
    set->chromosomes = xcalloc(set->nchromosomes, sizeof(*set->chromosomes));
    for (size_t i = 0; i < set->nchromosomes; ++i) {
        set->chromosomes[i].tid = (int)i;
        const char *name = sam_hdr_tid2name(header, (int)i);
        hts_pos_t length = sam_hdr_tid2len(header, (int)i);
        if (name == NULL || length < 0) {
            xerror("invalid SAM header reference at index %zu", i);
            intervals_destroy(set);
            return -1;
        }
        set->chromosomes[i].name = xstrdup(name);
        set->chromosomes[i].length = (uint64_t)length;
    }

    gzFile file = gzopen(path, "rb");
    if (file == NULL) {
        xerror("cannot open BED file '%s': %s", path, strerror(errno));
        intervals_destroy(set);
        return -1;
    }
    raw_vec_t *raw = xcalloc(set->nchromosomes, sizeof(*raw));
    char *line = NULL;
    size_t capacity = 0;
    size_t line_number = 0;
    size_t line_length = 0;
    int status;
    while ((status = read_gz_line(file, &line, &capacity, &line_length)) > 0) {
        ++line_number;
        char *cursor = line;
        while (isspace((unsigned char)*cursor)) ++cursor;
        if (*cursor == '\0' || *cursor == '#')
            continue;
        char *tokens[3] = {NULL, NULL, NULL};
        size_t token_count = 0;
        while (*cursor != '\0' && token_count < 3) {
            while (isspace((unsigned char)*cursor)) ++cursor;
            if (*cursor == '\0' || *cursor == '#') break;
            tokens[token_count++] = cursor;
            while (*cursor != '\0' && !isspace((unsigned char)*cursor)) ++cursor;
            if (*cursor != '\0') *cursor++ = '\0';
        }
        if (token_count < 3) {
            xerror("BED line %zu must contain chromosome, start and end", line_number);
            status = -1;
            break;
        }
        size_t chromosome = chromosome_index(set, tokens[0]);
        if (chromosome == SIZE_MAX) {
            xerror("BED line %zu references unknown chromosome '%s'", line_number, tokens[0]);
            status = -1;
            break;
        }
        uint64_t begin, end;
        if (parse_u64_token(tokens[1], &begin) || parse_u64_token(tokens[2], &end)) {
            xerror("BED line %zu has invalid coordinates", line_number);
            status = -1;
            break;
        }
        if (one_based) {
            if (begin == 0 || end < begin) {
                xerror("BED line %zu has invalid 1-based interval", line_number);
                status = -1;
                break;
            }
            --begin;
        } else if (end <= begin) {
            xerror("BED line %zu must satisfy end > start", line_number);
            status = -1;
            break;
        }
        if (end > set->chromosomes[chromosome].length) {
            xerror("BED line %zu exceeds chromosome '%s' length %llu", line_number,
                   tokens[0], (unsigned long long)set->chromosomes[chromosome].length);
            status = -1;
            break;
        }
        raw_push(&raw[chromosome], begin, end);
    }
    free(line);
    int close_status = gzclose(file);
    if (status >= 0 && close_status != Z_OK) {
        xerror("failed while closing BED file '%s' (zlib error %d)", path, close_status);
        status = -1;
    }
    if (status < 0) {
        for (size_t i = 0; i < set->nchromosomes; ++i) free(raw[i].items);
        free(raw);
        intervals_destroy(set);
        return -1;
    }
    for (size_t i = 0; i < set->nchromosomes; ++i) {
        merge_raw(&raw[i]);
        raw_to_regions(&raw[i], &set->chromosomes[i].target);
        for (size_t j = 0; j < set->chromosomes[i].target.count; ++j) {
            uint64_t length = set->chromosomes[i].target.items[j].length;
            if (u64_add(set->target_bases, length, &set->target_bases)) {
                xerror("target length overflow");
                for (size_t k = 0; k < set->nchromosomes; ++k) free(raw[k].items);
                free(raw);
                intervals_destroy(set);
                return -1;
            }
        }
        if (set->chromosomes[i].target.count > 0)
            build_flank_for_chromosome(&set->chromosomes[i], flank);
        if (build_coverage_windows(&set->chromosomes[i]) != 0) {
            for (size_t k = 0; k < set->nchromosomes; ++k) free(raw[k].items);
            free(raw);
            intervals_destroy(set);
            return -1;
        }
        for (size_t j = 0; j < set->chromosomes[i].flank.count; ++j)
            if (u64_add(set->flank_bases, set->chromosomes[i].flank.items[j].length,
                        &set->flank_bases)) {
                xerror("flank length overflow");
                for (size_t k = 0; k < set->nchromosomes; ++k) free(raw[k].items);
                free(raw);
                intervals_destroy(set);
                return -1;
            }
    }
    for (size_t i = 0; i < set->nchromosomes; ++i)
        free(raw[i].items);
    free(raw);
    return 0;
}

int region_ensure_buffers(depth_region_t *region, int summary_only)
{
    if (region == NULL)
        return -1;
    if (region->raw_diff != NULL && region->rmdup_diff != NULL &&
        region->cov_diff != NULL && (!summary_only || region->dirty_bits != NULL))
        return 0;

    /* A region should normally have either all three arrays or none.  Be
     * defensive if a caller hands us a partially initialized object: discard
     * the partial state before rebuilding it, so update_regions never writes
     * through a stale or mismatched buffer. */
    free(region->raw_diff);
    free(region->rmdup_diff);
    free(region->cov_diff);
    free(region->dirty_bits);
    region->raw_diff = NULL;
    region->rmdup_diff = NULL;
    region->cov_diff = NULL;
    region->dirty_bits = NULL;
    region->dirty_words = 0;

    /* A diff array has one extra sentinel slot for the end coordinate.  Do
     * all validation before allocating so a malformed/oversized interval can
     * never wrap the size conversion.  xcalloc() terminates on allocation
     * failure; the explicit checks below cover arithmetic overflow. */
    if (region->length > (uint64_t)SIZE_MAX - 1) {
        xerror("region length is too large to allocate safely");
        return -1;
    }
    size_t slots = (size_t)region->length + 1;
    if (slots > SIZE_MAX / sizeof(*region->raw_diff)) {
        xerror("region length is too large to allocate safely");
        return -1;
    }
    size_t dirty_words = 0;
    if (summary_only) {
        if (slots > SIZE_MAX - 63) {
            xerror("region length is too large to allocate safely");
            return -1;
        }
        dirty_words = (slots + 63) / 64;
        if (dirty_words > SIZE_MAX / sizeof(uint64_t)) {
            xerror("region length is too large to allocate safely");
            return -1;
        }
    }

    int64_t *raw = xcalloc(slots, sizeof(*raw));
    int64_t *rmdup = xcalloc(slots, sizeof(*rmdup));
    int64_t *coverage = xcalloc(slots, sizeof(*coverage));
    uint64_t *dirty = summary_only ? xcalloc(dirty_words, sizeof(*dirty)) : NULL;
    region->raw_diff = raw;
    region->rmdup_diff = rmdup;
    region->cov_diff = coverage;
    region->dirty_bits = dirty;
    region->dirty_words = dirty_words;
    return 0;
}

void regions_release_buffers(region_vec_t *regions)
{
    for (size_t i = 0; i < regions->count; ++i) {
        free(regions->items[i].raw_diff);
        free(regions->items[i].rmdup_diff);
        free(regions->items[i].cov_diff);
        free(regions->items[i].dirty_bits);
        regions->items[i].raw_diff = NULL;
        regions->items[i].rmdup_diff = NULL;
        regions->items[i].cov_diff = NULL;
        regions->items[i].dirty_bits = NULL;
        regions->items[i].dirty_words = 0;
    }
}

void intervals_destroy(interval_set_t *set)
{
    if (set == NULL)
        return;
    for (size_t i = 0; i < set->nchromosomes; ++i) {
        free(set->chromosomes[i].name);
        regions_release_buffers(&set->chromosomes[i].target);
        regions_release_buffers(&set->chromosomes[i].flank);
        free(set->chromosomes[i].windows);
        free(set->chromosomes[i].target.items);
        free(set->chromosomes[i].flank.items);
    }
    free(set->chromosomes);
    memset(set, 0, sizeof(*set));
}
