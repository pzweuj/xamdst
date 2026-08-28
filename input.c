#include "input.h"

#include "util.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <htslib/thread_pool.h>

static int headers_compatible(const sam_hdr_t *expected, const sam_hdr_t *actual,
                              const char *path)
{
    int expected_count = sam_hdr_nref(expected);
    int actual_count = sam_hdr_nref(actual);
    if (expected_count < 0 || actual_count < 0) {
        xerror("input '%s' has an invalid reference dictionary", path);
        return -1;
    }
    if (expected_count != actual_count) {
        xerror("input '%s' has %d references; expected %d", path, actual_count, expected_count);
        return -1;
    }
    for (int i = 0; i < expected_count; ++i) {
        const char *expected_name = sam_hdr_tid2name(expected, i);
        const char *actual_name = sam_hdr_tid2name(actual, i);
        hts_pos_t expected_length = sam_hdr_tid2len(expected, i);
        hts_pos_t actual_length = sam_hdr_tid2len(actual, i);
        if (expected_name == NULL || actual_name == NULL ||
            strcmp(expected_name, actual_name) != 0 || expected_length != actual_length) {
            xerror("input '%s' has an incompatible reference dictionary at index %d", path, i);
            return -1;
        }
    }
    return 0;
}

static int open_one(htsFile **file, sam_hdr_t **header, const char *path,
                    const xamdst_config_t *config, htsThreadPool *thread_pool)
{
    *file = hts_open(path, "r");
    if (*file == NULL) {
        xerror("cannot open input '%s'", path);
        return -1;
    }
    if (thread_pool != NULL && hts_set_thread_pool(*file, thread_pool) != 0)
        xwarn("HTSlib could not enable the shared I/O thread pool for '%s'", path);

    const htsFormat *format = hts_get_format(*file);
    if (format != NULL && format->format == cram && config->reference == NULL) {
        xerror("CRAM input '%s' requires --reference/-T", path);
        hts_close(*file);
        *file = NULL;
        return -1;
    }
    if (config->reference != NULL && hts_set_fai_filename(*file, config->reference) != 0) {
        xerror("cannot set reference '%s' for input '%s'", config->reference, path);
        hts_close(*file);
        *file = NULL;
        return -1;
    }
    if (config->reference != NULL && access(config->reference, R_OK) != 0) {
        xerror("cannot read reference '%s': %s", config->reference, strerror(errno));
        hts_close(*file);
        *file = NULL;
        return -1;
    }
    *header = sam_hdr_read(*file);
    if (*header == NULL) {
        xerror("cannot read SAM header from '%s'", path);
        hts_close(*file);
        *file = NULL;
        return -1;
    }
    /* For a stream (notably CRAM on stdin) HTSlib cannot always identify the
     * format until the header has been consumed.  Re-check after
     * sam_hdr_read() so every CRAM input, including "-", obeys the explicit
     * reference-file contract. */
    format = hts_get_format(*file);
    if (format != NULL && format->format == cram && config->reference == NULL) {
        xerror("CRAM input '%s' requires --reference/-T", path);
        sam_hdr_destroy(*header);
        *header = NULL;
        hts_close(*file);
        *file = NULL;
        return -1;
    }
    return 0;
}

int inputs_open(input_set_t *inputs, const xamdst_config_t *config)
{
    memset(inputs, 0, sizeof(*inputs));
    inputs->count = config->ninputs;
    inputs->files = xcalloc(inputs->count, sizeof(*inputs->files));
    if (config->threads > 0) {
        /* HTSlib creates one native worker per requested thread.  Bound the
         * allocation defensively for accidentally huge command-line values;
         * the CLI remains backwards compatible with non-negative integers. */
        int thread_count = config->threads > 256 ? 256 : config->threads;
        if (thread_count != config->threads)
            xwarn("capping HTSlib I/O threads at %d", thread_count);
        inputs->thread_pool.pool = hts_tpool_init(thread_count);
        if (inputs->thread_pool.pool == NULL) {
            xerror("failed to initialize %d HTSlib I/O threads", config->threads);
            inputs_close(inputs);
            return -1;
        }
        inputs->thread_pool.qsize = thread_count > INT_MAX / 2 ? INT_MAX : thread_count * 2;
        inputs->thread_pool_active = 1;
    }
    for (size_t i = 0; i < inputs->count; ++i) {
        sam_hdr_t *header = NULL;
        if (open_one(&inputs->files[i], &header, config->inputs[i], config,
                     inputs->thread_pool_active ? &inputs->thread_pool : NULL)) {
            sam_hdr_destroy(header);
            inputs_close(inputs);
            return -1;
        }
        if (i == 0) {
            inputs->header = header;
        } else {
            int compatible = headers_compatible(inputs->header, header, config->inputs[i]);
            sam_hdr_destroy(header);
            if (compatible) {
                inputs_close(inputs);
                return -1;
            }
        }
    }
    /* Open output streams use the very same pool; xamdst keeps inputs alive
     * until BAMout has been closed. */
    return 0;
}

int inputs_close(input_set_t *inputs)
{
    if (inputs == NULL)
        return 0;
    int status = 0;
    if (inputs->files != NULL) {
        for (size_t i = 0; i < inputs->count; ++i)
            if (inputs->files[i] != NULL) {
                if (hts_close(inputs->files[i]) < 0)
                    status = -1;
            }
        free(inputs->files);
    }
    if (inputs->thread_pool_active) {
        hts_tpool_destroy(inputs->thread_pool.pool);
        inputs->thread_pool.pool = NULL;
        inputs->thread_pool_active = 0;
    }
    sam_hdr_destroy(inputs->header);
    memset(inputs, 0, sizeof(*inputs));
    return status;
}
