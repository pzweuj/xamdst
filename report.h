#ifndef XAMDST_REPORT_H
#define XAMDST_REPORT_H

#include <stddef.h>

#include <stdio.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "config.h"
#include "engine.h"
#include "intervals.h"

#define REPORT_DNA_OUTPUTS 8
#define REPORT_MAX_OUTPUTS 10

/* All report file names.  The first REPORT_DNA_OUTPUTS form the DNA-mode
 * set; RNA mode appends splice.tsv.gz and distribution.tsv. */
extern const char *const report_output_names[REPORT_MAX_OUTPUTS];
size_t report_output_count(const xamdst_config_t *config);

typedef struct report_writer {
    char *final_paths[REPORT_MAX_OUTPUTS];
    char *temporary_paths[REPORT_MAX_OUTPUTS];
    char *backup_paths[REPORT_MAX_OUTPUTS];
    int temporary_created[REPORT_MAX_OUTPUTS];
    int backup_created[REPORT_MAX_OUTPUTS];
    size_t output_count;
    int depth_enabled;
    BGZF *depth;
    BGZF *region;
    FILE *uncovered;
    BGZF *splice;
    FILE *distribution;
    int depth_closed;
    int region_closed;
    int uncovered_closed;
    int splice_closed;
    int distribution_closed;
    int open;
    kstring_t depth_buffer;
    kstring_t region_buffer;
} report_writer_t;

int report_open(report_writer_t *writer, const char *outdir,
                const xamdst_config_t *config);
analysis_sink_t report_sink(report_writer_t *writer);
int report_finish(report_writer_t *writer, const xamdst_config_t *config,
                  const interval_set_t *intervals, const analysis_result_t *result);
int report_commit(report_writer_t *writer);
void report_abort(report_writer_t *writer);

#endif
