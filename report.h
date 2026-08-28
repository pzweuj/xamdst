#ifndef XAMDST_REPORT_H
#define XAMDST_REPORT_H

#include <stddef.h>

#include <stdio.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "config.h"
#include "engine.h"
#include "intervals.h"

typedef struct report_writer {
    char *final_paths[8];
    char *temporary_paths[8];
    char *backup_paths[8];
    int temporary_created[8];
    int backup_created[8];
    int depth_enabled;
    BGZF *depth;
    BGZF *region;
    FILE *uncovered;
    int depth_closed;
    int region_closed;
    int uncovered_closed;
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
