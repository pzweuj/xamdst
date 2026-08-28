#ifndef XAMDST_INPUT_H
#define XAMDST_INPUT_H

#include <stddef.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "config.h"

typedef struct {
    htsFile **files;
    sam_hdr_t *header;
    size_t count;
    htsThreadPool thread_pool;
    int thread_pool_active;
} input_set_t;

int inputs_open(input_set_t *inputs, const xamdst_config_t *config);
int inputs_close(input_set_t *inputs);

#endif
