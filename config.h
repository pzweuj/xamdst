#ifndef XAMDST_CONFIG_H
#define XAMDST_CONFIG_H

#include <stddef.h>

#define XAMDST_VERSION "3.1.0"
#define XAMDST_MAX_CUSTOM 10
#define XAMDST_MAX_COMPUTE_THREADS 256

typedef struct {
    char *bed_path;
    char *outdir;
    char *reference;
    char *bamout;
    char **inputs;
    size_t ninputs;
    int flank;
    int mapq;
    int maxdepth;
    int isize;
    int uncover;
    int threads;
    int compute_threads;
    int one_based;
    int fragment_mode;
    int summary_only;
    size_t ncutoffs;
    int cutoffs[XAMDST_MAX_CUSTOM];
    size_t nratios;
    double ratios[XAMDST_MAX_CUSTOM];
} xamdst_config_t;

void config_init(xamdst_config_t *config);
void config_destroy(xamdst_config_t *config);
void config_usage(const char *program, int full);

/* Returns 0 for a valid configuration, 1 for --help/--version, -1 on error. */
int config_parse(xamdst_config_t *config, int argc, char **argv);

#endif
