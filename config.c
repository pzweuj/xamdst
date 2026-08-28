#include "config.h"

#include "util.h"

#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

enum {
    OPT_MAXDEPTH = 1000,
    OPT_CUTOFF,
    OPT_ISIZE,
    OPT_UNCOVER,
    OPT_BAMOUT,
    OPT_DEPTHRATIO,
    OPT_THREADS,
    OPT_COMPUTE_THREADS,
    OPT_FRAGMENT_MODE,
    OPT_SUMMARY_ONLY
};

static const struct option long_options[] = {
    {"outdir", required_argument, NULL, 'o'},
    {"bed", required_argument, NULL, 'p'},
    {"flank", required_argument, NULL, 'f'},
    {"mapthres", required_argument, NULL, 'q'},
    {"reference", required_argument, NULL, 'T'},
    {"maxdepth", required_argument, NULL, OPT_MAXDEPTH},
    {"cutoffdepth", required_argument, NULL, OPT_CUTOFF},
    {"isize", required_argument, NULL, OPT_ISIZE},
    {"uncover", required_argument, NULL, OPT_UNCOVER},
    {"bamout", required_argument, NULL, OPT_BAMOUT},
    {"depthratio", required_argument, NULL, OPT_DEPTHRATIO},
    {"threads", required_argument, NULL, OPT_THREADS},
    {"compute-threads", required_argument, NULL, OPT_COMPUTE_THREADS},
    {"fragment-mode", no_argument, NULL, OPT_FRAGMENT_MODE},
    {"summary-only", no_argument, NULL, OPT_SUMMARY_ONLY},
    {"help", no_argument, NULL, 'h'},
    {"version", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}
};

void config_init(xamdst_config_t *config)
{
    memset(config, 0, sizeof(*config));
    config->flank = 200;
    config->mapq = 20;
    config->isize = 2000;
    config->uncover = 5;
    config->compute_threads = 1;
    config->nratios = 2;
    config->ratios[0] = 0.2;
    config->ratios[1] = 0.5;
}

void config_destroy(xamdst_config_t *config)
{
    free(config->bed_path);
    free(config->outdir);
    free(config->reference);
    free(config->bamout);
    if (config->inputs != NULL) {
        for (size_t i = 0; i < config->ninputs; ++i)
            free(config->inputs[i]);
        free(config->inputs);
    }
    memset(config, 0, sizeof(*config));
}

void config_usage(const char *program, int full)
{
    FILE *out = full ? stdout : stderr;
    fprintf(out,
            "xamdst %s\n"
            "Usage: %s [OPTION] -p probe.bed -o OUTPUT [FILE ...]\n"
            "       %s [OPTION] -p probe.bed -o OUTPUT -\n\n",
            XAMDST_VERSION, program, program);
    fprintf(out, "Required:\n"
                 "  -p, --bed FILE          target regions (BED, merged before counting)\n"
                 "  -o, --outdir DIR        output directory (created if needed)\n\n");
    fprintf(out, "Options:\n"
                 "  -T, --reference FILE    reference FASTA for CRAM\n"
                 "  -f, --flank N           outside flank size (default: 200)\n"
                 "  -q, --mapthres N        mapQ cutoff (default: 20)\n"
                 "      --maxdepth N        cap rows in cumu.plot (default: 0)\n"
                 "      --cutoffdepth LIST  comma-separated >= depth cutoffs (max 10)\n"
                 "      --depthratio LIST   comma-separated ratios in (0,1] (max 10)\n"
                 "      --isize N           insert-size limit (default: 2000)\n"
                 "      --uncover N         raw-depth uncover threshold, depth < N (default: 5)\n"
                 "      --bamout FILE       write target reads to BAM\n"
                 "      --threads N         HTSlib I/O threads (default: 0)\n"
                 "      --compute-threads N coverage workers (default: 1, max: 256)\n"
                 "      --fragment-mode     remove overlapping bases of matched read pairs\n"
                 "      --summary-only      omit per-base depth.tsv.gz output\n"
                 "  -1                       input BED is 1-based inclusive\n"
                 "  -h, --help               show this help\n"
                 "  -v, --version            show version\n\n"
                 "Input files must use the same coordinate-sorted header.\n"
                 "Outputs: coverage.report, coverage.report.json, cumu.plot, insert.plot,\n"
                 "         chromosome.report, region.tsv.gz, depth.tsv.gz, uncover.bed\n");
}

static int parse_long(const char *text, long min, long max, long *result, const char *name)
{
    if (text == NULL || *text == '\0' || isspace((unsigned char)*text)) {
        xerror("%s requires an integer", name);
        return -1;
    }
    errno = 0;
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' || value < min || value > max) {
        xerror("invalid %s: '%s' (expected %ld..%ld)", name, text, min, max);
        return -1;
    }
    *result = value;
    return 0;
}

static int parse_double(const char *text, double *result, const char *name)
{
    if (text == NULL || *text == '\0' || isspace((unsigned char)*text)) {
        xerror("%s requires a number", name);
        return -1;
    }
    errno = 0;
    char *end = NULL;
    double value = strtod(text, &end);
    if (errno == ERANGE || end == text || *end != '\0' || !isfinite(value) ||
        value <= 0.0 || value > 1.0) {
        xerror("invalid %s: '%s' (expected a value in (0,1])", name, text);
        return -1;
    }
    *result = value;
    return 0;
}

static int parse_cutoff_list(xamdst_config_t *config, const char *text)
{
    size_t text_length = text != NULL ? strlen(text) : 0;
    if (text_length == 0 || text[0] == ',' || text[text_length - 1] == ',' ||
        strstr(text, ",,") != NULL) {
        xerror("--cutoffdepth must contain non-empty comma-separated values");
        return -1;
    }
    char *copy = xstrdup(text);
    char *save = NULL;
    char *token = strtok_r(copy, ",", &save);
    config->ncutoffs = 0;
    if (token == NULL) {
        free(copy);
        xerror("--cutoffdepth cannot be empty");
        return -1;
    }
    while (token != NULL) {
        if (config->ncutoffs == XAMDST_MAX_CUSTOM) {
            free(copy);
            xerror("--cutoffdepth accepts at most %d values", XAMDST_MAX_CUSTOM);
            return -1;
        }
        long value;
        if (parse_long(token, 1, INT_MAX, &value, "cutoff depth")) {
            free(copy);
            return -1;
        }
        config->cutoffs[config->ncutoffs++] = (int)value;
        token = strtok_r(NULL, ",", &save);
    }
    free(copy);
    return 0;
}

static int parse_ratio_list(xamdst_config_t *config, const char *text)
{
    size_t text_length = text != NULL ? strlen(text) : 0;
    if (text_length == 0 || text[0] == ',' || text[text_length - 1] == ',' ||
        strstr(text, ",,") != NULL) {
        xerror("--depthratio must contain non-empty comma-separated values");
        return -1;
    }
    char *copy = xstrdup(text);
    char *save = NULL;
    char *token = strtok_r(copy, ",", &save);
    config->nratios = 0;
    if (token == NULL) {
        free(copy);
        xerror("--depthratio cannot be empty");
        return -1;
    }
    while (token != NULL) {
        if (config->nratios == XAMDST_MAX_CUSTOM) {
            free(copy);
            xerror("--depthratio accepts at most %d values", XAMDST_MAX_CUSTOM);
            return -1;
        }
        double value;
        if (parse_double(token, &value, "depth ratio")) {
            free(copy);
            return -1;
        }
        config->ratios[config->nratios++] = value;
        token = strtok_r(NULL, ",", &save);
    }
    free(copy);
    return 0;
}

static int has_bam_extension(const char *path)
{
    size_t length = strlen(path);
    return length >= 4 && strcasecmp(path + length - 4, ".bam") == 0;
}

int config_parse(xamdst_config_t *config, int argc, char **argv)
{
    int option;
    int option_index = 0;
    const char *program = argc > 0 && argv[0] != NULL ? argv[0] : "xamdst";
    optind = 1;
    opterr = 0;
    while ((option = getopt_long(argc, argv, "o:p:f:q:T:h1v", long_options, &option_index)) != -1) {
        long value;
        switch (option) {
        case 'o':
            free(config->outdir);
            config->outdir = xstrdup(optarg);
            break;
        case 'p':
            free(config->bed_path);
            config->bed_path = xstrdup(optarg);
            break;
        case 'f':
            if (parse_long(optarg, 0, INT_MAX, &value, "flank")) return -1;
            config->flank = (int)value;
            break;
        case 'q':
            if (parse_long(optarg, 0, 255, &value, "mapQ")) return -1;
            config->mapq = (int)value;
            break;
        case 'T':
            if (optarg == NULL || optarg[0] == '\0') {
                xerror("--reference/-T requires a non-empty path");
                return -1;
            }
            free(config->reference);
            config->reference = xstrdup(optarg);
            break;
        case OPT_MAXDEPTH:
            if (parse_long(optarg, 0, INT_MAX, &value, "maxdepth")) return -1;
            config->maxdepth = (int)value;
            break;
        case OPT_CUTOFF:
            if (parse_cutoff_list(config, optarg)) return -1;
            break;
        case OPT_ISIZE:
            if (parse_long(optarg, 1, INT_MAX, &value, "isize")) return -1;
            config->isize = (int)value;
            break;
        case OPT_UNCOVER:
            if (parse_long(optarg, 1, INT_MAX, &value, "uncover")) return -1;
            config->uncover = (int)value;
            break;
        case OPT_BAMOUT:
            free(config->bamout);
            config->bamout = xstrdup(optarg);
            break;
        case OPT_DEPTHRATIO:
            if (parse_ratio_list(config, optarg)) return -1;
            break;
        case OPT_THREADS:
            if (parse_long(optarg, 0, INT_MAX, &value, "threads")) return -1;
            config->threads = (int)value;
            break;
        case OPT_COMPUTE_THREADS:
            if (parse_long(optarg, 1, XAMDST_MAX_COMPUTE_THREADS, &value,
                           "compute threads")) return -1;
            config->compute_threads = (int)value;
            break;
        case OPT_FRAGMENT_MODE:
            config->fragment_mode = 1;
            break;
        case OPT_SUMMARY_ONLY:
            config->summary_only = 1;
            break;
        case '1':
            config->one_based = 1;
            break;
        case 'h':
            config_usage(program, 1);
            return 1;
        case 'v':
            printf("%s\n", XAMDST_VERSION);
            return 1;
        case '?':
        default:
            if (optopt != 0)
                xerror("unknown or incomplete option '-%c'", optopt);
            else
                xerror("unknown option near '%s'", argv[optind - 1]);
            config_usage(program, 0);
            return -1;
        }
    }

    if (config->bed_path == NULL || config->bed_path[0] == '\0') {
        xerror("--bed/-p is required");
        config_usage(program, 0);
        return -1;
    }
    if (config->outdir == NULL || config->outdir[0] == '\0') {
        xerror("--outdir/-o is required");
        config_usage(program, 0);
        return -1;
    }
    if (config->bamout != NULL && !has_bam_extension(config->bamout)) {
        xerror("--bamout must have a .bam extension: %s", config->bamout);
        return -1;
    }

    config->ninputs = (size_t)(argc - optind);
    if (config->ninputs == 0) {
        config->ninputs = 1;
        config->inputs = xcalloc(1, sizeof(*config->inputs));
        config->inputs[0] = xstrdup("-");
    } else {
        config->inputs = xcalloc(config->ninputs, sizeof(*config->inputs));
        int stdin_count = 0;
        for (size_t i = 0; i < config->ninputs; ++i) {
            config->inputs[i] = xstrdup(argv[(size_t)optind + i]);
            if (strcmp(config->inputs[i], "-") == 0)
                ++stdin_count;
        }
        if (stdin_count > 1) {
            xerror("stdin ('-') may appear at most once");
            return -1;
        }
    }
    return 0;
}
