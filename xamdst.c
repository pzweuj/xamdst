#include "config.h"
#include "engine.h"
#include "input.h"
#include "intervals.h"
#include "report.h"
#include "util.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#include <htslib/sam.h>

static char *canonical_path(const char *path)
{
    char *resolved = realpath(path, NULL);
    if (resolved != NULL)
        return resolved;

    /* The final report may not exist yet.  Resolve its parent when possible so
     * equivalent relative paths (for example ./out/../out/file) still match. */
    const char *slash = strrchr(path, '/');
    const char *basename = slash != NULL ? slash + 1 : path;
    char *parent = NULL;
    if (slash == NULL) {
        parent = xstrdup(".");
    } else if (slash == path) {
        parent = xstrdup("/");
    } else {
        size_t length = (size_t)(slash - path);
        parent = xmalloc(length + 1);
        memcpy(parent, path, length);
        parent[length] = '\0';
    }
    char *resolved_parent = realpath(parent, NULL);
    if (resolved_parent != NULL) {
        char *result = path_join(resolved_parent, basename);
        free(resolved_parent);
        free(parent);
        return result;
    }
    free(parent);

    char cwd[PATH_MAX];
    if (path[0] != '/' && getcwd(cwd, sizeof(cwd)) != NULL)
        return path_join(cwd, path);
    return xstrdup(path);
}

static int paths_equal(const char *left, const char *right)
{
    struct stat left_stat;
    struct stat right_stat;
    if (stat(left, &left_stat) == 0 && stat(right, &right_stat) == 0 &&
        left_stat.st_dev == right_stat.st_dev && left_stat.st_ino == right_stat.st_ino)
        return 1;
    char *left_canonical = canonical_path(left);
    char *right_canonical = canonical_path(right);
    int equal = strcmp(left_canonical, right_canonical) == 0;
    free(left_canonical);
    free(right_canonical);
    return equal;
}

static int path_equals_output(const xamdst_config_t *config, const char *path)
{
    static const char *const names[] = {
        "coverage.report", "coverage.report.json", "cumu.plot", "insert.plot",
        "chromosome.report", "region.tsv.gz", "depth.tsv.gz", "uncover.bed"
    };
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        char *output = path_join(config->outdir, names[i]);
        int equal = paths_equal(output, path);
        free(output);
        if (equal)
            return 1;
    }
    return 0;
}

static int path_equals_input(const xamdst_config_t *config, const char *path)
{
    for (size_t i = 0; i < config->ninputs; ++i) {
        if (strcmp(config->inputs[i], "-") == 0)
            continue;
        if (paths_equal(config->inputs[i], path))
            return 1;
    }
    return 0;
}

static int path_equals_protected_file(const xamdst_config_t *config, const char *path)
{
    if (path_equals_input(config, path))
        return 1;
    if (config->bed_path != NULL && paths_equal(config->bed_path, path))
        return 1;
    if (config->reference != NULL && paths_equal(config->reference, path))
        return 1;
    return 0;
}

static int report_path_conflicts(const xamdst_config_t *config)
{
    static const char *const names[] = {
        "coverage.report", "coverage.report.json", "cumu.plot", "insert.plot",
        "chromosome.report", "region.tsv.gz", "depth.tsv.gz", "uncover.bed"
    };
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        char *output = path_join(config->outdir, names[i]);
        int conflict = path_equals_protected_file(config, output);
        free(output);
        if (conflict) {
            xerror("output report '%s' would overwrite a protected input file", names[i]);
            return -1;
        }
    }
    return 0;
}

typedef struct {
    char *backup;
    const char *final;
    int moved_old;
    int installed;
} bamout_transaction_t;

static void rollback_bamout(bamout_transaction_t *transaction);

static int stage_bamout(bamout_transaction_t *transaction,
                        const char *temporary, const char *final)
{
    memset(transaction, 0, sizeof(*transaction));
    transaction->final = final;
    struct stat existing;
    if (stat(final, &existing) == 0 && !S_ISREG(existing.st_mode)) {
        xerror("target BAM path is not a regular file: '%s'", final);
        return -1;
    }
    if (stat(final, &existing) != 0 && errno != ENOENT) {
        xerror("cannot inspect target BAM '%s': %s", final, strerror(errno));
        return -1;
    }
    char *backup = temporary_path(final, 12000U);
    if (access(backup, F_OK) == 0) {
        xerror("stale target BAM backup exists '%s'", backup);
        free(backup);
        return -1;
    }
    int moved_old = 0;
    if (rename(final, backup) == 0) {
        moved_old = 1;
    } else if (errno != ENOENT) {
        xerror("cannot stage existing target BAM '%s': %s", final, strerror(errno));
        free(backup);
        return -1;
    }
    transaction->backup = backup;
    transaction->moved_old = moved_old;
    if (rename(temporary, final) != 0) {
        int saved_errno = errno;
        rollback_bamout(transaction);
        xerror("cannot install target BAM '%s': %s", final, strerror(saved_errno));
        return -1;
    }
    transaction->installed = 1;
    return 0;
}

static void rollback_bamout(bamout_transaction_t *transaction)
{
    if (transaction == NULL)
        return;
    if (transaction->installed && transaction->final != NULL) {
        if (remove(transaction->final) != 0 && errno != ENOENT)
            xwarn("cannot remove failed target BAM '%s': %s", transaction->final,
                  strerror(errno));
    }
    if (transaction->moved_old && transaction->backup != NULL) {
        if (rename(transaction->backup, transaction->final) != 0) {
            xerror("cannot restore previous target BAM '%s': %s", transaction->final,
                   strerror(errno));
        } else {
            transaction->moved_old = 0;
        }
    }
    if (transaction->backup != NULL && transaction->moved_old == 0)
        (void)remove(transaction->backup);
    free(transaction->backup);
    memset(transaction, 0, sizeof(*transaction));
}

static void commit_bamout(bamout_transaction_t *transaction)
{
    if (transaction == NULL)
        return;
    if (transaction->moved_old && transaction->backup != NULL &&
        remove(transaction->backup) != 0) {
        xwarn("cannot remove target BAM backup '%s': %s", transaction->backup,
              strerror(errno));
    }
    free(transaction->backup);
    memset(transaction, 0, sizeof(*transaction));
}

static int reserve_path(const char *path)
{
    int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0600);
    if (fd < 0)
        return -1;
    if (close(fd) != 0) {
        int saved_errno = errno;
        (void)remove(path);
        errno = saved_errno;
        return -1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    xamdst_config_t config;
    config_init(&config);
    int parse_status = config_parse(&config, argc, argv);
    if (parse_status != 0) {
        config_destroy(&config);
        return parse_status > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    /* Create/validate the output directory before comparing paths.  This
     * lets realpath() resolve the report parent and reliably catches aliases
     * such as ./out/../out/depth.tsv.gz before any input is opened. */
    if (mkdir_p(config.outdir, 0755) != 0) {
        xerror("cannot create output directory '%s': %s", config.outdir,
               strerror(errno));
        config_destroy(&config);
        return EXIT_FAILURE;
    }
    if (report_path_conflicts(&config) != 0) {
        config_destroy(&config);
        return EXIT_FAILURE;
    }
    if (config.bamout != NULL && path_equals_output(&config, config.bamout)) {
        xerror("--bamout must not overwrite one of the report outputs");
        config_destroy(&config);
        return EXIT_FAILURE;
    }
    if (config.bamout != NULL && path_equals_protected_file(&config, config.bamout)) {
        xerror("--bamout must not overwrite an input, BED or reference file");
        config_destroy(&config);
        return EXIT_FAILURE;
    }

    input_set_t inputs;
    if (inputs_open(&inputs, &config)) {
        config_destroy(&config);
        return EXIT_FAILURE;
    }
    interval_set_t intervals;
    if (intervals_load(&intervals, config.bed_path, inputs.header,
                       config.one_based, config.flank)) {
        inputs_close(&inputs);
        config_destroy(&config);
        return EXIT_FAILURE;
    }

    analysis_result_t result;
    analysis_result_init(&result, &intervals);
    report_writer_t writer;
    if (report_open(&writer, config.outdir, &config)) {
        analysis_result_destroy(&result);
        intervals_destroy(&intervals);
        inputs_close(&inputs);
        config_destroy(&config);
        return EXIT_FAILURE;
    }
    analysis_sink_t sink = report_sink(&writer);

    htsFile *bamout = NULL;
    char *bamout_temp = NULL;
    bamout_transaction_t bamout_transaction = {0};
    if (config.bamout != NULL) {
        bamout_temp = temporary_path(config.bamout, 9000);
        if (reserve_path(bamout_temp) != 0) {
            xerror("cannot reserve temporary target BAM '%s': %s", bamout_temp,
                   strerror(errno));
            report_abort(&writer);
            free(bamout_temp);
            analysis_result_destroy(&result);
            intervals_destroy(&intervals);
            inputs_close(&inputs);
            config_destroy(&config);
            return EXIT_FAILURE;
        }
        bamout = hts_open(bamout_temp, "wb");
        if (bamout == NULL) {
            xerror("cannot create target BAM '%s': %s", config.bamout, strerror(errno));
            (void)remove(bamout_temp);
            report_abort(&writer);
            free(bamout_temp);
            analysis_result_destroy(&result);
            intervals_destroy(&intervals);
            inputs_close(&inputs);
            config_destroy(&config);
            return EXIT_FAILURE;
        }
        if (inputs.thread_pool_active) {
            if (hts_set_thread_pool(bamout, &inputs.thread_pool) != 0)
                xwarn("HTSlib could not enable the shared I/O thread pool for BAMout");
        } else if (config.threads > 0 && hts_set_threads(bamout, config.threads) != 0) {
            xwarn("HTSlib could not enable %d output threads", config.threads);
        }
        if (sam_hdr_write(bamout, inputs.header) < 0) {
            xerror("cannot write target BAM header");
            hts_close(bamout);
            remove(bamout_temp);
            free(bamout_temp);
            report_abort(&writer);
            analysis_result_destroy(&result);
            intervals_destroy(&intervals);
            inputs_close(&inputs);
            config_destroy(&config);
            return EXIT_FAILURE;
        }
    }

    int status = analysis_run(&config, &inputs, &intervals, &result, &sink, bamout);
    if (bamout != NULL && hts_close(bamout) < 0)
        status = -1;
    bamout = NULL;
    if (inputs_close(&inputs) != 0)
        status = -1;
    if (status == 0 && report_finish(&writer, &config, &intervals, &result) != 0)
        status = -1;
    if (status == 0 && bamout_temp != NULL &&
        stage_bamout(&bamout_transaction, bamout_temp, config.bamout) != 0)
        status = -1;
    if (status == 0 && report_commit(&writer) != 0)
        status = -1;
    if (status != 0) {
        report_abort(&writer);
        if (bamout_transaction.installed)
            rollback_bamout(&bamout_transaction);
        if (bamout_temp != NULL)
            remove(bamout_temp);
    } else if (bamout_transaction.installed) {
        commit_bamout(&bamout_transaction);
    }
    free(bamout_temp);
    analysis_result_destroy(&result);
    intervals_destroy(&intervals);
    config_destroy(&config);
    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
