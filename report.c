#include "report.h"

#include "util.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

enum { OUTPUT_COUNT = 8 };
enum { REPORT_BUFFER_SIZE = 64 * 1024 };

static const char *const output_names[OUTPUT_COUNT] = {
    "coverage.report",
    "coverage.report.json",
    "cumu.plot",
    "insert.plot",
    "chromosome.report",
    "region.tsv.gz",
    "depth.tsv.gz",
    "uncover.bed"
};

static int write_bgzf(BGZF *file, const char *text, size_t length)
{
    ssize_t written = bgzf_write(file, text, length);
    if (written < 0 || (size_t)written != length) {
        xerror("failed to write compressed output");
        return -1;
    }
    return 0;
}

static int flush_buffer(BGZF *file, kstring_t *buffer)
{
    if (buffer->l == 0)
        return 0;
    if (write_bgzf(file, buffer->s, buffer->l) != 0)
        return -1;
    buffer->l = 0;
    return 0;
}

static int depth_row(void *opaque, const char *chromosome, uint64_t position,
                     uint64_t raw, uint64_t rmdup, uint64_t coverage)
{
    report_writer_t *writer = (report_writer_t *)opaque;
    if (!writer->depth_enabled)
        return 0;
    if (ksprintf(&writer->depth_buffer,
                 "%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                 chromosome, position, raw, rmdup, coverage) < 0) {
        xerror("failed to format depth output row");
        return -1;
    }
    if (writer->depth_buffer.l >= REPORT_BUFFER_SIZE)
        return flush_buffer(writer->depth, &writer->depth_buffer);
    return 0;
}

static int region_row(void *opaque, const char *chromosome, const depth_region_t *region)
{
    report_writer_t *writer = (report_writer_t *)opaque;
    if (ksprintf(&writer->region_buffer,
                 "%s\t%" PRIu64 "\t%" PRIu64
                 "\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                 chromosome, region->start, region->end, region->mean, region->median,
                 region->length ? (double)region->raw_covered * 100.0 / (double)region->length : 0.0,
                 region->cov_mean, region->cov_median,
                 region->length ? (double)region->covered * 100.0 / (double)region->length : 0.0) < 0) {
        xerror("failed to format region output row");
        return -1;
    }
    if (writer->region_buffer.l >= REPORT_BUFFER_SIZE)
        return flush_buffer(writer->region, &writer->region_buffer);
    return 0;
}

static int uncovered_row(void *opaque, const char *chromosome, uint64_t start, uint64_t end)
{
    report_writer_t *writer = (report_writer_t *)opaque;
    if (fprintf(writer->uncovered, "%s\t%" PRIu64 "\t%" PRIu64 "\n",
                chromosome, start, end) < 0) {
        xerror("failed to write uncover.bed");
        return -1;
    }
    return 0;
}

int report_open(report_writer_t *writer, const char *outdir,
                const xamdst_config_t *config)
{
    memset(writer, 0, sizeof(*writer));
    writer->depth_enabled = config == NULL || !config->summary_only;
    if (mkdir_p(outdir, 0755) != 0) {
        xerror("cannot create output directory '%s': %s", outdir, strerror(errno));
        return -1;
    }
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        writer->final_paths[i] = path_join(outdir, output_names[i]);
        writer->temporary_paths[i] = temporary_path(writer->final_paths[i], (unsigned)i);
        if (access(writer->temporary_paths[i], F_OK) == 0) {
            xerror("stale temporary output exists '%s'", writer->temporary_paths[i]);
            report_abort(writer);
            return -1;
        }
        if (i == 6 && !writer->depth_enabled)
            continue;
        int fd = open(writer->temporary_paths[i], O_WRONLY | O_CREAT | O_EXCL, 0600);
        if (fd < 0) {
            xerror("cannot reserve temporary output '%s': %s",
                   writer->temporary_paths[i], strerror(errno));
            report_abort(writer);
            return -1;
        }
        if (close(fd) != 0) {
            xerror("cannot close temporary output '%s': %s",
                   writer->temporary_paths[i], strerror(errno));
            (void)remove(writer->temporary_paths[i]);
            report_abort(writer);
            return -1;
        }
        writer->temporary_created[i] = 1;
    }
    if (writer->depth_enabled)
        writer->depth = bgzf_open(writer->temporary_paths[6], "w");
    writer->region = bgzf_open(writer->temporary_paths[5], "w");
    writer->uncovered = fopen(writer->temporary_paths[7], "wb");
    if ((writer->depth_enabled && writer->depth == NULL) || writer->region == NULL || writer->uncovered == NULL) {
        xerror("cannot create compressed output files in '%s'", outdir);
        report_abort(writer);
        return -1;
    }
    (void)setvbuf(writer->uncovered, NULL, _IOFBF, REPORT_BUFFER_SIZE);
    const char *depth_header = "#Chr\tPos\tRaw depth\tRmdup depth\tCoverage with deletions\n";
    const char *region_header = "#Chr\tStart\tStop\tRaw mean\tRaw median\tRaw coverage\tCoverage with deletions mean\tCoverage with deletions median\tCoverage with deletions coverage\n";
    if ((writer->depth_enabled && write_bgzf(writer->depth, depth_header, strlen(depth_header))) ||
        write_bgzf(writer->region, region_header, strlen(region_header)) ||
        fprintf(writer->uncovered, "#Chr\tStart\tEnd\n") < 0) {
        report_abort(writer);
        return -1;
    }
    writer->open = 1;
    return 0;
}

analysis_sink_t report_sink(report_writer_t *writer)
{
    analysis_sink_t sink;
    sink.opaque = writer;
    sink.depth_row = depth_row;
    sink.region_row = region_row;
    sink.uncovered_row = uncovered_row;
    return sink;
}

static double percentage(uint64_t numerator, uint64_t denominator)
{
    return denominator == 0 ? 0.0 : (double)numerator * 100.0 / (double)denominator;
}

static int close_text_file(FILE *file, const char *path)
{
    int write_error = ferror(file);
    int close_error = fclose(file);
    if (write_error || close_error != 0) {
        xerror("failed to write '%s': %s", path,
               close_error != 0 ? strerror(errno) : "stream error");
        return -1;
    }
    return 0;
}

static uint64_t count_greater_than(const depth_histogram_t *histogram, double threshold)
{
    size_t count = 0;
    const histogram_pair_t *pairs = histogram_sorted(histogram, &count);
    uint64_t result = 0;
    for (size_t i = 0; i < count; ++i)
        if ((double)pairs[i].depth > threshold) {
            if (UINT64_MAX - result < pairs[i].count) {
                xerror("histogram count overflow");
                exit(EXIT_FAILURE);
            }
            result += pairs[i].count;
        }
    return result;
}

static int write_plot(const char *path, const depth_histogram_t *histogram, int maxdepth)
{
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        xerror("cannot create '%s': %s", path, strerror(errno));
        return -1;
    }
    (void)setvbuf(file, NULL, _IOFBF, REPORT_BUFFER_SIZE);
    size_t count = 0;
    const histogram_pair_t *pairs = histogram_sorted(histogram, &count);
    uint64_t total = histogram_total(histogram);
    uint64_t cumulative = total;
    for (size_t i = 0; i < count; ++i) {
        cumulative -= pairs[i].count;
        if (maxdepth > 0 && pairs[i].depth > (uint64_t)maxdepth)
            continue;
        if (fprintf(file, "%" PRIu64 "\t%" PRIu64 "\t%.9f\t%" PRIu64 "\t%.9f\n",
                    pairs[i].depth, pairs[i].count, percentage(pairs[i].count, total) / 100.0,
                    cumulative, percentage(cumulative, total) / 100.0) < 0) {
            xerror("failed to write '%s'", path);
            fclose(file);
            return -1;
        }
    }
    return close_text_file(file, path);
}

static int write_chromosome_report(const char *path, const xamdst_config_t *config,
                                   const analysis_result_t *result)
{
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        xerror("cannot create '%s': %s", path, strerror(errno));
        return -1;
    }
    (void)setvbuf(file, NULL, _IOFBF, REPORT_BUFFER_SIZE);
    fprintf(file, "#Chromosome\tRAW_DATA(%%)\tRaw avg depth\tRaw median\tRaw coverage%%\tRaw cov 4x %%\tRaw cov 10x %%\tRaw cov 30x %%\tRaw cov 100x %%\tCoverage with deletions avg depth\tCoverage with deletions median\tCoverage with deletions%%\tCoverage with deletions 4x %%\tCoverage with deletions 10x %%\tCoverage with deletions 30x %%\tCoverage with deletions 100x %%");
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "\tCov %dx %%", config->cutoffs[i]);
    fputc('\n', file);
    for (size_t i = 0; i < result->nchromosomes; ++i) {
        const chromosome_summary_t *summary = &result->chromosomes[i];
        if (summary->length == 0)
            continue;
        double average = summary->length ? (double)summary->data / (double)summary->length : 0.0;
        double coverage_average = summary->length ? (double)summary->coverage_data / (double)summary->length : 0.0;
        double data_fraction = percentage(summary->data, result->target_data);
        fprintf(file, "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f",
                summary->name, data_fraction, average, histogram_median(summary->depth),
                percentage(histogram_count_at_least(summary->depth, 1), summary->length),
                percentage(histogram_count_at_least(summary->depth, 4), summary->length),
                percentage(histogram_count_at_least(summary->depth, 10), summary->length),
                percentage(histogram_count_at_least(summary->depth, 30), summary->length),
                percentage(histogram_count_at_least(summary->depth, 100), summary->length),
                coverage_average, histogram_median(summary->coverage_depth),
                percentage(histogram_count_at_least(summary->coverage_depth, 1), summary->length),
                percentage(histogram_count_at_least(summary->coverage_depth, 4), summary->length),
                percentage(histogram_count_at_least(summary->coverage_depth, 10), summary->length),
                percentage(histogram_count_at_least(summary->coverage_depth, 30), summary->length),
                percentage(histogram_count_at_least(summary->coverage_depth, 100), summary->length));
        for (size_t j = 0; j < config->ncutoffs; ++j)
            fprintf(file, "\t%.6f", percentage(histogram_count_at_least(summary->depth,
                                                                        (uint64_t)config->cutoffs[j]),
                                                summary->length));
        fputc('\n', file);
    }
    return close_text_file(file, path);
}

static void json_escaped(FILE *file, const char *text)
{
    fputc('"', file);
    for (const unsigned char *p = (const unsigned char *)text; *p != '\0'; ++p) {
        switch (*p) {
        case '"': fputs("\\\"", file); break;
        case '\\': fputs("\\\\", file); break;
        case '\b': fputs("\\b", file); break;
        case '\f': fputs("\\f", file); break;
        case '\n': fputs("\\n", file); break;
        case '\r': fputs("\\r", file); break;
        case '\t': fputs("\\t", file); break;
        default:
            if (*p < 0x20) fprintf(file, "\\u%04x", (unsigned int)*p);
            else fputc(*p, file);
        }
    }
    fputc('"', file);
}

static void json_uint(FILE *file, const char *key, uint64_t value, int comma)
{
    fprintf(file, "    \"%s\": %" PRIu64 "%s\n", key, value, comma ? "," : "");
}

static void json_double(FILE *file, const char *key, double value, int comma)
{
    if (!isfinite(value)) value = 0.0;
    fprintf(file, "    \"%s\": %.6f%s\n", key, value, comma ? "," : "");
}

static void json_pct(FILE *file, const char *key, uint64_t numerator, uint64_t denominator, int comma)
{
    json_double(file, key, percentage(numerator, denominator), comma);
}

static void json_coverage(FILE *file, const char *indent, const depth_histogram_t *histogram,
                          uint64_t length, const xamdst_config_t *config, double average,
                          int include_ratio)
{
    uint64_t total = histogram_total(histogram);
    uint64_t denominator = length != 0 ? length : total;
    uint64_t ge0 = histogram_count_at_least(histogram, 1);
    uint64_t ge4 = histogram_count_at_least(histogram, 4);
    uint64_t ge10 = histogram_count_at_least(histogram, 10);
    uint64_t ge30 = histogram_count_at_least(histogram, 30);
    uint64_t ge100 = histogram_count_at_least(histogram, 100);
    fprintf(file, "%s\"gt_0x\": %.6f,\n", indent, percentage(ge0, denominator));
    fprintf(file, "%s\"gte_4x\": %.6f,\n", indent, percentage(ge4, denominator));
    fprintf(file, "%s\"gte_10x\": %.6f,\n", indent, percentage(ge10, denominator));
    fprintf(file, "%s\"gte_30x\": %.6f,\n", indent, percentage(ge30, denominator));
    fprintf(file, "%s\"gte_100x\": %.6f", indent, percentage(ge100, denominator));
    if (include_ratio) {
        /* These two keys are part of the historical JSON contract.  Keep them
         * even when the caller adds custom ratios below. */
        fprintf(file, ",\n%s\"gt_0_2_avg\": %.6f,\n%s\"gt_0_5_avg\": %.6f",
                indent, percentage(count_greater_than(histogram, 0.2 * average), denominator),
                indent, percentage(count_greater_than(histogram, 0.5 * average), denominator));
    }
    if (include_ratio && config->nratios > 0) {
        fprintf(file, ",\n%s\"depth_ratio\": [", indent);
        for (size_t i = 0; i < config->nratios; ++i) {
            if (i) fputs(", ", file);
            fprintf(file, "{\"ratio\": %.6f, \"coverage\": %.6f}",
                    config->ratios[i],
                    percentage(count_greater_than(histogram, config->ratios[i] * average), denominator));
        }
        fputs("],\n", file);
        fprintf(file, "%s\"ratio_based\": {", indent);
        for (size_t i = 0; i < config->nratios; ++i) {
            if (i) fputs(", ", file);
            fprintf(file, "\"gt_%.6g_avg\": %.6f", config->ratios[i],
                    percentage(count_greater_than(histogram, config->ratios[i] * average), denominator));
        }
        fputc('}', file);
    }
    if (config->ncutoffs > 0) {
        fprintf(file, ",\n%s\"custom\": {", indent);
        for (size_t i = 0; i < config->ncutoffs; ++i) {
            if (i) fputs(", ", file);
            fprintf(file, "\"gte_%dx\": %.6f", config->cutoffs[i],
                    percentage(histogram_count_at_least(histogram, (uint64_t)config->cutoffs[i]), denominator));
        }
        fputc('}', file);
    }
    fputc('\n', file);
}

static int write_json(const char *path, const xamdst_config_t *config,
                      const analysis_result_t *result)
{
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        xerror("cannot create '%s': %s", path, strerror(errno));
        return -1;
    }
    (void)setvbuf(file, NULL, _IOFBF, REPORT_BUFFER_SIZE);
    const read_stats_t *stats = &result->reads;
    double average = percentage(result->target_data, result->target_bases) / 100.0;
    double average_rmdup = percentage(result->target_rmdup_data, result->target_bases) / 100.0;
    double flank_average = percentage(result->flank_data, result->flank_bases) / 100.0;
    fprintf(file, "{\n  \"schema_version\": \"3.1\",\n  \"version\": \"%s\",\n"
                 "  \"depth_output\": %s,\n  \"files\": [",
            XAMDST_VERSION, config->summary_only ? "false" : "true");
    for (size_t i = 0; i < config->ninputs; ++i) {
        if (i) fputs(", ", file);
        json_escaped(file, config->inputs[i]);
    }
    fputs("],\n  \"total\": {\n", file);
    json_uint(file, "raw_reads", stats->n_reads, 1);
    json_uint(file, "qc_fail_reads", stats->n_qcfail, 1);
    json_double(file, "raw_data_mb", (double)stats->n_data / 1e6, 1);
    json_uint(file, "paired_reads", stats->n_pair_all, 1);
    json_uint(file, "mapped_reads", stats->n_mapped, 1);
    json_pct(file, "mapped_reads_fraction", stats->n_mapped, stats->n_reads, 1);
    json_double(file, "mapped_data_mb", (double)stats->n_mdata / 1e6, 1);
    json_pct(file, "mapped_data_fraction", stats->n_mdata, stats->n_data, 1);
    json_uint(file, "properly_paired", stats->n_pair_good, 1);
    json_pct(file, "properly_paired_fraction", stats->n_pair_good, stats->n_reads, 1);
    json_uint(file, "read_mate_paired", stats->n_pair_map, 1);
    json_pct(file, "read_mate_paired_fraction", stats->n_pair_map, stats->n_reads, 1);
    json_uint(file, "singletons", stats->n_sgltn, 1);
    json_uint(file, "diff_chr", stats->n_diffchr, 1);
    json_uint(file, "read1", stats->n_read1, 1);
    json_uint(file, "read2", stats->n_read2, 1);
    json_uint(file, "read1_rmdup", stats->n_rmdup1, 1);
    json_uint(file, "read2_rmdup", stats->n_rmdup2, 1);
    json_uint(file, "forward_strand", stats->n_pstrand, 1);
    json_uint(file, "backward_strand", stats->n_mstrand, 1);
    json_uint(file, "pcr_duplicates", stats->n_dup, 1);
    json_pct(file, "pcr_duplicates_fraction", stats->n_dup, stats->n_mapped, 1);
    fprintf(file, "    \"mapq_cutoff\": %d,\n", config->mapq);
    json_uint(file, "mapq_reads", stats->n_qual, 1);
    json_pct(file, "mapq_reads_fraction_all", stats->n_qual, stats->n_reads, 1);
    json_pct(file, "mapq_reads_fraction_mapped", stats->n_qual, stats->n_mapped, 0);
    fputs("\n  },\n  \"fragment\": {\n", file);
    fprintf(file, "    \"enabled\": %s,\n", config->fragment_mode ? "true" : "false");
    json_uint(file, "paired", result->fragments.paired, 1);
    json_uint(file, "clipped", result->fragments.clipped, 1);
    json_uint(file, "unmatched", result->fragments.unmatched, 1);
    json_uint(file, "ambiguous", result->fragments.ambiguous, 0);
    fputs("  },\n  \"insert_size\": {\n", file);
    json_double(file, "average", histogram_total(result->insert_sizes) ?
                (double)histogram_sum(result->insert_sizes) / histogram_total(result->insert_sizes) : 0.0, 1);
    json_double(file, "median", histogram_median(result->insert_sizes), 0);
    fputs("\n  },\n  \"target\": {\n", file);
    json_uint(file, "target_reads", stats->n_tgt, 1);
    json_pct(file, "target_reads_fraction_all", stats->n_tgt, stats->n_reads, 1);
    json_pct(file, "target_reads_fraction_mapped", stats->n_tgt, stats->n_mapped, 1);
    json_double(file, "target_data_mb", (double)result->target_data / 1e6, 1);
    json_double(file, "target_data_rmdup_mb", (double)result->target_rmdup_data / 1e6, 1);
    json_double(file, "target_data_with_deletions_mb",
                (double)result->target_coverage_data / 1e6, 1);
    json_pct(file, "target_data_fraction_all", result->target_data, stats->n_data, 1);
    json_pct(file, "target_data_fraction_mapped", result->target_data, stats->n_mdata, 1);
    json_pct(file, "target_data_with_deletions_fraction_all",
             result->target_coverage_data, stats->n_data, 1);
    json_pct(file, "target_data_with_deletions_fraction_mapped",
             result->target_coverage_data, stats->n_mdata, 1);
    json_uint(file, "region_length", result->target_bases, 1);
    json_double(file, "average_depth", average, 1);
    json_double(file, "average_depth_rmdup", average_rmdup, 1);
    json_double(file, "average_depth_with_deletions",
                result->target_bases ? (double)result->target_coverage_data /
                                      (double)result->target_bases : 0.0, 1);
    fputs("    \"coverage\": {\n", file);
    json_coverage(file, "      ", result->target_depth, result->target_bases, config, average, 1);
    fputs("    },\n    \"coverage_rmdup\": {\n", file);
    json_coverage(file, "      ", result->target_rmdup_depth, result->target_bases, config, average_rmdup, 0);
    fputs("    },\n    \"coverage_with_deletions\": {\n", file);
    json_coverage(file, "      ", result->target_coverage_depth, result->target_bases,
                  config, result->target_bases ? (double)result->target_coverage_data /
                  (double)result->target_bases : 0.0, 0);
    fputs("    },\n", file);
    json_uint(file, "region_count", result->target_region_count, 1);
    fputs("    \"region_coverage\": {\n", file);
    json_uint(file, "gt_0x_count", histogram_count_at_least(result->region_means, 1), 1);
    json_pct(file, "gt_0x_fraction", histogram_count_at_least(result->region_means, 1),
             histogram_total(result->region_means), 1);
    json_pct(file, "gte_4x_fraction", histogram_count_at_least(result->region_means, 4),
             histogram_total(result->region_means), 1);
    json_pct(file, "gte_10x_fraction", histogram_count_at_least(result->region_means, 10),
             histogram_total(result->region_means), 1);
    json_pct(file, "gte_30x_fraction", histogram_count_at_least(result->region_means, 30),
             histogram_total(result->region_means), 1);
    json_pct(file, "gte_100x_fraction", histogram_count_at_least(result->region_means, 100),
             histogram_total(result->region_means), 0);
    fputs("\n    }\n  },\n  \"flank\": {\n", file);
    fprintf(file, "    \"flank_size\": %d,\n", config->flank);
    json_uint(file, "region_length", result->flank_bases, 1);
    json_double(file, "average_depth", flank_average, 1);
    json_double(file, "average_depth_with_deletions",
                result->flank_bases ? (double)result->flank_coverage_data /
                                      (double)result->flank_bases : 0.0, 1);
    json_uint(file, "flank_reads", stats->n_flk, 1);
    json_pct(file, "flank_reads_fraction_all", stats->n_flk, stats->n_reads, 1);
    json_pct(file, "flank_reads_fraction_mapped", stats->n_flk, stats->n_mapped, 1);
    json_double(file, "flank_data_mb", (double)result->flank_data / 1e6, 1);
    json_double(file, "flank_data_with_deletions_mb",
                (double)result->flank_coverage_data / 1e6, 1);
    json_pct(file, "flank_data_fraction_all", result->flank_data, stats->n_data, 1);
    json_pct(file, "flank_data_fraction_mapped", result->flank_data, stats->n_mdata, 1);
    json_pct(file, "flank_data_with_deletions_fraction_all",
             result->flank_coverage_data, stats->n_data, 1);
    json_pct(file, "flank_data_with_deletions_fraction_mapped",
             result->flank_coverage_data, stats->n_mdata, 1);
    fputs("    \"coverage\": {\n", file);
    json_coverage(file, "      ", result->flank_depth, result->flank_bases, config, flank_average, 0);
    fputs("    },\n    \"coverage_with_deletions\": {\n", file);
    json_coverage(file, "      ", result->flank_coverage_depth, result->flank_bases, config,
                  result->flank_bases ? (double)result->flank_coverage_data /
                  (double)result->flank_bases : 0.0, 0);
    fputs("    }\n  }\n}\n", file);
    return close_text_file(file, path);
}

static int write_text_report(const char *path, const xamdst_config_t *config,
                             const analysis_result_t *result)
{
    FILE *file = fopen(path, "wb");
    if (file == NULL) {
        xerror("cannot create '%s': %s", path, strerror(errno));
        return -1;
    }
    (void)setvbuf(file, NULL, _IOFBF, REPORT_BUFFER_SIZE);
    const read_stats_t *s = &result->reads;
    uint64_t inserts = histogram_total(result->insert_sizes);
    double insert_average = inserts ? (double)histogram_sum(result->insert_sizes) / (double)inserts : 0.0;
    fprintf(file, "## The file was created by xamdst\n## Version : %s\n## Files :", XAMDST_VERSION);
    for (size_t i = 0; i < config->ninputs; ++i) fprintf(file, " %s", config->inputs[i]);
    fprintf(file, "\n[Total] Primary reads\t%" PRIu64 "\n[Total] QC Fail reads\t%" PRIu64
                 "\n[Total] Raw Data(Mb)\t%.6f\n[Total] Mapped Reads\t%" PRIu64
                 "\n[Total] Fraction of Mapped Reads\t%.6f%%\n[Total] Mapped Data(Mb)\t%.6f\n"
                 "[Total] Fraction of Mapped Data(Mb)\t%.6f%%\n[Total] Paired Reads\t%" PRIu64
                 "\n[Total] Properly paired\t%" PRIu64 "\n[Total] Read and mate paired\t%" PRIu64
                 "\n[Total] Singletons\t%" PRIu64 "\n[Total] Read and mate map to diff chr\t%" PRIu64
                 "\n[Total] Read1\t%" PRIu64 "\n[Total] Read2\t%" PRIu64
                 "\n[Total] Read1(rmdup)\t%" PRIu64 "\n[Total] Read2(rmdup)\t%" PRIu64
                 "\n[Total] forward strand reads\t%" PRIu64 "\n[Total] backward strand reads\t%" PRIu64
                 "\n[Total] PCR duplicate reads\t%" PRIu64 "\n[Total] Fraction of PCR duplicate reads\t%.6f%%\n"
                 "[Total] Map quality cutoff value\t%d\n[Total] MapQuality above cutoff reads\t%" PRIu64
                 "\n[Insert size] Average\t%.6f\n[Insert size] Median\t%.6f\n"
                 "[Target] Target Reads\t%" PRIu64 "\n[Target] Target Data(Mb)\t%.6f\n"
                 "[Target] Target Data Rmdup(Mb)\t%.6f\n[Target] Len of region\t%" PRIu64
                 "\n[Target] Average depth\t%.6f\n[Target] Average depth(rmdup)\t%.6f\n"
                 "[Target] Target Region Count\t%" PRIu64 "\n[flank] flank size\t%d\n[flank] Len of region (outside target)\t%" PRIu64
                 "\n[flank] Average depth\t%.6f\n[flank] flank Reads\t%" PRIu64 "\n",
            s->n_reads, s->n_qcfail, (double)s->n_data / 1e6, s->n_mapped,
            percentage(s->n_mapped, s->n_reads), (double)s->n_mdata / 1e6,
            percentage(s->n_mdata, s->n_data), s->n_pair_all, s->n_pair_good, s->n_pair_map,
            s->n_sgltn, s->n_diffchr, s->n_read1, s->n_read2, s->n_rmdup1, s->n_rmdup2,
            s->n_pstrand, s->n_mstrand, s->n_dup, percentage(s->n_dup, s->n_mapped), config->mapq,
            s->n_qual, insert_average, histogram_median(result->insert_sizes), s->n_tgt,
            (double)result->target_data / 1e6, (double)result->target_rmdup_data / 1e6,
            result->target_bases, result->target_bases ? (double)result->target_data / result->target_bases : 0.0,
             result->target_bases ? (double)result->target_rmdup_data / result->target_bases : 0.0,
             result->target_region_count, config->flank, result->flank_bases,
             result->flank_bases ? (double)result->flank_data / result->flank_bases : 0.0, s->n_flk);
    fprintf(file, "[Target] Target Data with deletions(Mb)\t%.6f\n"
                 "[Target] Average depth with deletions\t%.6f\n"
                 "[flank] Average depth with deletions\t%.6f\n"
                 "[flank] flank Data with deletions(Mb)\t%.6f\n",
            (double)result->target_coverage_data / 1e6,
            result->target_bases ? (double)result->target_coverage_data / result->target_bases : 0.0,
            result->flank_bases ? (double)result->flank_coverage_data / result->flank_bases : 0.0,
            (double)result->flank_coverage_data / 1e6);
    fprintf(file, "[Fragment] Enabled\t%s\n[Fragment] Matched pairs\t%" PRIu64
                 "\n[Fragment] Overlap intervals removed\t%" PRIu64
                 "\n[Fragment] Unmatched pairs\t%" PRIu64
                 "\n[Fragment] Ambiguous pairs\t%" PRIu64 "\n",
            config->fragment_mode ? "true" : "false", result->fragments.paired,
            result->fragments.clipped, result->fragments.unmatched,
            result->fragments.ambiguous);
    fprintf(file, "[Total] Fraction of Properly paired primary reads\t%.6f%%\n"
                 "[Total] Fraction of Read and mate paired primary reads\t%.6f%%\n"
                 "[Total] Fraction of MapQ reads in all primary reads\t%.6f%%\n"
                 "[Total] Fraction of MapQ reads in mapped reads\t%.6f%%\n"
                 "[Target] Fraction of Target Reads in all primary reads\t%.6f%%\n"
                 "[Target] Fraction of Target Reads in mapped reads\t%.6f%%\n"
                 "[Target] Fraction of Target Data in all data\t%.6f%%\n"
                 "[Target] Fraction of Target Data in mapped data\t%.6f%%\n",
             percentage(s->n_pair_good, s->n_reads),
             percentage(s->n_pair_map, s->n_reads),
             percentage(s->n_qual, s->n_reads),
             percentage(s->n_qual, s->n_mapped),
             percentage(s->n_tgt, s->n_reads),
             percentage(s->n_tgt, s->n_mapped),
             percentage(result->target_data, s->n_data),
             percentage(result->target_data, s->n_mdata));
    fprintf(file, "[Target] Coverage >0x\t%.6f%%\n[Target] Coverage >=4x\t%.6f%%\n"
                 "[Target] Coverage >=10x\t%.6f%%\n[Target] Coverage >=30x\t%.6f%%\n"
                 "[Target] Coverage >=100x\t%.6f%%\n[flank] Coverage >0x\t%.6f%%\n"
                 "[flank] Coverage >=4x\t%.6f%%\n[flank] Coverage >=10x\t%.6f%%\n"
                 "[flank] Coverage >=30x\t%.6f%%\n[flank] Coverage >=100x\t%.6f%%\n",
            percentage(histogram_count_at_least(result->target_depth, 1), result->target_bases),
            percentage(histogram_count_at_least(result->target_depth, 4), result->target_bases),
            percentage(histogram_count_at_least(result->target_depth, 10), result->target_bases),
            percentage(histogram_count_at_least(result->target_depth, 30), result->target_bases),
            percentage(histogram_count_at_least(result->target_depth, 100), result->target_bases),
            percentage(histogram_count_at_least(result->flank_depth, 1), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_depth, 4), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_depth, 10), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_depth, 30), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_depth, 100), result->flank_bases));
    fprintf(file, "[Target] Coverage with deletions >0x\t%.6f%%\n"
                 "[Target] Coverage with deletions >=4x\t%.6f%%\n"
                 "[Target] Coverage with deletions >=10x\t%.6f%%\n"
                 "[Target] Coverage with deletions >=30x\t%.6f%%\n"
                 "[Target] Coverage with deletions >=100x\t%.6f%%\n"
                 "[flank] Coverage with deletions >0x\t%.6f%%\n"
                 "[flank] Coverage with deletions >=4x\t%.6f%%\n"
                 "[flank] Coverage with deletions >=10x\t%.6f%%\n"
                 "[flank] Coverage with deletions >=30x\t%.6f%%\n"
                 "[flank] Coverage with deletions >=100x\t%.6f%%\n",
            percentage(histogram_count_at_least(result->target_coverage_depth, 1), result->target_bases),
            percentage(histogram_count_at_least(result->target_coverage_depth, 4), result->target_bases),
            percentage(histogram_count_at_least(result->target_coverage_depth, 10), result->target_bases),
            percentage(histogram_count_at_least(result->target_coverage_depth, 30), result->target_bases),
            percentage(histogram_count_at_least(result->target_coverage_depth, 100), result->target_bases),
            percentage(histogram_count_at_least(result->flank_coverage_depth, 1), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_coverage_depth, 4), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_coverage_depth, 10), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_coverage_depth, 30), result->flank_bases),
            percentage(histogram_count_at_least(result->flank_coverage_depth, 100), result->flank_bases));
    fprintf(file, "[Target] Coverage(rmdup) >0x\t%.6f%%\n"
                 "[Target] Coverage(rmdup) >=4x\t%.6f%%\n"
                 "[Target] Coverage(rmdup) >=10x\t%.6f%%\n"
                 "[Target] Coverage(rmdup) >=30x\t%.6f%%\n"
                 "[Target] Coverage(rmdup) >=100x\t%.6f%%\n"
                 "[Target] Region covered >0x\t%" PRIu64 "\n"
                 "[Target] Fraction Region covered >0x\t%.6f%%\n"
                 "[Target] Fraction Region covered >=4x\t%.6f%%\n"
                 "[Target] Fraction Region covered >=10x\t%.6f%%\n"
                 "[Target] Fraction Region covered >=30x\t%.6f%%\n"
                 "[Target] Fraction Region covered >=100x\t%.6f%%\n"
                 "[flank] Fraction of flank Reads in all primary reads\t%.6f%%\n"
                 "[flank] Fraction of flank Reads in mapped reads\t%.6f%%\n"
                 "[flank] flank Data(Mb)\t%.6f\n"
                 "[flank] Fraction of flank Data in all data\t%.6f%%\n"
                 "[flank] Fraction of flank Data in mapped data\t%.6f%%\n",
             percentage(histogram_count_at_least(result->target_rmdup_depth, 1), result->target_bases),
             percentage(histogram_count_at_least(result->target_rmdup_depth, 4), result->target_bases),
             percentage(histogram_count_at_least(result->target_rmdup_depth, 10), result->target_bases),
             percentage(histogram_count_at_least(result->target_rmdup_depth, 30), result->target_bases),
             percentage(histogram_count_at_least(result->target_rmdup_depth, 100), result->target_bases),
             histogram_count_at_least(result->region_means, 1),
             percentage(histogram_count_at_least(result->region_means, 1),
                        histogram_total(result->region_means)),
             percentage(histogram_count_at_least(result->region_means, 4),
                        histogram_total(result->region_means)),
             percentage(histogram_count_at_least(result->region_means, 10),
                        histogram_total(result->region_means)),
             percentage(histogram_count_at_least(result->region_means, 30),
                        histogram_total(result->region_means)),
             percentage(histogram_count_at_least(result->region_means, 100),
                        histogram_total(result->region_means)),
             percentage(s->n_flk, s->n_reads),
             percentage(s->n_flk, s->n_mapped),
             (double)result->flank_data / 1e6,
             percentage(result->flank_data, s->n_data),
             percentage(result->flank_data, s->n_mdata));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[Target] Coverage >=%dx\t%.6f%%\n", config->cutoffs[i],
                percentage(histogram_count_at_least(result->target_depth, (uint64_t)config->cutoffs[i]),
                           result->target_bases));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[Target] Coverage(rmdup) >=%dx\t%.6f%%\n", config->cutoffs[i],
                percentage(histogram_count_at_least(result->target_rmdup_depth,
                                                     (uint64_t)config->cutoffs[i]),
                           result->target_bases));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[Target] Coverage with deletions >=%dx\t%.6f%%\n", config->cutoffs[i],
                percentage(histogram_count_at_least(result->target_coverage_depth,
                                                     (uint64_t)config->cutoffs[i]),
                           result->target_bases));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[Target] Fraction Region covered >=%dx\t%.6f%%\n",
                config->cutoffs[i],
                percentage(histogram_count_at_least(result->region_means,
                                                     (uint64_t)config->cutoffs[i]),
                           histogram_total(result->region_means)));
    for (size_t i = 0; i < config->nratios; ++i)
        fprintf(file, "[Target] Coverage >%.6f*Avg\t%.6f%%\n", config->ratios[i],
                percentage(count_greater_than(result->target_depth, config->ratios[i] *
                                              (result->target_bases ? (double)result->target_data / result->target_bases : 0.0)),
                           result->target_bases));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[flank] Coverage >=%dx\t%.6f%%\n", config->cutoffs[i],
                percentage(histogram_count_at_least(result->flank_depth,
                                                     (uint64_t)config->cutoffs[i]),
                           result->flank_bases));
    for (size_t i = 0; i < config->ncutoffs; ++i)
        fprintf(file, "[flank] Coverage with deletions >=%dx\t%.6f%%\n",
                config->cutoffs[i],
                percentage(histogram_count_at_least(result->flank_coverage_depth,
                                                     (uint64_t)config->cutoffs[i]),
                           result->flank_bases));
    return close_text_file(file, path);
}

int report_finish(report_writer_t *writer, const xamdst_config_t *config,
                  const interval_set_t *intervals, const analysis_result_t *result)
{
    (void)intervals;
    int status = 0;
    if (writer->depth != NULL && flush_buffer(writer->depth, &writer->depth_buffer) != 0)
        status = -1;
    if (writer->region != NULL && flush_buffer(writer->region, &writer->region_buffer) != 0)
        status = -1;
    if (writer->depth != NULL && !writer->depth_closed) {
        int close_status = bgzf_close(writer->depth);
        writer->depth = NULL;
        writer->depth_closed = 1;
        if (close_status != 0) status = -1;
    }
    if (writer->region != NULL && !writer->region_closed) {
        int close_status = bgzf_close(writer->region);
        writer->region = NULL;
        writer->region_closed = 1;
        if (close_status != 0) status = -1;
    }
    if (writer->uncovered != NULL && !writer->uncovered_closed) {
        int write_error = ferror(writer->uncovered);
        int close_error = fclose(writer->uncovered);
        if (write_error || close_error != 0) {
            xerror("failed to close uncover.bed: %s",
                   close_error != 0 ? strerror(errno) : "stream error");
            writer->uncovered = NULL;
            writer->uncovered_closed = 1;
            status = -1;
        }
        writer->uncovered = NULL;
        writer->uncovered_closed = 1;
    }
    if (status != 0)
        return -1;
    if (write_text_report(writer->temporary_paths[0], config, result) ||
        write_json(writer->temporary_paths[1], config, result) ||
        write_plot(writer->temporary_paths[2], result->target_depth, config->maxdepth) ||
        write_plot(writer->temporary_paths[3], result->insert_sizes, 0) ||
        write_chromosome_report(writer->temporary_paths[4], config, result))
        return -1;
    return 0;
}

int report_commit(report_writer_t *writer)
{
    int moved_old[OUTPUT_COUNT] = {0};
    int installed[OUTPUT_COUNT] = {0};

    /* Reserve backup names before moving anything.  This avoids overwriting a
     * stale backup left by an interrupted process (or a user file). */
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        struct stat existing;
        if (stat(writer->final_paths[i], &existing) == 0 && !S_ISREG(existing.st_mode)) {
            xerror("output path is not a regular file '%s'", writer->final_paths[i]);
            report_abort(writer);
            return -1;
        }
        writer->backup_paths[i] = temporary_path(writer->final_paths[i], 10000U + (unsigned)i);
        if (access(writer->backup_paths[i], F_OK) == 0) {
            xerror("stale output backup exists '%s'", writer->backup_paths[i]);
            for (size_t j = 0; j < OUTPUT_COUNT; ++j) {
                free(writer->backup_paths[j]);
                writer->backup_paths[j] = NULL;
            }
            report_abort(writer);
            return -1;
        }
    }

    /* Move existing results out of the way.  If any subsequent installation
     * fails, the rollback below restores the complete previous set instead of
     * leaving a mixture of old and new reports. */
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        if (rename(writer->final_paths[i], writer->backup_paths[i]) == 0) {
            moved_old[i] = 1;
            writer->backup_created[i] = 1;
        } else if (errno != ENOENT) {
            xerror("cannot stage existing output '%s': %s", writer->final_paths[i],
                   strerror(errno));
            for (size_t j = 0; j < OUTPUT_COUNT; ++j) {
                if (moved_old[j]) {
                    if (rename(writer->backup_paths[j], writer->final_paths[j]) != 0) {
                        xerror("cannot restore previous output '%s': %s", writer->final_paths[j],
                               strerror(errno));
                    } else {
                        writer->backup_created[j] = 0;
                    }
                }
            }
            report_abort(writer);
            return -1;
        }
    }
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        if (writer->temporary_paths[i] == NULL || !writer->temporary_created[i])
            continue;
        if (rename(writer->temporary_paths[i], writer->final_paths[i]) != 0) {
            xerror("cannot install output '%s': %s", writer->final_paths[i], strerror(errno));
            for (size_t j = 0; j < OUTPUT_COUNT; ++j) {
                if (installed[j])
                    (void)remove(writer->final_paths[j]);
                else if (writer->temporary_paths[j] != NULL && writer->temporary_created[j])
                    (void)remove(writer->temporary_paths[j]);
            }
            for (size_t j = 0; j < OUTPUT_COUNT; ++j) {
                if (moved_old[j]) {
                    if (rename(writer->backup_paths[j], writer->final_paths[j]) != 0) {
                        xerror("cannot restore previous output '%s': %s", writer->final_paths[j],
                               strerror(errno));
                    } else {
                        writer->backup_created[j] = 0;
                    }
                }
            }
            report_abort(writer);
            return -1;
        }
        installed[i] = 1;
        writer->temporary_created[i] = 0;
    }
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        if (moved_old[i] && remove(writer->backup_paths[i]) != 0)
            xwarn("cannot remove output backup '%s': %s", writer->backup_paths[i],
                  strerror(errno));
        writer->backup_created[i] = 0;
    }
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        free(writer->final_paths[i]);
        free(writer->temporary_paths[i]);
        free(writer->backup_paths[i]);
        writer->final_paths[i] = writer->temporary_paths[i] = NULL;
        writer->backup_paths[i] = NULL;
    }
    free(writer->depth_buffer.s);
    free(writer->region_buffer.s);
    writer->depth_buffer = (kstring_t){0, 0, NULL};
    writer->region_buffer = (kstring_t){0, 0, NULL};
    writer->open = 0;
    return 0;
}

void report_abort(report_writer_t *writer)
{
    if (writer == NULL)
        return;
    if (writer->depth != NULL) bgzf_close(writer->depth);
    if (writer->region != NULL) bgzf_close(writer->region);
    if (writer->uncovered != NULL) fclose(writer->uncovered);
    writer->depth = NULL;
    writer->region = NULL;
    writer->uncovered = NULL;
    free(writer->depth_buffer.s);
    free(writer->region_buffer.s);
    writer->depth_buffer = (kstring_t){0, 0, NULL};
    writer->region_buffer = (kstring_t){0, 0, NULL};
    for (size_t i = 0; i < OUTPUT_COUNT; ++i) {
        if (writer->temporary_paths[i] != NULL) {
            if (writer->temporary_created[i])
                (void)remove(writer->temporary_paths[i]);
            free(writer->temporary_paths[i]);
        }
        if (writer->backup_paths[i] != NULL) {
            if (writer->backup_created[i]) {
                xwarn("preserving output backup '%s' after rollback failure",
                      writer->backup_paths[i]);
            }
            free(writer->backup_paths[i]);
        }
        free(writer->final_paths[i]);
        writer->temporary_paths[i] = writer->final_paths[i] = NULL;
        writer->backup_paths[i] = NULL;
        writer->temporary_created[i] = 0;
        writer->backup_created[i] = 0;
    }
    writer->open = 0;
}
