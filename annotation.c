#include "annotation.h"
#include "util.h"
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

typedef struct { uint64_t start, end; } range_t;
typedef struct { range_t *items; size_t count, capacity; } ranges_t;
typedef struct { ranges_t exon, intron, plus, minus, body; } features_t;
typedef struct { int tid; range_t span; } gene_t;
KHASH_MAP_INIT_STR(annotation_genes, gene_t)
struct annotation {
    features_t *chromosomes;
    size_t count;
    uint64_t exon_bases, intron_bases, span_bases;
};

static void append(ranges_t *v, uint64_t start, uint64_t end)
{
    if (start >= end) return;
    if (v->count == v->capacity) {
        size_t next = v->capacity ? v->capacity * 2 : 16;
        if (next < v->capacity) { xerror("annotation too large"); exit(EXIT_FAILURE); }
        v->items = xreallocarray(v->items, next, sizeof(*v->items));
        v->capacity = next;
    }
    v->items[v->count++] = (range_t){start, end};
}

static int compare_range(const void *left, const void *right)
{
    const range_t *a = left, *b = right;
    if (a->start != b->start) return a->start < b->start ? -1 : 1;
    return (a->end > b->end) - (a->end < b->end);
}

static void merge_ranges(ranges_t *v)
{
    if (v->count == 0) return;
    qsort(v->items, v->count, sizeof(*v->items), compare_range);
    size_t n = 1;
    for (size_t i = 1; i < v->count; ++i) {
        range_t *last = &v->items[n-1];
        if (v->items[i].start <= last->end) {
            if (v->items[i].end > last->end) last->end = v->items[i].end;
        } else v->items[n++] = v->items[i];
    }
    v->count = n;
}

static uint64_t overlap(const ranges_t *v, uint64_t start, uint64_t end)
{
    size_t lo = 0, hi = v->count;
    while (lo < hi) {
        size_t mid = lo + (hi-lo)/2;
        if (v->items[mid].end <= start) lo = mid + 1; else hi = mid;
    }
    uint64_t bases = 0;
    for (size_t i = lo; i < v->count && v->items[i].start < end; ++i) {
        uint64_t a = start > v->items[i].start ? start : v->items[i].start;
        uint64_t b = end < v->items[i].end ? end : v->items[i].end;
        bases += b-a; /* disjoint intervals, bounded by end-start */
    }
    return bases;
}

/* Parse GTF and GFF2 tag/value attributes (quoted strings may contain ';').
 * gene_id/gene take precedence over transcript_id/Transcript. */
static int group_id(const char *text, char **out)
{
    if (text[0] == '\0' || !strcmp(text, ".")) { *out = NULL; return 0; }
    const char *p = text;
    char *gene = NULL, *transcript = NULL;
    while (*p) {
        while (*p && (isspace((unsigned char)*p) || *p == ';')) ++p;
        if (!*p || *p == '#') break;
        const char *key = p;
        while (*p && !isspace((unsigned char)*p) && *p != ';' && *p != '=') ++p;
        size_t nk = (size_t)(p-key);
        if (*p == '=' || *p == ';' || !*p) goto fail;
        while (isspace((unsigned char)*p)) ++p;
        char *value = xmalloc(strlen(p)+1);
        size_t n = 0;
        if (*p == '"') {
            ++p;
            while (*p && *p != '"') {
                if (*p == '\\' && p[1]) ++p;
                value[n++] = *p++;
            }
            if (*p != '"') { free(value); goto fail; }
            ++p;
        } else {
            while (*p && *p != ';' && !isspace((unsigned char)*p)) value[n++] = *p++;
        }
        value[n] = '\0';
        while (isspace((unsigned char)*p)) ++p;
        /* GFF2 permits additional unquoted values in attributes we ignore. */
        while (*p && *p != ';') ++p;
        int is_gene = (nk == 7 && !strncmp(key,"gene_id",7)) ||
                      (nk == 4 && !strncmp(key,"gene",4));
        int is_transcript = (nk == 13 && !strncmp(key,"transcript_id",13)) ||
                            (nk == 10 && !strncmp(key,"Transcript",10));
        char **dest = is_gene ? &gene : is_transcript ? &transcript : NULL;
        if (dest && n) {
            if (*dest && strcmp(*dest,value)) { free(value); goto fail; }
            free(*dest); *dest = value;
        } else free(value);
    }
    *out = gene ? gene : transcript;
    if (gene) free(transcript);
    return 0;
fail:
    free(gene); free(transcript);
    return -1;
}

static int coordinate(const char *s, uint64_t *out)
{
    if (!isdigit((unsigned char)*s)) return -1;
    char *end;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno || *end || v == 0 || v > INT64_MAX) return -1;
    *out = (uint64_t)v;
    return 0;
}

int annotation_load(annotation_t **out, const char *path, const sam_hdr_t *header)
{
    *out = NULL;
    htsFile *file = hts_open(path, "r");
    if (!file) { xerror("cannot open annotation '%s'", path); return -1; }
    annotation_t *a = xcalloc(1, sizeof(*a));
    a->count = (size_t)sam_hdr_nref(header);
    a->chromosomes = xcalloc(a->count ? a->count : 1, sizeof(*a->chromosomes));
    khash_t(annotation_genes) *genes = kh_init(annotation_genes);
    kstring_t line = {0,0,NULL};
    uint64_t line_number = 0, matched = 0, skipped = 0;
    int status = -1, read_status = 0;
    if (!genes) goto done;
    for (size_t i = 0; i < a->count; ++i) {
        hts_pos_t length = sam_hdr_tid2len(header, (int)i);
        if (length < 0 || u64_add(a->span_bases, (uint64_t)length, &a->span_bases))
            goto done;
    }
    while ((read_status = hts_getline(file, '\n', &line)) >= 0) {
        ++line_number;
        if (line.l && line.s[line.l-1] == '\r') line.s[--line.l] = 0;
        if (!line.l || line.s[0] == '#') continue;
        char *fields[9]; fields[0] = line.s;
        char *p = line.s;
        for (size_t i = 1; i < 9; ++i) {
            p = strchr(p, '\t');
            if (!p) goto malformed;
            *p++ = 0; fields[i] = p;
        }
        if (strchr(fields[8], '\t')) goto malformed;
        uint64_t start, end;
        if (coordinate(fields[3], &start) || coordinate(fields[4], &end) ||
            start > end || strlen(fields[6]) != 1 || !strchr("+-.?",fields[6][0]))
            goto malformed;
        int exon = !strcmp(fields[2], "exon");
        int body = !strcmp(fields[2], "gene") || !strcmp(fields[2], "transcript");
        if (!exon && !body) continue;
        int tid = sam_hdr_name2tid((sam_hdr_t *)header, fields[0]);
        if (tid < 0) { ++skipped; continue; }
        if (end > (uint64_t)sam_hdr_tid2len(header,tid)) goto malformed;
        --start;
        char *id = NULL;
        if (group_id(fields[8], &id)) goto malformed;
        features_t *f = &a->chromosomes[tid];
        if (id) {
            size_t size = strlen(id) + 32;
            char *key = xmalloc(size);
            snprintf(key,size,"%d:%s",tid,id);
            free(id);
            int ret;
            khiter_t k = kh_put(annotation_genes,genes,key,&ret);
            if (ret < 0) { free(key); goto done; }
            if (ret) kh_val(genes,k) = (gene_t){tid,{start,end}};
            else {
                free(key);
                range_t *r = &kh_val(genes,k).span;
                if (start < r->start) r->start = start;
                if (end > r->end) r->end = end;
            }
        }
        if (body) append(&f->body,start,end);
        if (exon) {
            append(&f->exon,start,end);
            if (fields[6][0] == '+') append(&f->plus,start,end);
            if (fields[6][0] == '-') append(&f->minus,start,end);
        }
        ++matched;
        continue;
malformed:
        xerror("invalid GTF/GFF2 annotation at %s:%llu",path,(unsigned long long)line_number);
        goto done;
    }
    if (read_status < -1 || !matched) {
        xerror("annotation read failed or no exon/gene/transcript features match the BAM header");
        goto done;
    }
    for (khiter_t k = kh_begin(genes); k != kh_end(genes); ++k)
        if (kh_exist(genes,k)) {
            gene_t g = kh_val(genes,k);
            append(&a->chromosomes[g.tid].body,g.span.start,g.span.end);
        }
    for (size_t c = 0; c < a->count; ++c) {
        features_t *f = &a->chromosomes[c];
        merge_ranges(&f->exon); merge_ranges(&f->body);
        merge_ranges(&f->plus); merge_ranges(&f->minus);
        size_t e = 0;
        for (size_t b = 0; b < f->body.count; ++b) {
            uint64_t pos = f->body.items[b].start, end = f->body.items[b].end;
            while (e < f->exon.count && f->exon.items[e].end <= pos) ++e;
            for (size_t j = e; j < f->exon.count && f->exon.items[j].start < end; ++j) {
                append(&f->intron,pos,f->exon.items[j].start < end ? f->exon.items[j].start : end);
                if (f->exon.items[j].end > pos) pos = f->exon.items[j].end;
            }
            append(&f->intron,pos,end);
        }
        for (size_t j = 0; j < f->exon.count; ++j)
            if (u64_add(a->exon_bases,f->exon.items[j].end-f->exon.items[j].start,&a->exon_bases)) goto done;
        for (size_t j = 0; j < f->intron.count; ++j)
            if (u64_add(a->intron_bases,f->intron.items[j].end-f->intron.items[j].start,&a->intron_bases)) goto done;
        free(f->body.items); memset(&f->body,0,sizeof(f->body));
    }
    if (skipped) xwarn("ignored %llu annotation features on references absent from the BAM header",(unsigned long long)skipped);
    status = 0;
done:
    if (genes) {
        for (khiter_t k = kh_begin(genes); k != kh_end(genes); ++k)
            if (kh_exist(genes,k)) free((char *)kh_key(genes,k));
        kh_destroy(annotation_genes,genes);
    }
    free(line.s);
    if (hts_close(file) != 0) status = -1;
    if (status != 0) annotation_destroy(a); else *out = a;
    return status;
}

void annotation_destroy(annotation_t *a)
{
    if (!a) return;
    for (size_t i = 0; i < a->count; ++i) {
        features_t *f = &a->chromosomes[i];
        free(f->exon.items); free(f->intron.items); free(f->plus.items);
        free(f->minus.items); free(f->body.items);
    }
    free(a->chromosomes); free(a);
}

int annotation_overlap(const annotation_t *a, int32_t tid, uint64_t start,
                       uint64_t end, annotation_overlap_t *result)
{
    memset(result,0,sizeof(*result));
    if (!a || tid < 0 || (size_t)tid >= a->count || start > end) return -1;
    const features_t *f = &a->chromosomes[tid];
    result->exon = overlap(&f->exon,start,end);
    result->intron = overlap(&f->intron,start,end);
    result->exon_plus = overlap(&f->plus,start,end);
    result->exon_minus = overlap(&f->minus,start,end);
    return 0;
}

uint64_t annotation_exon_bases(const annotation_t *a) { return a ? a->exon_bases : 0; }
uint64_t annotation_intron_bases(const annotation_t *a) { return a ? a->intron_bases : 0; }
uint64_t annotation_span_bases(const annotation_t *a) { return a ? a->span_bases : 0; }
