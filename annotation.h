#ifndef XAMDST_ANNOTATION_H
#define XAMDST_ANNOTATION_H

#include <stddef.h>
#include <stdint.h>

#include <htslib/sam.h>

typedef struct annotation annotation_t;

typedef struct {
    uint64_t exon;
    uint64_t intron;
    uint64_t exon_plus;
    uint64_t exon_minus;
} annotation_overlap_t;

int annotation_load(annotation_t **out, const char *path, const sam_hdr_t *header);
void annotation_destroy(annotation_t *annotation);
int annotation_overlap(const annotation_t *annotation, int32_t tid,
                      uint64_t start, uint64_t end,
                      annotation_overlap_t *overlap);
uint64_t annotation_exon_bases(const annotation_t *annotation);
uint64_t annotation_intron_bases(const annotation_t *annotation);
uint64_t annotation_span_bases(const annotation_t *annotation);

#endif
