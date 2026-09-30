#ifndef XAMDST_RNA_DEDUP_H
#define XAMDST_RNA_DEDUP_H
#include <stddef.h>
#include <stdint.h>
typedef struct rna_dedup rna_dedup_t;
rna_dedup_t *rna_dedup_create(const char *outdir, size_t memory_limit);
int rna_dedup_add(rna_dedup_t *state, size_t input, unsigned end, const char *name);
int rna_dedup_finish(rna_dedup_t *state, uint64_t *count);
void rna_dedup_destroy(rna_dedup_t *state);
#endif
