#include "rna_dedup.h"
#include "util.h"
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Bounded sorted partitions, compacted like a binary counter. At most 64
 * partitions survive between flushes; merging uses two record buffers and
 * never relies on a hash collision or an external sorting executable. */
struct rna_dedup {
    const char *outdir;
    size_t limit, bytes, count, capacity;
    char **keys;
    char *runs[64];
};

rna_dedup_t *rna_dedup_create(const char *outdir, size_t memory_limit)
{
    rna_dedup_t *s = calloc(1, sizeof(*s));
    if (s == NULL) { xerror("cannot allocate RNA dedup state"); return NULL; }
    s->outdir = outdir;
    s->limit = memory_limit;
    return s;
}

static FILE *scratch(rna_dedup_t *s, char **path)
{
    size_t n = strlen(s->outdir);
    const char suffix[] = "/.xamdst-rna-dedup-XXXXXX";
    if (n > SIZE_MAX-sizeof(suffix)) return NULL;
    *path = malloc(n+sizeof(suffix));
    if (*path == NULL) { xerror("cannot allocate RNA dedup path"); return NULL; }
    memcpy(*path,s->outdir,n);
    memcpy(*path+n,suffix,sizeof(suffix));
    int fd = mkstemp(*path);
    FILE *f = fd >= 0 ? fdopen(fd, "w+b") : NULL;
    if (f == NULL) {
        xerror("cannot create RNA dedup partition: %s", strerror(errno));
        if (fd >= 0) { close(fd); unlink(*path); }
        free(*path);
        *path = NULL;
    }
    return f;
}

static void discard(char **path)
{
    if (*path != NULL) {
        if (unlink(*path) != 0 && errno != ENOENT)
            xwarn("cannot remove RNA dedup partition '%s': %s", *path, strerror(errno));
        free(*path);
        *path = NULL;
    }
}

static int compare(const void *a, const void *b)
{
    return strcmp(*(char *const *)a, *(char *const *)b);
}

static int merge(rna_dedup_t *s, const char *left, const char *right, char **output)
{
    FILE *a = fopen(left, "rb"), *b = fopen(right, "rb"), *o = NULL;
    char *la = NULL, *lb = NULL;
    size_t ca = 0, cb = 0;
    int status = -1;
    if (a == NULL || b == NULL || (o = scratch(s, output)) == NULL)
        goto done;
    ssize_t na = getline(&la, &ca, a), nb = getline(&lb, &cb, b);
    while (na >= 0 || nb >= 0) {
        int cmp = na < 0 ? 1 : nb < 0 ? -1 : strcmp(la, lb);
        if (fputs(cmp <= 0 ? la : lb, o) == EOF)
            goto done;
        if (cmp <= 0) na = getline(&la, &ca, a);
        if (cmp >= 0) nb = getline(&lb, &cb, b);
    }
    if (feof(a) && feof(b) && !ferror(a) && !ferror(b)) status = 0;
done:
    free(la); free(lb);
    if (a != NULL && fclose(a) != 0) status = -1;
    if (b != NULL && fclose(b) != 0) status = -1;
    if (o != NULL && fclose(o) != 0) status = -1;
    if (status != 0) {
        xerror("failed to merge RNA dedup partitions");
        discard(output);
    }
    return status;
}

static int flush(rna_dedup_t *s)
{
    if (s->count == 0) return 0;
    char *path = NULL;
    FILE *f = scratch(s, &path);
    if (f == NULL) return -1;
    qsort(s->keys, s->count, sizeof(*s->keys), compare);
    int status = 0;
    for (size_t i = 0; i < s->count; ++i)
        if ((i == 0 || strcmp(s->keys[i], s->keys[i-1]) != 0) &&
            fprintf(f, "%s\n", s->keys[i]) < 0) { status = -1; break; }
    if (fclose(f) != 0) status = -1;
    if (status != 0) {
        xerror("failed to write RNA dedup partition");
        discard(&path);
        return -1;
    }
    for (size_t i = 0; i < s->count; ++i) free(s->keys[i]);
    free(s->keys);
    s->keys = NULL;
    s->capacity = s->count = s->bytes = 0;
    for (size_t level = 0; level < 64; ++level) {
        if (s->runs[level] == NULL) { s->runs[level] = path; return 0; }
        char *combined = NULL;
        status = merge(s, s->runs[level], path, &combined);
        discard(&path);
        if (status != 0) return -1;
        discard(&s->runs[level]);
        path = combined;
    }
    discard(&path);
    xerror("too many RNA dedup partitions");
    return -1;
}

int rna_dedup_add(rna_dedup_t *s, size_t input, unsigned end, const char *name)
{
    size_t n = strlen(name);
    if (n > (SIZE_MAX - 64) / 2) return -1;
    size_t length = n * 2 + 64;
    /* Include the next pointer-array growth before accepting a key. A single
     * oversized key is allowed so every legal BAM name remains supported. */
    size_t next_capacity = s->capacity;
    if (s->count == s->capacity) {
        next_capacity = s->capacity ? s->capacity * 2 : 4;
        if (next_capacity < s->capacity) return -1;
    }
    size_t pointer_bytes;
    if (size_mul(next_capacity,sizeof(*s->keys),&pointer_bytes)) return -1;
    int exceeds = length > s->limit || s->bytes > s->limit-length ||
                  pointer_bytes > s->limit-length-s->bytes;
    if (s->count && exceeds) {
        if (flush(s) != 0) return -1;
        next_capacity = 4;
        if (size_mul(next_capacity,sizeof(*s->keys),&pointer_bytes)) return -1;
    }
    if (s->count == s->capacity) {
        size_t bytes;
        if (size_mul(next_capacity,sizeof(*s->keys),&bytes)) return -1;
        char **keys = realloc(s->keys,bytes);
        if (keys == NULL) { xerror("cannot grow RNA dedup buffer"); return -1; }
        s->keys = keys;
        s->capacity = next_capacity;
    }
    char *key = malloc(length);
    if (key == NULL) { xerror("cannot allocate RNA dedup key"); return -1; }
    int prefix = snprintf(key, length, "%zu:%u:", input, end);
    if (prefix < 0 || (size_t)prefix >= 64) { free(key); return -1; }
    const char *hex = "0123456789abcdef";
    for (size_t i = 0; i < n; ++i) {
        unsigned c = (unsigned char)name[i];
        key[prefix + i * 2] = hex[c >> 4];
        key[prefix + i * 2 + 1] = hex[c & 15];
    }
    key[prefix + n * 2] = '\0';
    s->keys[s->count++] = key;
    s->bytes += length;
    return 0;
}

int rna_dedup_finish(rna_dedup_t *s, uint64_t *count)
{
    *count = 0;
    if (s == NULL) return 0;
    int spilled = 0;
    for (size_t i = 0; i < 64; ++i) spilled |= s->runs[i] != NULL;
    if (!spilled) {
        if (s->count == 0) return 0;
        qsort(s->keys, s->count, sizeof(*s->keys), compare);
        for (size_t i = 0; i < s->count; ++i)
            if (i == 0 || strcmp(s->keys[i], s->keys[i-1]) != 0) ++*count;
        return 0;
    }
    if (flush(s) != 0) return -1;
    size_t first = 0;
    while (first < 64 && s->runs[first] == NULL) ++first;
    if (first == 64) return 0;
    for (size_t i = first + 1; i < 64; ++i) {
        if (s->runs[i] == NULL) continue;
        char *combined = NULL;
        if (merge(s, s->runs[first], s->runs[i], &combined) != 0) return -1;
        discard(&s->runs[first]); discard(&s->runs[i]);
        s->runs[first] = combined;
    }
    FILE *f = fopen(s->runs[first], "rb");
    if (f == NULL) { xerror("cannot read RNA dedup partition"); return -1; }
    char *line = NULL;
    size_t capacity = 0;
    int status = 0;
    while (getline(&line, &capacity, f) >= 0) {
        if (*count == UINT64_MAX) { status = -1; break; }
        ++*count;
    }
    if (!feof(f) || ferror(f)) status = -1;
    free(line);
    if (fclose(f) != 0) status = -1;
    if (status != 0) xerror("failed to count RNA dedup partitions");
    return status;
}

void rna_dedup_destroy(rna_dedup_t *s)
{
    if (s == NULL) return;
    for (size_t i = 0; i < s->count; ++i) free(s->keys[i]);
    free(s->keys);
    for (size_t i = 0; i < 64; ++i) discard(&s->runs[i]);
    free(s);
}
