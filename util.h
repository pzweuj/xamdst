#ifndef XAMDST_UTIL_H
#define XAMDST_UTIL_H

#include <stddef.h>
#include <stdint.h>

void xwarn(const char *fmt, ...);
void xerror(const char *fmt, ...);

void *xmalloc(size_t size);
void *xcalloc(size_t n, size_t size);
void *xreallocarray(void *ptr, size_t n, size_t size);
char *xstrdup(const char *value);

/* Create path recursively. Returns 0 on success and -1 on error. */
int mkdir_p(const char *path, unsigned mode);

/* Allocate a path consisting of directory + basename. */
char *path_join(const char *directory, const char *basename);

/* Return a unique temporary path beside final_path. */
char *temporary_path(const char *final_path, unsigned sequence);

/* Checked integer helpers used for coordinates and allocation sizes. */
int u64_add(uint64_t a, uint64_t b, uint64_t *result);
int size_mul(size_t a, size_t b, size_t *result);

#endif
