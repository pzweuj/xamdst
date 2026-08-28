#include "util.h"

#include <errno.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif

static void vmessage(const char *prefix, const char *fmt, va_list ap)
{
    fprintf(stderr, "%s: ", prefix);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

void xwarn(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vmessage("xamdst: warning", fmt, ap);
    va_end(ap);
}

void xerror(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vmessage("xamdst: error", fmt, ap);
    va_end(ap);
}

void *xmalloc(size_t size)
{
    if (size == 0)
        size = 1;
    void *ptr = malloc(size);
    if (ptr == NULL) {
        xerror("out of memory while allocating %zu bytes", size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void *xcalloc(size_t n, size_t size)
{
    size_t total;
    if (size_mul(n, size, &total)) {
        xerror("allocation size overflow (%zu x %zu)", n, size);
        exit(EXIT_FAILURE);
    }
    void *ptr = calloc(n, size);
    if (ptr == NULL && total != 0) {
        xerror("out of memory while allocating %zu bytes", total);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void *xreallocarray(void *ptr, size_t n, size_t size)
{
    size_t total;
    if (size_mul(n, size, &total)) {
        xerror("allocation size overflow (%zu x %zu)", n, size);
        exit(EXIT_FAILURE);
    }
    void *result = realloc(ptr, total == 0 ? 1 : total);
    if (result == NULL) {
        xerror("out of memory while allocating %zu bytes", total);
        exit(EXIT_FAILURE);
    }
    return result;
}

char *xstrdup(const char *value)
{
    if (value == NULL)
        return NULL;
    size_t length = strlen(value);
    if (length == SIZE_MAX) {
        xerror("string is too long to duplicate");
        exit(EXIT_FAILURE);
    }
    char *copy = xmalloc(length + 1);
    memcpy(copy, value, length + 1);
    return copy;
}

int size_mul(size_t a, size_t b, size_t *result)
{
    if (a != 0 && b > SIZE_MAX / a)
        return -1;
    *result = a * b;
    return 0;
}

int u64_add(uint64_t a, uint64_t b, uint64_t *result)
{
    if (b > UINT64_MAX - a)
        return -1;
    *result = a + b;
    return 0;
}

int mkdir_p(const char *path, unsigned mode)
{
    if (path == NULL || path[0] == '\0') {
        errno = EINVAL;
        return -1;
    }

    char *copy = xstrdup(path);
    size_t length = strlen(copy);
    while (length > 1 && (copy[length - 1] == '/' || copy[length - 1] == '\\'))
        copy[--length] = '\0';

    for (char *p = copy + 1; *p != '\0'; ++p) {
        if (*p != '/' && *p != '\\')
            continue;
        char saved = *p;
        *p = '\0';
        if (copy[0] != '\0' && mkdir(copy, (mode_t)mode) != 0 && errno != EEXIST) {
            int saved_errno = errno;
            free(copy);
            errno = saved_errno;
            return -1;
        }
        *p = saved;
    }
    if (mkdir(copy, (mode_t)mode) != 0 && errno != EEXIST) {
        int saved_errno = errno;
        free(copy);
        errno = saved_errno;
        return -1;
    }
    struct stat st;
    int result = stat(copy, &st);
    int saved_errno = errno;
    free(copy);
    if (result != 0)
        return -1;
    if (!S_ISDIR(st.st_mode)) {
        errno = ENOTDIR;
        return -1;
    }
    errno = saved_errno;
    return 0;
}

char *path_join(const char *directory, const char *basename)
{
    if (directory == NULL || basename == NULL)
        return NULL;
    size_t dl = strlen(directory);
    size_t bl = strlen(basename);
    int need_separator = dl > 0 && directory[dl - 1] != '/' && directory[dl - 1] != '\\';
    size_t total;
    if (dl > SIZE_MAX - (size_t)need_separator ||
        dl + (size_t)need_separator > SIZE_MAX - bl ||
        dl + (size_t)need_separator + bl == SIZE_MAX) {
        xerror("path is too long");
        exit(EXIT_FAILURE);
    }
    total = dl + (size_t)need_separator + bl + 1;
    char *result = xmalloc(total);
    memcpy(result, directory, dl);
    size_t pos = dl;
    if (need_separator)
        result[pos++] = '/';
    memcpy(result + pos, basename, bl);
    result[pos + bl] = '\0';
    return result;
}

char *temporary_path(const char *final_path, unsigned sequence)
{
    size_t length = strlen(final_path);
    char suffix[64];
    int written = snprintf(suffix, sizeof(suffix), ".tmp.%ld.%u", (long)getpid(), sequence);
    if (written < 0 || (size_t)written >= sizeof(suffix)) {
        xerror("failed to construct temporary output path");
        exit(EXIT_FAILURE);
    }
    if (length > SIZE_MAX - (size_t)written - 1) {
        xerror("temporary output path is too long");
        exit(EXIT_FAILURE);
    }
    char *result = xmalloc(length + (size_t)written + 1);
    memcpy(result, final_path, length);
    memcpy(result + length, suffix, (size_t)written + 1);
    return result;
}
