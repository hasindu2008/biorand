#ifndef F5CMISC_H
#define F5CMISC_H

#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>

// taken from minimap2/misc
static inline double realtime(void) {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

// taken from minimap2/misc
static inline double cputime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

// Prints to the provided buffer a nice number of bytes (KB, MB, GB, etc)
//from https://www.mbeckler.org/blog/?p=114
static inline void print_size(const char* name, uint64_t bytes)
{
    const char* suffixes[7];
    suffixes[0] = "B";
    suffixes[1] = "KB";
    suffixes[2] = "MB";
    suffixes[3] = "GB";
    suffixes[4] = "TB";
    suffixes[5] = "PB";
    suffixes[6] = "EB";
    uint64_t s = 0; // which suffix to use
    double count = bytes;
    while (count >= 1024 && s < 7)
    {
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
        fprintf(stderr, "%s : %d %s\n", name, (int)count, suffixes[s]);
    else
        fprintf(stderr, "%s : %.1f %s\n", name, count, suffixes[s]);
}


#define WARN "[%s::WARNING]\033[1;33m "
#define ERR "[%s::ERROR]\033[1;31m "
#define CEND "\033[0m\n"

#define WARNING(arg, ...)                                                      \
    fprintf(stderr, "[%s::WARNING]\033[1;33m " arg "\033[0m\n", __func__,      \
            __VA_ARGS__)
#define ERROR(arg, ...)                                                        \
    fprintf(stderr, "[%s::ERROR]\033[1;31m " arg "\033[0m\n", __func__,        \
            __VA_ARGS__)
#define INFO(arg, ...)                                                         \
    fprintf(stderr, "[%s::INFO]\033[1;34m " arg "\033[0m\n", __func__,         \
            __VA_ARGS__)
#define SUCCESS(arg, ...)                                                      \
    fprintf(stderr, "[%s::SUCCESS]\033[1;32m " arg "\033[0m\n", __func__,      \
            __VA_ARGS__)
#define DEBUG(arg, ...)                                                        \
    fprintf(stderr,                                                            \
            "[%s::DEBUG]\033[1;35m Error occured at %s:%d. " arg "\033[0m\n",  \
            __func__, __FILE__, __LINE__ - 2, __VA_ARGS__)

#define MALLOC_CHK(ret) malloc_chk((void*)ret, __func__, __FILE__, __LINE__ - 1)
#define F_CHK(ret, filename)                                                   \
    f_chk((void*)ret, __func__, __FILE__, __LINE__ - 1, filename);
#define NULL_CHK(ret) null_chk((void*)ret, __func__, __FILE__, __LINE__ - 1)
#define NEG_CHK(ret) neg_chk(ret, __func__, __FILE__, __LINE__ - 1)

static inline void malloc_chk(void* ret, const char* func, const char* file,
                              int line) {
    if (ret != NULL)
        return;
    fprintf(
        stderr,
        "[%s::ERROR]\033[1;31m Failed to allocate memory : "
        "%s.\033[0m\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n",
        func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

static inline void f_chk(void* ret, const char* func, const char* file,
                         int line, const char* fopen_f) {
    if (ret != NULL)
        return;
    fprintf(
        stderr,
        "[%s::ERROR]\033[1;31m Failed to open %s : "
        "%s.\033[0m\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n",
        func, fopen_f, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

// Die on error. Print the error and exit if the return value of the previous function NULL
static inline void null_chk(void* ret, const char* func, const char* file,
                            int line) {
    if (ret != NULL)
        return;
    fprintf(stderr,
            "[%s::ERROR]\033[1;31m %s.\033[0m\n[%s::DEBUG]\033[1;35m Error "
            "occured at %s:%d.\033[0m\n\n",
            func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

// Die on error. Print the error and exit if the return value of the previous function is -1
static inline void neg_chk(int ret, const char* func, const char* file,
                           int line) {
    if (ret >= 0)
        return;
    fprintf(stderr,
            "[%s::ERROR]\033[1;31m %s.\033[0m\n[%s::DEBUG]\033[1;35m Error "
            "occured at %s:%d.\033[0m\n\n",
            func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

#endif
