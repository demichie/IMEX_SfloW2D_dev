#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#else
#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __APPLE__
#include <sys/syscall.h>
#endif
#include <sys/uio.h>
#include <unistd.h>

#ifdef __APPLE__
#define DYLD_INTERPOSE(replacement, replacee) \
    __attribute__((used)) static struct { \
        const void *replacement; \
        const void *replacee; \
    } interpose_##replacee __attribute__((section("__DATA,__interpose"))) = { \
        (const void *)(unsigned long)&replacement, \
        (const void *)(unsigned long)&replacee \
    }
#endif

typedef ssize_t (*write_function)(int, const void *, size_t);
typedef int (*put_vara_function)(int, int, const size_t *,
                                 const size_t *, const void *);
typedef int (*inquire_varname_function)(int, int, char *);
typedef int (*inquire_vartype_function)(int, int, nc_type *);
typedef int (*sync_function)(int);

static int write_matches;
static int netcdf_time_records;
static int resolving_write;
#ifdef __APPLE__
extern ssize_t system_write_nocancel(int fd, const void *buffer, size_t count)
    __asm("_write$NOCANCEL");
static void *netcdf_library;
#endif

static write_function next_write(void)
{
    static write_function function;
#ifdef __APPLE__
    if (function == NULL) function = system_write_nocancel;
#else
    if (function == NULL && !resolving_write) {
        resolving_write = 1;
        function = (write_function)dlsym(RTLD_NEXT, "write");
        resolving_write = 0;
    }
#endif
    return function;
}

#ifdef __APPLE__
static void *next_netcdf_symbol(const char *name)
{
    if (netcdf_library == NULL) {
        Dl_info info;
        void *symbol;
        resolving_write = 1;
        symbol = dlsym(RTLD_DEFAULT, "nc_open");
        if (symbol != NULL && dladdr(symbol, &info) != 0 && info.dli_fname != NULL) {
            netcdf_library = dlopen(info.dli_fname, RTLD_LAZY | RTLD_LOCAL);
        }
        resolving_write = 0;
    }
    return netcdf_library == NULL ? NULL : dlsym(netcdf_library, name);
}
#else
static void *next_netcdf_symbol(const char *name)
{
    return dlsym(RTLD_NEXT, name);
}
#endif

static ssize_t raw_write(int fd, const void *buffer, size_t count)
{
#ifdef __APPLE__
    return system_write_nocancel(fd, buffer, count);
#else
    return syscall(SYS_write, fd, buffer, count);
#endif
}

static put_vara_function next_put_vara(void)
{
    static put_vara_function function;
    if (function == NULL) {
        function = (put_vara_function)next_netcdf_symbol("nc_put_vara");
    }
    return function;
}

static int descriptor_path(int fd, char *path, size_t capacity)
{
#ifdef __APPLE__
    return fcntl(fd, F_GETPATH, path) == 0 && path[0] != '\0';
#else
    char link[64];
    int length = snprintf(link, sizeof(link), "/proc/self/fd/%d", fd);
    ssize_t written;
    if (length < 0 || (size_t)length >= sizeof(link)) return 0;
    written = readlink(link, path, capacity - 1);
    if (written < 0 || (size_t)written >= capacity - 1) return 0;
    path[written] = '\0';
    return 1;
#endif
}

static ssize_t fault_write(int fd, const void *buffer, size_t count)
{
    write_function write_file;
    const char *match;
    if (resolving_write) return raw_write(fd, buffer, count);
    write_file = next_write();
    if (write_file == NULL) return raw_write(fd, buffer, count);
    match = getenv("IMEX_FAULT_WRITE_MATCH");
    if (match != NULL && match[0] != '\0') {
        char path[PATH_MAX];
        if (descriptor_path(fd, path, sizeof(path)) && strstr(path, match) != NULL) {
            int occurrence = __sync_add_and_fetch(&write_matches, 1);
            const char *selected = getenv("IMEX_FAULT_WRITE_OCCURRENCE");
            int wanted = selected == NULL ? 1 : atoi(selected);
            if (occurrence == wanted) {
                size_t partial_count = count / 2;
                ssize_t written;
                if (partial_count == 0) {
                    dprintf(STDERR_FILENO, "fault injector: target write is too small\n");
                    _exit(87);
                }
                written = write_file(fd, buffer, partial_count);
                if (written <= 0 || fsync(fd) != 0) {
                    dprintf(STDERR_FILENO, "fault injector: partial write failed: %s\n",
                            strerror(errno));
                    _exit(87);
                }
                dprintf(STDERR_FILENO,
                        "fault injector: interrupted write to %s after %zd of %zu bytes\n",
                        path, written, count);
                _exit(86);
            }
        }
    }
    return write_file(fd, buffer, count);
}

static int fault_put_vara_double(int ncid, int varid, const size_t *start,
                                 const size_t *count, const double *values)
{
    put_vara_function put = next_put_vara();
    inquire_vartype_function inquire_type =
        (inquire_vartype_function)next_netcdf_symbol("nc_inq_vartype");
    inquire_varname_function inquire =
        (inquire_varname_function)next_netcdf_symbol("nc_inq_varname");
    sync_function sync = (sync_function)next_netcdf_symbol("nc_sync");
    nc_type type;
    char name[NC_MAX_NAME + 1] = {0};
    int status;
    const char *selected = getenv("IMEX_FAULT_NETCDF_TIME_OCCURRENCE");
    if (put == NULL || inquire_type == NULL || inquire == NULL || sync == NULL) {
        dprintf(STDERR_FILENO, "fault injector: cannot resolve NetCDF functions\n");
        _exit(87);
    }
    status = inquire_type(ncid, varid, &type);
    if (status != NC_NOERR || type != NC_DOUBLE) {
        dprintf(STDERR_FILENO,
                "fault injector: generic NetCDF write requires a double variable\n");
        _exit(87);
    }
    status = inquire(ncid, varid, name);
    if (status != NC_NOERR) {
        dprintf(STDERR_FILENO, "fault injector: cannot inspect NetCDF variable\n");
        _exit(87);
    }
    /* The solver's NetCDF schema stores every variable as NC_DOUBLE. */
    status = put(ncid, varid, start, count, values);
    if (status == NC_NOERR && selected != NULL) {
        int wanted = atoi(selected);
        if (strcmp(name, "time") == 0 &&
            __sync_add_and_fetch(&netcdf_time_records, 1) == wanted) {
            if (sync(ncid) != NC_NOERR) {
                dprintf(STDERR_FILENO, "fault injector: NetCDF sync failed\n");
                _exit(87);
            }
            dprintf(STDERR_FILENO, "fault injector: interrupted after time record %d\n",
                    wanted);
            _exit(86);
        }
    }
    return status;
}

#ifdef __APPLE__
DYLD_INTERPOSE(fault_write, write);
DYLD_INTERPOSE(fault_put_vara_double, nc_put_vara_double);
#else
ssize_t write(int fd, const void *buffer, size_t count)
{
    return fault_write(fd, buffer, count);
}

int nc_put_vara_double(int ncid, int varid, const size_t *start,
                       const size_t *count, const double *values)
{
    return fault_put_vara_double(ncid, varid, start, count, values);
}
#endif
