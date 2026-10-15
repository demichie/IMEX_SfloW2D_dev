#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE

#include <errno.h>
#include <dirent.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int set_error(char *buffer, size_t capacity, const char *message)
{
    if (capacity > 0) {
        snprintf(buffer, capacity, "%s", message);
    }
    return -1;
}

static int validate_directory_contents(int directory_fd, unsigned int depth,
                                       char *error_message, size_t error_capacity)
{
    DIR *entries;
    int scan_fd;
    struct dirent *entry;

    if (depth > 64) {
        return set_error(error_message, error_capacity,
                         "existing run output directory is nested too deeply");
    }
    scan_fd = dup(directory_fd);
    if (scan_fd < 0 || (entries = fdopendir(scan_fd)) == NULL) {
        if (scan_fd >= 0) {
            close(scan_fd);
        }
        return set_error(error_message, error_capacity,
                         "cannot inspect existing run output directory");
    }

    while ((entry = readdir(entries)) != NULL) {
        struct stat info;
        int status = 0;

        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }
        if (fstatat(directory_fd, entry->d_name, &info, AT_SYMLINK_NOFOLLOW) != 0) {
            status = set_error(error_message, error_capacity,
                               "cannot inspect existing run output entry");
        } else if (S_ISLNK(info.st_mode)) {
            status = set_error(error_message, error_capacity,
                               "existing run output contains a symlink");
        } else if (S_ISDIR(info.st_mode)) {
            int child_fd = openat(directory_fd, entry->d_name,
                                  O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW);
            if (child_fd < 0) {
                status = set_error(error_message, error_capacity,
                                   "cannot safely inspect existing run output directory");
            } else {
                status = validate_directory_contents(child_fd, depth + 1,
                                                     error_message, error_capacity);
                close(child_fd);
            }
        } else if (!S_ISREG(info.st_mode) || info.st_nlink > 1) {
            status = set_error(error_message, error_capacity,
                               "existing run output contains a non-regular or multiply-linked file");
        }
        if (status != 0) {
            closedir(entries);
            return status;
        }
    }

    closedir(entries);
    return 0;
}

static int ensure_directory(const char *path, int *existed,
                            int inspect_existing_contents,
                            char *error_message, size_t error_capacity)
{
    char copy[PATH_MAX];
    char *component;
    char *saveptr = NULL;
    int current_fd;
    int absolute;

    if (path == NULL || path[0] == '\0' || strlen(path) >= sizeof(copy)) {
        return set_error(error_message, error_capacity, "invalid or overlong directory path");
    }
    strcpy(copy, path);
    absolute = copy[0] == '/';
    current_fd = open(absolute ? "/" : ".", O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (current_fd < 0) {
        return set_error(error_message, error_capacity, "cannot open starting directory");
    }

    *existed = 1;
    component = strtok_r(copy, "/", &saveptr);
    while (component != NULL) {
        int next_fd;
        struct stat info;

        if (strcmp(component, ".") == 0) {
            component = strtok_r(NULL, "/", &saveptr);
            continue;
        }
        if (strcmp(component, "..") == 0) {
            close(current_fd);
            return set_error(error_message, error_capacity,
                             "parent-directory components are not allowed");
        }
        if (fstatat(current_fd, component, &info, AT_SYMLINK_NOFOLLOW) != 0) {
            if (errno != ENOENT) {
                close(current_fd);
                return set_error(error_message, error_capacity,
                                 "cannot inspect output directory component");
            }
            *existed = 0;
            if (mkdirat(current_fd, component, 0775) != 0 && errno != EEXIST) {
                close(current_fd);
                return set_error(error_message, error_capacity,
                                 "cannot create output directory component");
            }
        } else if (S_ISLNK(info.st_mode) || !S_ISDIR(info.st_mode)) {
            close(current_fd);
            return set_error(error_message, error_capacity,
                             "output path contains a symlink or non-directory");
        }

        next_fd = openat(current_fd, component,
                         O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW);
        if (next_fd < 0) {
            close(current_fd);
            return set_error(error_message, error_capacity,
                             "cannot safely open output directory component");
        }
        close(current_fd);
        current_fd = next_fd;
        component = strtok_r(NULL, "/", &saveptr);
    }
    if (inspect_existing_contents && *existed &&
        validate_directory_contents(current_fd, 0, error_message,
                                     error_capacity) != 0) {
        close(current_fd);
        return -1;
    }
    close(current_fd);
    return 0;
}

int imex_prepare_run_directory(const char *base, const char *run_name,
                               int overwrite, int allow_existing,
                               char *error_message, size_t error_capacity)
{
    char run_path[PATH_MAX];
    size_t run_length;
    size_t i;
    int existed;
    static const char *subdirectories[] = {
        "metadata", "fields", "snapshots", "diagnostics", "probes",
        "vulnerability", "terrain", "restart"
    };

    if (base == NULL || run_name == NULL || base[0] == '\0' ||
        run_name[0] == '\0') {
        return set_error(error_message, error_capacity, "base directory and RUN_NAME must be nonempty");
    }
    for (i = 0; run_name[i] != '\0'; ++i) {
        unsigned char c = (unsigned char)run_name[i];
        if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') || c == '_' || c == '-')) {
            return set_error(error_message, error_capacity,
                             "RUN_NAME may contain only ASCII letters, digits, '_' and '-'");
        }
    }
    if (strcmp(run_name, ".") == 0 || strcmp(run_name, "..") == 0) {
        return set_error(error_message, error_capacity, "invalid RUN_NAME");
    }

    run_length = strlen(base);
    if (run_length + strlen(run_name) + 2 > sizeof(run_path)) {
        return set_error(error_message, error_capacity, "output run path is too long");
    }
    memcpy(run_path, base, run_length);
    if (run_length > 0 && base[run_length - 1] != '/') {
        run_path[run_length++] = '/';
    }
    strcpy(run_path + run_length, run_name);

    if (ensure_directory(base, &existed, 0, error_message, error_capacity) != 0) {
        return -1;
    }
    if (ensure_directory(run_path, &existed, 1, error_message, error_capacity) != 0) {
        return -1;
    }
    if (existed && !overwrite && !allow_existing) {
        return set_error(error_message, error_capacity,
                         "run output directory already exists; set OUTPUT_OVERWRITE=T or choose a new RUN_NAME");
    }

    for (i = 0; i < sizeof(subdirectories) / sizeof(subdirectories[0]); ++i) {
        char subdirectory[PATH_MAX];
        if (snprintf(subdirectory, sizeof(subdirectory), "%s/%s",
                     run_path, subdirectories[i]) >= (int)sizeof(subdirectory)) {
            return set_error(error_message, error_capacity, "output directory path is too long");
        }
        if (ensure_directory(subdirectory, &existed, 0,
                             error_message, error_capacity) != 0) {
            return -1;
        }
    }
    return 0;
}

int imex_ensure_directory(const char *path, int reject_existing,
                          char *error_message, size_t error_capacity)
{
    int existed;
    int status = ensure_directory(path, &existed, 0, error_message,
                                  error_capacity);
    if (status != 0) {
        return status;
    }
    if (reject_existing && existed) {
        return set_error(error_message, error_capacity,
                         "refusing to reuse an existing snapshot directory during binary restart");
    }
    return 0;
}

int imex_snapshot_prefix_exists(const char *directory, const char *prefix,
                                char *error_message, size_t error_capacity)
{
    DIR *entries;
    int directory_fd;
    struct dirent *entry;
    size_t prefix_length;

    if (directory == NULL || prefix == NULL || directory[0] == '\0' ||
        prefix[0] == '\0') {
        return set_error(error_message, error_capacity,
                         "invalid snapshot directory or index prefix");
    }
    directory_fd = open(directory, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW);
    if (directory_fd < 0) {
        return set_error(error_message, error_capacity,
                         "cannot inspect snapshot output directory");
    }
    entries = fdopendir(directory_fd);
    if (entries == NULL) {
        close(directory_fd);
        return set_error(error_message, error_capacity,
                         "cannot scan snapshot output directory");
    }

    prefix_length = strlen(prefix);
    while ((entry = readdir(entries)) != NULL) {
        if (strncmp(entry->d_name, prefix, prefix_length) == 0) {
            closedir(entries);
            return 1;
        }
    }
    closedir(entries);
    return 0;
}
