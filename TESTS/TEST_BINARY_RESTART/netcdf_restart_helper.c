#include <math.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int report_netcdf_error(int status, const char *operation)
{
    if (status == NC_NOERR) {
        return 0;
    }
    fprintf(stderr, "%s: %s\n", operation, nc_strerror(status));
    return 1;
}

static int append_partial_time_record(const char *path, double time_value)
{
    int ncid, dimid, varid, status;
    size_t length;
    size_t start[1], count[1] = {1};

    status = nc_open(path, NC_WRITE, &ncid);
    if (report_netcdf_error(status, "open")) return 1;
    status = nc_inq_dimid(ncid, "time", &dimid);
    if (!report_netcdf_error(status, "find time dimension")) {
        status = nc_inq_dimlen(ncid, dimid, &length);
    }
    if (report_netcdf_error(status, "read time dimension")) {
        nc_close(ncid);
        return 1;
    }
    status = nc_inq_varid(ncid, "time", &varid);
    if (report_netcdf_error(status, "find time variable")) {
        nc_close(ncid);
        return 1;
    }

    start[0] = length;
    status = nc_put_vara_double(ncid, varid, start, count, &time_value);
    if (report_netcdf_error(status, "write partial time record")) {
        nc_close(ncid);
        return 1;
    }
    status = nc_close(ncid);
    return report_netcdf_error(status, "close");
}

static int compare_variable_data(int left_ncid, int right_ncid,
                                 const char *name, size_t value_count)
{
    int left_varid, right_varid, status;
    double *left_values, *right_values;
    size_t bytes = value_count * sizeof(double);

    status = nc_inq_varid(left_ncid, name, &left_varid);
    if (report_netcdf_error(status, "find variable in reference")) return 1;
    status = nc_inq_varid(right_ncid, name, &right_varid);
    if (report_netcdf_error(status, "find variable in resumed output")) return 1;

    left_values = malloc(bytes);
    right_values = malloc(bytes);
    if (left_values == NULL || right_values == NULL) {
        fprintf(stderr, "cannot allocate NetCDF comparison buffers\n");
        free(left_values);
        free(right_values);
        return 1;
    }
    status = nc_get_var_double(left_ncid, left_varid, left_values);
    if (report_netcdf_error(status, "read reference variable")) {
        free(left_values);
        free(right_values);
        return 1;
    }
    status = nc_get_var_double(right_ncid, right_varid, right_values);
    if (report_netcdf_error(status, "read resumed variable")) {
        free(left_values);
        free(right_values);
        return 1;
    }

    if (memcmp(left_values, right_values, bytes) != 0) {
        size_t i;
        for (i = 0; i < value_count; ++i) {
            if (memcmp(&left_values[i], &right_values[i], sizeof(double)) != 0) {
                fprintf(stderr, "NetCDF variable differs: %s at %zu: %.17g vs %.17g\n",
                        name, i, left_values[i], right_values[i]);
                break;
            }
        }
        free(left_values);
        free(right_values);
        return 1;
    }
    free(left_values);
    free(right_values);
    return 0;
}

static int compare_outputs(const char *reference_path, const char *resumed_path)
{
    int reference_ncid, resumed_ncid, status, result = 1;
    int x_dimid, y_dimid, time_dimid;
    size_t x_length, y_length, time_length, reference_time_length;
    size_t field_count;

    status = nc_open(reference_path, NC_NOWRITE, &reference_ncid);
    if (report_netcdf_error(status, "open reference")) return 1;
    status = nc_open(resumed_path, NC_NOWRITE, &resumed_ncid);
    if (report_netcdf_error(status, "open resumed output")) {
        nc_close(reference_ncid);
        return 1;
    }

    status = nc_inq_dimid(reference_ncid, "x", &x_dimid);
    if (report_netcdf_error(status, "find x dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(reference_ncid, x_dimid, &x_length),
                            "read x dimension")) goto done;
    status = nc_inq_dimid(reference_ncid, "y", &y_dimid);
    if (report_netcdf_error(status, "find y dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(reference_ncid, y_dimid, &y_length),
                            "read y dimension")) goto done;
    status = nc_inq_dimid(reference_ncid, "time", &time_dimid);
    if (report_netcdf_error(status, "find time dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(reference_ncid, time_dimid,
                                           &reference_time_length),
                            "read reference time dimension")) goto done;
    if (reference_time_length == 0) {
        fprintf(stderr, "reference NetCDF has no time records\n");
        goto done;
    }

    status = nc_inq_dimid(resumed_ncid, "x", &x_dimid);
    if (report_netcdf_error(status, "find resumed x dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(resumed_ncid, x_dimid, &time_length),
                            "read resumed x dimension")) goto done;
    if (time_length != x_length) {
        fprintf(stderr, "NetCDF x dimensions differ\n");
        goto done;
    }
    status = nc_inq_dimid(resumed_ncid, "y", &y_dimid);
    if (report_netcdf_error(status, "find resumed y dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(resumed_ncid, y_dimid, &time_length),
                            "read resumed y dimension")) goto done;
    if (time_length != y_length) {
        fprintf(stderr, "NetCDF y dimensions differ\n");
        goto done;
    }
    status = nc_inq_dimid(resumed_ncid, "time", &time_dimid);
    if (report_netcdf_error(status, "find resumed time dimension")) goto done;
    if (report_netcdf_error(nc_inq_dimlen(resumed_ncid, time_dimid, &time_length),
                            "read resumed time dimension")) goto done;
    if (time_length != reference_time_length) {
        fprintf(stderr, "NetCDF time dimensions differ: %zu vs %zu\n",
                reference_time_length, time_length);
        goto done;
    }

    if (x_length > SIZE_MAX / y_length ||
        x_length * y_length > SIZE_MAX / time_length) {
        fprintf(stderr, "NetCDF dimensions exceed comparison capacity\n");
        goto done;
    }
    field_count = x_length * y_length * time_length;
    if (compare_variable_data(reference_ncid, resumed_ncid, "time", time_length) ||
        compare_variable_data(reference_ncid, resumed_ncid, "h", field_count)) {
        goto done;
    }
    result = 0;

done:
    nc_close(resumed_ncid);
    nc_close(reference_ncid);
    return result;
}

int main(int argc, char **argv)
{
    char *end;
    double time_value;

    if (argc == 4 && strcmp(argv[1], "partial") == 0) {
        time_value = strtod(argv[3], &end);
        if (end == argv[3] || *end != '\0' || !isfinite(time_value)) {
            fprintf(stderr, "invalid partial time value\n");
            return 2;
        }
        return append_partial_time_record(argv[2], time_value);
    }
    if (argc == 4 && strcmp(argv[1], "compare") == 0) {
        return compare_outputs(argv[2], argv[3]);
    }
    fprintf(stderr, "Usage: %s partial FILE TIME | compare REFERENCE RESUMED\n",
            argv[0]);
    return 2;
}
