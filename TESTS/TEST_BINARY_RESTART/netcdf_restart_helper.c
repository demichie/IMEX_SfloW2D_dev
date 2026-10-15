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

static int verify_partial_time_record(const char *path, double expected_time)
{
    int ncid, time_dimid, time_varid, field_varid, status;
    int field_rank, field_dimids[NC_MAX_VAR_DIMS], time_axis = -1;
    size_t time_length, start[NC_MAX_VAR_DIMS], count[NC_MAX_VAR_DIMS];
    size_t record_values = 1;
    size_t axis, i;
    double *times = NULL, *field_values = NULL;
    double fill_value = NC_FILL_DOUBLE;
    char dimension_name[NC_MAX_NAME + 1];
    int result = 1, has_written_field_value = 0;

    status = nc_open(path, NC_NOWRITE, &ncid);
    if (report_netcdf_error(status, "open partially written dataset")) return 1;
    status = nc_inq_dimid(ncid, "time", &time_dimid);
    if (report_netcdf_error(status, "find partial time dimension")) goto done;
    status = nc_inq_dimlen(ncid, time_dimid, &time_length);
    if (report_netcdf_error(status, "read partial time dimension")) goto done;
    if (time_length != 3) {
        fprintf(stderr, "expected 3 time records after injected interruption, got %zu\n",
                time_length);
        goto done;
    }
    status = nc_inq_varid(ncid, "time", &time_varid);
    if (report_netcdf_error(status, "find partial time variable")) goto done;
    times = malloc(time_length * sizeof(*times));
    if (times == NULL) {
        fprintf(stderr, "cannot allocate partial time buffer\n");
        goto done;
    }
    status = nc_get_var_double(ncid, time_varid, times);
    if (report_netcdf_error(status, "read partial time records")) goto done;
    if (!isfinite(times[2]) || fabs(times[2] - expected_time) > 1.e-12) {
        fprintf(stderr, "unexpected interrupted time record: %.17g vs %.17g\n",
                times[2], expected_time);
        goto done;
    }

    status = nc_inq_varid(ncid, "h", &field_varid);
    if (report_netcdf_error(status, "find partial field variable")) goto done;
    status = nc_inq_var(ncid, field_varid, NULL, NULL, &field_rank,
                        field_dimids, NULL);
    if (report_netcdf_error(status, "inspect partial field variable")) goto done;
    if (field_rank < 1 || field_rank > NC_MAX_VAR_DIMS) {
        fprintf(stderr, "invalid partial field rank: %d\n", field_rank);
        goto done;
    }
    for (axis = 0; axis < (size_t)field_rank; ++axis) {
        size_t dimension_length;
        status = nc_inq_dim(ncid, field_dimids[axis], dimension_name,
                            &dimension_length);
        if (report_netcdf_error(status, "inspect partial field dimension")) goto done;
        start[axis] = 0;
        count[axis] = dimension_length;
        if (strcmp(dimension_name, "time") == 0) {
            if (time_axis >= 0) {
                fprintf(stderr, "partial field has multiple time dimensions\n");
                goto done;
            }
            time_axis = (int)axis;
            count[axis] = 1;
        } else {
            if (dimension_length != 0 &&
                record_values > SIZE_MAX / dimension_length) {
                fprintf(stderr, "partial field record exceeds addressable size\n");
                goto done;
            }
            record_values *= dimension_length;
        }
    }
    if (time_axis < 0 || record_values == 0 ||
        record_values > SIZE_MAX / sizeof(*field_values)) {
        fprintf(stderr, "partial field has no usable time slice\n");
        goto done;
    }
    status = nc_get_att_double(ncid, field_varid, "_FillValue", &fill_value);
    if (status != NC_NOERR && status != NC_ENOTATT) {
        report_netcdf_error(status, "read partial field fill value");
        goto done;
    }
    field_values = malloc(record_values * sizeof(*field_values));
    if (field_values == NULL) {
        fprintf(stderr, "cannot allocate partial field buffer\n");
        goto done;
    }
    status = nc_get_vara_double(ncid, field_varid, start, count, field_values);
    if (report_netcdf_error(status, "read completed field record")) goto done;
    for (i = 0; i < record_values; ++i) {
        if (!isfinite(field_values[i])) {
            fprintf(stderr, "completed field record contains a non-finite value\n");
            goto done;
        }
        if (field_values[i] != fill_value) has_written_field_value = 1;
    }
    if (!has_written_field_value) {
        fprintf(stderr, "completed field record contains only fill values\n");
        goto done;
    }

    start[time_axis] = time_length - 1;
    status = nc_get_vara_double(ncid, field_varid, start, count, field_values);
    if (report_netcdf_error(status, "read interrupted field record")) goto done;
    for (i = 0; i < record_values; ++i) {
        if (field_values[i] != fill_value) {
            fprintf(stderr, "interrupted field record contains written data at %zu\n", i);
            goto done;
        }
    }
    result = 0;

done:
    free(field_values);
    free(times);
    nc_close(ncid);
    return result;
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
    if (argc == 4 && strcmp(argv[1], "partial-check") == 0) {
        time_value = strtod(argv[3], &end);
        if (end == argv[3] || *end != '\0' || !isfinite(time_value)) {
            fprintf(stderr, "invalid expected partial time value\n");
            return 2;
        }
        return verify_partial_time_record(argv[2], time_value);
    }
    if (argc == 4 && strcmp(argv[1], "compare") == 0) {
        return compare_outputs(argv[2], argv[3]);
    }
    fprintf(stderr, "Usage: %s partial FILE TIME | partial-check FILE TIME | "
            "compare REFERENCE RESUMED\n",
            argv[0]);
    return 2;
}
