//
//
// What You Need
// ALGLIB C source files: alglib.h, ap.c, interpolation.c, alglibinternal.c
// NetCDF C library: Install via libnetcdf-dev (Debian/Ubuntu) or equivalent
// Compile with:
//    gcc rbf_netCDF.c ap.c interpolation.c alglibinternal.c -lnetcdf -lm -o rbf_interp
//
//
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <math.h>
#include "alglib.h"

#define DEP_SCALE_FACTOR 111000.0

// Error handling macro
#define NC_CHECK(call) \
    if ((call) != NC_NOERR) { \
        fprintf(stderr, "NetCDF error: %s\n", nc_strerror(call)); \
        exit(EXIT_FAILURE); \
    }

// Load 1D variable from NetCDF
double* load_nc_1d(const char* filepath, const char* varname, size_t* len) {
    int ncid, varid;
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));
    NC_CHECK(nc_inq_dimlen(ncid, varid, len));

    double* data = (double*)malloc(*len * sizeof(double));
    NC_CHECK(nc_get_var_double(ncid, varid, data));
    NC_CHECK(nc_close(ncid));
    return data;
}

// Load 3D variable from NetCDF
double* load_nc_3d(const char* filepath, const char* varname, size_t nx, size_t ny, size_t nz) {
    int ncid, varid;
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));

    double* data = (double*)malloc(nx * ny * nz * sizeof(double));
    NC_CHECK(nc_get_var_double(ncid, varid, data));
    NC_CHECK(nc_close(ncid));
    return data;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <netcdf_file> <lon> <lat> <depth>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* ncfpath = argv[1];
    double lonq = atof(argv[2]);
    double latq = atof(argv[3]);
    double depq = atof(argv[4]) / DEP_SCALE_FACTOR;

    size_t nx, ny, nz;
    double* loni = load_nc_1d(ncfpath, "longitude", &nx);
    double* lati = load_nc_1d(ncfpath, "latitude", &ny);
    double* depi = load_nc_1d(ncfpath, "depth", &nz);

    double* vpi = load_nc_3d(ncfpath, "vp", nx, ny, nz);
    double* vsi = load_nc_3d(ncfpath, "vs", nx, ny, nz);

    ae_state _state;
    ae_init_pool(&_state);

    size_t nxyz = nx * ny * nz;
    ae_matrix xy, y;
    ae_matrix_init(&xy, nxyz, 3, DT_REAL, &_state);
    ae_matrix_init(&y, nxyz, 2, DT_REAL, &_state);

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                size_t idx = (k * ny + j) * nx + i;
                xy.ptr.p_double[idx][0] = loni[i];
                xy.ptr.p_double[idx][1] = lati[j];
                xy.ptr.p_double[idx][2] = depi[k] / DEP_SCALE_FACTOR;
                y.ptr.p_double[idx][0] = vpi[idx];
                y.ptr.p_double[idx][1] = vsi[idx];
            }
        }
    }

    rbfmodel model;
    rbfreport rep;
    rbfcreate(3, 2, &model, &_state);
    rbfsetpoints(&model, &xy, &y, &_state);
    rbfbuildmodel(&model, 1.0, &rep, &_state);

    ae_vector query, result;
    ae_vector_init(&query, 3, DT_REAL, &_state);
    ae_vector_init(&result, 2, DT_REAL, &_state);

    query.ptr.p_double[0] = lonq;
    query.ptr.p_double[1] = latq;
    query.ptr.p_double[2] = depq;

    rbfcalc(&model, &query, &result, &_state);
    printf("Interpolated Vp = %.2f m/s, Vs = %.2f m/s\n", result.ptr.p_double[0], result.ptr.p_double[1]);

    ae_free_pool(&_state);
    free(loni); free(lati); free(depi);
    free(vpi);  free(vsi);
    return 0;
}
