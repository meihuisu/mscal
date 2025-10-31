

// g++ interp.cpp alglibinternal.cpp interpolation.cpp alglibmisc.cpp ap.cpp -lnetcdf -lm -o interp
// 


#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "alglib.h"

using namespace alglib;

#define DEP_SCALE_FACTOR 111000.0

// Error handling macro
#define NC_CHECK(call) \
    if ((call) != NC_NOERR) { \
        fprintf(stderr, "NetCDF error: %s\n", nc_strerror(call)); \
        exit(EXIT_FAILURE); \
    }

// Load 1D variable from NetCDF
std::vector<double> load_nc_1d(const char* filepath, const char* varname, size_t* len) {
    int ncid, varid;
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));
    NC_CHECK(nc_inq_dimlen(ncid, varid, len));

    std::vector<double> data(*len);
    NC_CHECK(nc_get_var_double(ncid, varid, data.data()));
    NC_CHECK(nc_close(ncid));
    return data;
}

// Load 3D variable from NetCDF
std::vector<double> load_nc_3d(const char* filepath, const char* varname, size_t nx, size_t ny, size_t nz) {
    int ncid, varid;
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));

    std::vector<double> data(nx * ny * nz);
    NC_CHECK(nc_get_var_double(ncid, varid, data.data()));
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
    auto loni = load_nc_1d(ncfpath, "longitude", &nx);
    auto lati = load_nc_1d(ncfpath, "latitude", &ny);
    auto depi = load_nc_1d(ncfpath, "depth", &nz);

    auto vpi = load_nc_3d(ncfpath, "vp", nx, ny, nz);
    auto vsi = load_nc_3d(ncfpath, "vs", nx, ny, nz);

    size_t nxyz = nx * ny * nz;

    // Prepare input for ALGLIB
    real_2d_array xy;
    real_2d_array y;
    xy.setlength(nxyz, 3);
    y.setlength(nxyz, 2);

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                size_t idx = (k * ny + j) * nx + i;
                xy[idx][0] = loni[i];
                xy[idx][1] = lati[j];
                xy[idx][2] = depi[k] / DEP_SCALE_FACTOR;
                y[idx][0] = vpi[idx];
                y[idx][1] = vsi[idx];
            }
        }
    }

    // Build RBF model
    rbfmodel model;
    rbfreport rep;
    rbfcreate(3, 2, model);
    rbfsetpoints(model, xy, y);
    rbfbuildmodel(model, 1.0, rep);

    // Query point
    real_1d_array query;
    query.setlength(3);
    query[0] = lonq;
    query[1] = latq;
    query[2] = depq;

    real_1d_array result = rbfcalc(model, query);
    printf("Interpolated Vp = %.2f m/s, Vs = %.2f m/s\n", result[0], result[1]);

    return 0;
}

