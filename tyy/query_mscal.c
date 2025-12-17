//
// gcc -o interp interp.c -lnetcdf -lm
//


#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <math.h>

#define DEP_SCALE_FACTOR 111000.0

// Error handling macro
#define NC_CHECK(call) \
    if ((call) != NC_NOERR) { \
        fprintf(stderr, "NetCDF error: %s\n", nc_strerror(call)); \
        exit(EXIT_FAILURE); \
    }

// Load 1D variable from NetCDF file
double* load_ncdata(const char* filepath, const char* varname, size_t* len) {
    int ncid, varid;
    NC_CHECK(nc_open(filepath, NC_NOWRITE, &ncid));
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));
    NC_CHECK(nc_inq_dimlen(ncid, varid, len));

    double* data = (double*)malloc(*len * sizeof(double));
    NC_CHECK(nc_get_var_double(ncid, varid, data));
    NC_CHECK(nc_close(ncid));
    return data;
}

// Simple trilinear interpolation (placeholder)
double interpolate_mod(double* loni, double* lati, double* depi, double* vari,
                       size_t nx, size_t ny, size_t nz,
                       double lonq, double latq, double depq) {
// XXX 

    // This is a placeholder. Real implementation would require 3D grid interpolation.
    // For now, return the first value as a dummy.
    return vari[0];
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
    double* loni = load_ncdata(ncfpath, "longitude", &nx);
    double* lati = load_ncdata(ncfpath, "latitude", &ny);
    double* depi = load_ncdata(ncfpath, "depth", &nz);

    size_t nxyz = nx * ny * nz;
    double* vpi = load_ncdata(ncfpath, "vp", &nxyz);
    double* vsi = load_ncdata(ncfpath, "vs", &nxyz);
    double* rhoi = load_ncdata(ncfpath, "rho", &nxyz);

    double vpq = interpolate_mod(loni, lati, depi, vpi, nx, ny, nz, lonq, latq, depq);
    double vsq = interpolate_mod(loni, lati, depi, vsi, nx, ny, nz, lonq, latq, depq);
    double rhoq = interpolate_mod(loni, lati, depi, rhoi, nx, ny, nz, lonq, latq, depq);

    printf("Vp, Vs, Density = %.2f m/s, %.2f m/s, %.2f kg/m^3\n", vpq, vsq, rhoq);

    free(loni); free(lati); free(depi);
    free(vpi); free(vsi); free(rhoi);
    return 0;
}


double trilinear_interp(double* grid, size_t nx, size_t ny, size_t nz,
                        double lonq, double latq, double depq,
                        double lon_min, double lon_step,
                        double lat_min, double lat_step,
                        double dep_min, double dep_step) {
    // Compute indices
    int i = (int)((lonq - lon_min) / lon_step);
    int j = (int)((latq - lat_min) / lat_step);
    int k = (int)((depq - dep_min) / dep_step);

    // Clamp indices to valid range
    if (i < 0 || i >= nx - 1 || j < 0 || j >= ny - 1 || k < 0 || k >= nz - 1) {
        return NAN;
    }

    // Compute fractional distances
    double xd = (lonq - (lon_min + i * lon_step)) / lon_step;
    double yd = (latq - (lat_min + j * lat_step)) / lat_step;
    double zd = (depq - (dep_min + k * dep_step)) / dep_step;

    // Get values at the 8 corners of the cube
    size_t idx = (k * ny + j) * nx + i;
    double c000 = grid[idx];
    double c001 = grid[idx + nx * ny];
    double c010 = grid[idx + nx];
    double c011 = grid[idx + nx + nx * ny];
    double c100 = grid[idx + 1];
    double c101 = grid[idx + 1 + nx * ny];
    double c110 = grid[idx + 1 + nx];
    double c111 = grid[idx + 1 + nx + nx * ny];

    // Interpolate
    double c00 = c000 * (1 - xd) + c100 * xd;
    double c01 = c001 * (1 - xd) + c101 * xd;
    double c10 = c010 * (1 - xd) + c110 * xd;
    double c11 = c011 * (1 - xd) + c111 * xd;

    double c0 = c00 * (1 - yd) + c10 * yd;
    double c1 = c01 * (1 - yd) + c11 * yd;

    double c = c0 * (1 - zd) + c1 * zd;

    return c;
}
