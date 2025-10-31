
/***
   alglib_rbf_wrapper.cpp
***/

#include "interpolation.h"

extern "C" {

typedef void* RBFHandle;

// Create and build an RBF model
RBFHandle create_rbf_model(const double* xy_data, const double* y_data, int n_points, int dim, double rbf_radius) {
    try {
        alglib::real_2d_array xy;
        alglib::real_1d_array y;
        alglib::rbfmodel* model = new alglib::rbfmodel();
        alglib::rbfreport rep;

        xy.setlength(n_points, dim);
        y.setlength(n_points);

        for (int i = 0; i < n_points; ++i) {
            for (int j = 0; j < dim; ++j) {
                xy[i][j] = xy_data[i * dim + j];
            }
            y[i] = y_data[i];
        }

        alglib::rbfcreate(dim, 1, *model);
        alglib::rbfsetpoints(*model, xy, y);
        alglib::rbfbuildmodel(*model, rbf_radius, rep);

        return static_cast<void*>(model);
    } catch (...) {
        return nullptr;
    }
}

// Evaluate the RBF model at a query point
int eval_rbf_model(RBFHandle handle, const double* query, int dim, double* result_out) {
    if (!handle) return -1;
    try {
        alglib::rbfmodel* model = static_cast<alglib::rbfmodel*>(handle);
        alglib::real_1d_array query_point;
        query_point.setlength(dim);
        for (int i = 0; i < dim; ++i) {
            query_point[i] = query[i];
        }

        alglib::real_1d_array result = alglib::rbfcalc(*model, query_point);
        for (int i = 0; i < result.length(); ++i) {
            result_out[i] = result[i];
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

// Free the RBF model
void free_rbf_model(RBFHandle handle) {
    if (handle) {
        delete static_cast<alglib::rbfmodel*>(handle);
    }
}

}

