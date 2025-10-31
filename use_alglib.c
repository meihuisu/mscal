#include <stdio.h>

typedef void* RBFHandle;

extern RBFHandle create_rbf_model(const double* xy_data, const double* y_data, int n_points, int dim, double rbf_radius);
extern int eval_rbf_model(RBFHandle handle, const double* query, int dim, double* result_out);
extern void free_rbf_model(RBFHandle handle);

int main() {
    double xy[] = {
        0.0, 0.0, 0.0,
        1.0, 1.0, 1.0,
        2.0, 2.0, 2.0
    };
    double y[] = {1.0, 2.0, 3.0};
    double query[] = {1.5, 1.5, 1.5};
    double result[1];

    RBFHandle model = create_rbf_model(xy, y, 3, 3, 1.0);
    if (model) {
        if (eval_rbf_model(model, query, 3, result) == 0) {
            printf("RBF result: %f\n", result[0]);
        } else {
            printf("Failed to evaluate RBF model.\n");
        }
        free_rbf_model(model);
    } else {
        printf("Failed to create RBF model.\n");
    }

    return 0;
}

