
/***
   alglib_rbf_wrapper.h
***/

typedef void* RBFHandle;

extern RBFHandle create_rbf_model(const double* xy_data, const double* y_data, int n_points, int dim, double rbf_radius);
extern int eval_rbf_model(RBFHandle handle, const double* query, int dim, double* result_out);
extern void free_rbf_model(RBFHandle handle);
