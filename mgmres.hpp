void atx_cr ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[], 
  double w[] );
void atx_st ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[], 
  double w[] );
void ax_cr ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[], 
  double w[] );
void ax_st ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[], 
  double w[] );
void diagonal_pointer_cr ( long int n, long int nz_num, long int ia[], long int ja[], long int ua[] );
void lus_cr ( long int n, long int nz_num, long int ia[], long int ja[], double l[], long int ua[], 
  double r[], double z[] );
void mgmres_st ( long int n, long int nz_num, long int ia[], long int ja[], double a[], double x[],
  double rhs[], long int itr_max, long int mr, double tol_abs, double tol_rel );
void mult_givens ( double c, double s, long int k, double g[] );
void pmgmres_ilu_cr ( long int n, long int nz_num, long int ia[], long int ja[], double a[], 
  double x[], double rhs[], long int itr_max, long int mr, double tol_abs, 
  double tol_rel, short int verbose=0 );
double r8vec_dot ( long int n, double a1[], double a2[] );
double *r8vec_uniform_01 ( long int n, long int *seed );
void rearrange_cr ( long int n, long int nz_num, long int ia[], long int ja[], double a[] );
void timestamp ( );
