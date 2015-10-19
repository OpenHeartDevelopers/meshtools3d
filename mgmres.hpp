void atx_cr ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[],   double w[] );
void atx_st ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[],   double w[] );
void ax_cr ( const long int & n, const long int & nz_num, const long int * ia, const long int * ja, const double * a, const double *x,  double * w );
void ax_st ( long int n, long int nz_num, double a[], long int ia[], long int ja[], double x[],   double w[] );
void diagonal_pointer_cr ( const long int & n, const long int & nz_num, const long int * ia, const long int * ja, long int * ua );
void lus_cr ( const long int & n, const long int & nz_num, const long int *ia, const long int *ja, const double *l, const long int *ua,   const double *r, double * z );
void mgmres_st ( long int n, long int nz_num, long int ia[], long int ja[], double a[], double x[],  double rhs[], long int itr_max, long int mr, double tol_abs, double tol_rel );

void mult_givens ( const double & c, const double & s, const long int & k, double * g );

void pmgmres_ilu_cr ( const long int & n, const long int & nz_num, 
                      const long int * ia,  long int * ja,  double *a,   double * x, const double *rhs, 
                      const long int & itr_max, const long int & mr, const double & tol_abs, 
                      const double & tol_rel, short int verbose=0 );

double r8vec_dot ( const long int & n, const double * a1, const double * a2 );
double *r8vec_uniform_01 ( const long int & n, long int *seed );
void rearrange_cr ( const long int & n, const long int & nz_num, const long int * ia, long int * ja, double * a );
void timestamp ( );
