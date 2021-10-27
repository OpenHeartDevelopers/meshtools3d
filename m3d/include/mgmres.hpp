#ifdef CGAL_LINKED_WITH_TBB
#ifndef USE_TBB_PARALLEL
#define USE_TBB_PARALLEL 1
#endif
#endif
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
void rearrange_cr_row ( const long int & j1,const long int & j2, long int * ja, double * a );
void timestamp ( );




#ifdef USE_TBB_PARALLEL

#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <tbb/tbb.h>

class r8vec_dot_parallel
{
  private:
    const double *inputArray1;
    const double *inputArray2;
    double  sum;
  public:
    r8vec_dot_parallel()
    :inputArray1(NULL),
    inputArray2(NULL),
    sum(0.0)
    {}
    r8vec_dot_parallel(const double *iArray1, const double *iArray2)
    :inputArray1(iArray1),
    inputArray2(iArray2),
    sum(0.0)
    {}
	  r8vec_dot_parallel(r8vec_dot_parallel & pvecdot)
    :inputArray1(pvecdot.inputArray1),
    inputArray2(pvecdot.inputArray2),
    sum(0.0)
    {}
	  r8vec_dot_parallel(r8vec_dot_parallel & pvecdot,tbb::split)
    :inputArray1(pvecdot.inputArray1),
    inputArray2(pvecdot.inputArray2),
    sum(0.0)
    {}
	  double getSum(){return sum;}
	  
    void operator()(const tbb::blocked_range<size_t> &r)
	  {
        for(size_t count = r.begin(); count != r.end(); count++)
	      {
	          sum += inputArray1[count]*inputArray2[count];
	      }
	  }
    void join(r8vec_dot_parallel & pvecdotSub)
    {
        sum  += pvecdotSub.sum;
	  }
};// r8vec_dot_parallel

struct diagonal_pointer_cr_parallel {
    const long int * inputArrayRow;
    const long int * inputArrayCol;
    long int * outputArray;
    void operator()(tbb::blocked_range<size_t> &r) const
    {
      for(size_t count = r.begin(); count != r.end(); count++)
	    {
	      outputArray[count]=-1;
	      long int j=0;
        for(j=inputArrayRow[count]; j<inputArrayRow[count+1]; j++)
        {
          if (inputArrayCol[j] == count) 
          {
             outputArray[count] = j;
             break;
          }
        }
	    }
    }
};//diagonal_pointer_cr_parallel

template<typename T>
struct copy_vector_parallel{
  const T * inputArray;
  T * outputArray;
  void operator()(tbb::blocked_range<size_t> &r) const
  {
    for(size_t count = r.begin(); count != r.end(); count++)
    {
        outputArray[count]=inputArray[count];
    }
  }
};//copy_vector_parallel
	 
template<typename T> void vectorCopy(const T * source, T* target, long int size)
{
  copy_vector_parallel<T> vcp;
  vcp.inputArray=source;
  vcp.outputArray=target;
  tbb::parallel_for(tbb::blocked_range<size_t>(0,static_cast<size_t>(size)), vcp);
};

template<typename T>
struct initialize_vector_parallel{
  T * array;
  T value;
  void operator()(tbb::blocked_range<size_t> &r) const
  {
    for(size_t count = r.begin(); count != r.end(); count++)
    {
        array[count]=value;
    }
  }
};//initialize_vector_parallel

template<typename T> void initializeVector( T * Vec, T val, long int size)
{
  initialize_vector_parallel<T> ini_p;
  ini_p.array=Vec;
  ini_p.value=val;
  tbb::parallel_for(tbb::blocked_range<size_t>(0,static_cast<size_t>(size)), ini_p);
};

struct ax_cr_parallel{
  double * resArray;
  const double * aij;
  const double * xA;
  const long int * I;  
  const long int * J;
  
  void operator()(tbb::blocked_range<size_t> &r) const
  {
      for(size_t count = r.begin(); count != r.end(); count++)
      {
        double sum=0.0;
        for (long int k = I[count]; k < I[count+1]; k++ )      
        {
            sum=sum+ aij[k]*xA[J[k]];
        }
        resArray[count]=sum;
      }
  }
};//ax_cr_parallel

template<typename T>
struct vxmv{
  T * updatedArray;
  const T * updatingArray;
  void operator()(tbb::blocked_range<size_t> &r) const
  {
      for(size_t count = r.begin(); count != r.end(); count++)
      {
        updatedArray[count]=updatingArray[count]-updatedArray[count];
      }

  }
};//vxmv

template<typename T> void vxmvVectorUpdating(T * IOarray, const T * updatingArray, long int size)
{
    vxmv<T> vxmv_op;
    vxmv_op.updatedArray=IOarray;
    vxmv_op.updatingArray=updatingArray;
    tbb::parallel_for(tbb::blocked_range<size_t>(0,static_cast<size_t>(size)), vxmv_op);
};

struct xxvy{
  double * xV;
  const double * yV;
  const double * vCM;
  long int k;
  long int n;
  void operator()(tbb::blocked_range<size_t> &r) const
  {
      for(size_t count = r.begin(); count != r.end(); count++)
      {
          for (long int j = 0; j < k + 1; j++ )
          {
              xV[count] = xV[count] + vCM[j*n+count] * yV[j];
          }
      }
  }
};//xxvy


void xxvyFunc(double * xVArray, const double * yVArray, const double * vCMmat, long int k_index, long int size);

template<typename T>
struct yaxs{
  const T * sourceArray;
  T* targetArray;
  T alpha;
  void operator()(tbb::blocked_range<size_t> &r) const
  {
      for(size_t count = r.begin(); count != r.end(); count++)
      {
          targetArray[count]=alpha*sourceArray[count];
      }
  }
};

template<typename T> void yaxsFunc(const T * source, T* target, T scalval,long int size)
{
  yaxs<T> yaxs_op;
  yaxs_op.sourceArray=source;
  yaxs_op.targetArray=target;
  yaxs_op.alpha=scalval;
  tbb::parallel_for(tbb::blocked_range<size_t>(0,static_cast<size_t>(size)), yaxs_op);
}

#endif











