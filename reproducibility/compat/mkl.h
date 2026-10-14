// Narrow, audited qbasis compatibility layer, NOT an implementation of MKL.
// Dense BLAS/LAPACK calls go directly to OpenBLAS/LAPACKE (LP64).
// Only the zero-based, non-transposed CSR operations used by pinned qbasis
// are supported. Unsupported operations return an error, never a wrong answer.
#pragma once
#include <complex>
#include <cstdint>
#include <type_traits>
#include <stdexcept>
using MKL_INT = int;
using MKL_INT64 = long long;
#include <cblas.h>
#include <lapacke.h>
using CBLAS_LAYOUT = CBLAS_ORDER;
enum sparse_status_t { SPARSE_STATUS_SUCCESS, SPARSE_STATUS_NOT_SUPPORTED, SPARSE_STATUS_INVALID_VALUE };
enum sparse_index_base_t { SPARSE_INDEX_BASE_ZERO };
enum sparse_operation_t { SPARSE_OPERATION_NON_TRANSPOSE };
enum sparse_matrix_type_t { SPARSE_MATRIX_TYPE_GENERAL, SPARSE_MATRIX_TYPE_SYMMETRIC, SPARSE_MATRIX_TYPE_HERMITIAN };
enum sparse_fill_mode_t { SPARSE_FILL_MODE_UPPER };
enum sparse_diag_type_t { SPARSE_DIAG_NON_UNIT };
struct matrix_descr { sparse_matrix_type_t type; sparse_fill_mode_t mode; sparse_diag_type_t diag; };
struct mmhk_csr { int rows, cols; int *begin, *end, *index; void *values; bool complex; };
using sparse_matrix_t = mmhk_csr*;
template<class T> inline sparse_status_t mmhk_create(sparse_matrix_t* A, sparse_index_base_t base,
 int rows,int cols,int* begin,int* end,int* index,T* values) {
 if(base!=SPARSE_INDEX_BASE_ZERO || rows<0 || cols<0 || !A) return SPARSE_STATUS_INVALID_VALUE;
 *A=new mmhk_csr{rows,cols,begin,end,index,values,!std::is_same_v<T,double>};
 return SPARSE_STATUS_SUCCESS;
}
inline sparse_status_t mkl_sparse_d_create_csr(sparse_matrix_t* A,sparse_index_base_t b,int r,int c,int* s,int* e,int* j,double* v) { return mmhk_create(A,b,r,c,s,e,j,v); }
inline sparse_status_t mkl_sparse_z_create_csr(sparse_matrix_t* A,sparse_index_base_t b,int r,int c,int* s,int* e,int* j,std::complex<double>* v) { return mmhk_create(A,b,r,c,s,e,j,v); }
inline sparse_status_t mkl_sparse_destroy(sparse_matrix_t A) { delete A; return SPARSE_STATUS_SUCCESS; }
inline double mmhk_conj(double x) { return x; }
inline std::complex<double> mmhk_conj(std::complex<double> x) { return std::conj(x); }
template<class T> inline sparse_status_t mmhk_mv(sparse_operation_t op,T alpha,sparse_matrix_t A,
 matrix_descr d,const T* x,T beta,T* y) {
 if(!A || op!=SPARSE_OPERATION_NON_TRANSPOSE || A->complex!=!std::is_same_v<T,double>) return SPARSE_STATUS_INVALID_VALUE;
 if(d.diag!=SPARSE_DIAG_NON_UNIT || d.mode!=SPARSE_FILL_MODE_UPPER) return SPARSE_STATUS_NOT_SUPPORTED;
 if(d.type!=SPARSE_MATRIX_TYPE_GENERAL && A->rows!=A->cols) return SPARSE_STATUS_NOT_SUPPORTED;
 for(int r=0;r<A->rows;r++) y[r]*=beta;
 const T* v=static_cast<const T*>(A->values);
 // Serial traversal avoids data races in triangular Hermitian accumulation.
 for(int r=0;r<A->rows;r++) for(int p=A->begin[r];p<A->end[r];p++) {
  int c=A->index[p]; if(c<0 || c>=A->cols) return SPARSE_STATUS_INVALID_VALUE;
  if(d.type!=SPARSE_MATRIX_TYPE_GENERAL && c<r) return SPARSE_STATUS_INVALID_VALUE;
  y[r]+=alpha*v[p]*x[c];
  if(d.type!=SPARSE_MATRIX_TYPE_GENERAL && c!=r)
   y[c]+=alpha*(d.type==SPARSE_MATRIX_TYPE_HERMITIAN?mmhk_conj(v[p]):v[p])*x[r];
 }
 return SPARSE_STATUS_SUCCESS;
}
inline sparse_status_t mkl_sparse_d_mv(sparse_operation_t o,double a,sparse_matrix_t A,matrix_descr d,const double* x,double b,double* y) { return mmhk_mv(o,a,A,d,x,b,y); }
inline sparse_status_t mkl_sparse_z_mv(sparse_operation_t o,std::complex<double> a,sparse_matrix_t A,matrix_descr d,const std::complex<double>* x,std::complex<double> b,std::complex<double>* y) { return mmhk_mv(o,a,A,d,x,b,y); }
struct MKLVersion { int MajorVersion=0,MinorVersion=0,UpdateVersion=0; const char* ProductStatus="OpenBLAS compatibility"; const char* Build="LP64"; const char* Platform="portable"; const char* Processor="OpenBLAS"; };
inline void MKL_Get_Version(MKLVersion* v) { *v=MKLVersion{}; }
inline int mkl_get_max_threads() { return openblas_get_num_threads(); }
// FEAST is not used by the MMHK workflow. Never silently approximate it.
template<class... A> inline void feastinit(A...) { throw std::runtime_error("FEAST is unavailable in the portable OpenBLAS backend; use Lanczos/IRAM"); }
template<class... A> inline void zfeast_hcsrev(A...) { throw std::runtime_error("FEAST is unavailable in the portable OpenBLAS backend"); }
