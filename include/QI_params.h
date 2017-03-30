#ifndef _QI_PARAMS_H
#define _QI_PARAMS_H
#include <quda.h>
#include <Cdefs.h>
namespace qi_qcd{

  struct QI_params {
    QudaGaugeParam gauge_param;
    QudaInvertParam mg_inv_param;
    QudaMultigridParam mg_param;
    QudaInvertParam inv_param;
    QI_params(){
      gauge_param=newQudaGaugeParam();
      mg_inv_param=newQudaInvertParam();
      mg_param=newQudaMultigridParam();
      inv_param=newQudaInvertParam();
    }
  };

  struct QI_geo {
    int xdim; // lLx local dimensions
    int ydim; // lLy
    int zdim; // lLz
    int tdim; // lLt
    int gridsize[4]; // number of processors in each direction
    int L[4]; // total Lattice size
  };

EXTERN_C
  QudaTboundary get_boundary(char* s);
 const char* get_boundary_str(QudaTboundary boundary);
  QudaPrecision get_prec(char* s);
  const char* get_prec_str(QudaPrecision prec);
  QudaReconstructType get_recon(char* s);
  const char* get_recon_str(QudaReconstructType recon);
  QudaDslashType get_dslash_type(char* s);
  const char* get_dslash_str(QudaDslashType type);
  QudaInverterType get_solver_type(char* s);
  const char* get_solver_str(QudaInverterType type);
  QudaVerbosity get_verbosity(char* s);
EXTERN_C_END
}

#endif
