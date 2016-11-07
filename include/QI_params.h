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
      mg_inv_param=QudaInvertParam();
      mg_param=newQudaMultigridParam();
      inv_param=QudaInvertParam();
    }
  };

  struct QI_geo {
    int xdim;
    int ydim;
    int zdim;
    int tdim;
    int gridsize[4];
  };

EXTERN_C
  QudaPrecision get_prec(char* s);
  const char* get_prec_str(QudaPrecision prec);
  QudaReconstructType get_recon(char* s);
  const char* get_recon_str(QudaReconstructType recon);
  QudaDslashType get_dslash_type(char* s);
  const char* get_dslash_str(QudaDslashType type);
  QudaInverterType get_solver_type(char* s);
  const char* get_solver_str(QudaInverterType type);
EXTERN_C_END
}

#endif
