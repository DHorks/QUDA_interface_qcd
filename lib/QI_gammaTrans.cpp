#include <QI_gammaTrans.h>
#include <QI_params.h>
#include <util_quda.h>
#include <complex>
#include <cassert>
extern qi_qcd::QI_geo qi_geo;
extern qi_qcd::QI_params qi_params;
extern int V;

typedef std::complex<double> Complex;
EXTERN_C
void gammaRotate_dg_ch(void *vec){
  assert(qi_params.inv_param.dirac_order == QUDA_DIRAC_ORDER); //colors inside spins
  assert(qi_params.inv_param.gamma_basis == QUDA_DEGRAND_ROSSI_GAMMA_BASIS); // MG works only with DEGRAND_ROSSI
  Complex tmp[4];
  for(int i = 0 ; i < V ; i++)
    for(int c = 0 ; c < 3 ; c++){
      Complex *Cvec = ((Complex*)vec) + i*12;
      for(int s = 0 ; s < 4 ; s++)
	tmp[s] = Cvec[s*3+c];
      Cvec[0*3+c] = -tmp[3];
      Cvec[1*3+c] = tmp[2];
      Cvec[2*3+c] = tmp[1];
      Cvec[3*3+c] = -tmp[0];
    }
}
EXTERN_C_END
