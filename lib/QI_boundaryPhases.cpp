#include <QI_boundaryPhases.h>
#include <QI_params.h>
#include <util_quda.h>
#include <complex>
#include <cassert>
#include <comm_quda.h>
#include <face_quda.h>
extern qi_qcd::QI_geo qi_geo;
extern qi_qcd::QI_params qi_params;
typedef std::complex<double> Complex;
extern Topology *default_topo;
EXTERN_C
// warning this functions assumes that the data are Lexicographic order and in Quda order
void boundary_phases(void *vec, bool dagger){
  if(qi_params.gauge_param.t_boundary == QUDA_PERIODIC_T)return; // theta=0
  Complex phase;
  double arg;
  double pi=4*atan(1);
  for(int tt=0; tt<qi_geo.tdim; tt++){
    arg=(pi*(tt + comm_coords(default_topo)[3]*qi_geo.tdim))/qi_geo.L[3];
    if(dagger) phase = Complex(cos(arg),sin(arg));
    else phase = Complex(cos(arg),-sin(arg));
    for(int zz=0; zz<qi_geo.zdim; zz++)
      for(int yy=0; yy<qi_geo.ydim; yy++)
        for(int xx=0; xx<qi_geo.xdim; xx++)
          for(int dd=0; dd<4; dd++)
            for(int cc=0; cc<3; cc++){
              long int pt = (((tt*qi_geo.zdim + zz)*qi_geo.ydim + yy)*qi_geo.xdim + xx)*12 + dd*3 + cc;
	      Complex *elem = (Complex*)vec + pt;
	      *elem *= phase;
	    }
  }

}

EXTERN_C_END
