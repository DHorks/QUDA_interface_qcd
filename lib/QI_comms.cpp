#include <Cdefs.h>
#include <QI_params.h> 
#include <QI_comms.h>
#include <QI_qcd.h>
#include <comm_quda.h>
#include <cassert>
#include <stdlib.h>

extern qi_qcd::QI_geo qi_geo;
extern Topology *default_topo;

// default topology for Quda (t,z,y,x) where t is varying faster
/*
static int lex_rank_from_coords(const int *coords, void *fdata)
{
  int *md = (int*) (fdata);
  int rank = coords[0];
  for (int i = 1; i < 4; i++) {
    rank = md[i] * rank + coords[i];
  }
  return rank;
}
*/

// topology which is the same as qcd package which is (t,x,y,z) with t varying faster
static int lex_rank_from_coords_txyz(const int *coords, void *fdata){
  int *md = (int*) (fdata);
  return coords[3] + coords[0]*md[3] + coords[1]*md[0]*md[3] + coords[2]*md[1]*md[0]*md[3];
}


void initCommsQuda(char* params, int params_len){
  getGridInfo(params,params_len); // initialize qi_geo
  for(int i = 0 ; i < 4 ; i++)
    assert(qi_geo.gridsize[i] > 0);
  assert(qi_geo.xdim>0);assert(qi_geo.ydim>0);
  assert(qi_geo.zdim>0);assert(qi_geo.tdim>0);
  QudaCommsMap func=lex_rank_from_coords_txyz;
  initCommsGridQuda(4,qi_geo.gridsize,func,qi_geo.gridsize);
}

