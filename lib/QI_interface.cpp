#include <iostream>
#include <quda.h>
#include <QI_params.h>
#include <QI_qcd.h>
#include <QI_io.h>
#include <QI_mapping_parity.h>
#include <comm_quda.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <cassert>
#define gaugeSiteSize 18
#define spinorSiteSize 24
using namespace qi_qcd;

extern QI_params qi_params;
extern QI_geo qi_geo;
extern int V;
void *gauge[4];
void *mg_preconditioner;
EXTERN_C

static void initMG(){
  mg_preconditioner = newMultigridQuda(&qi_params.mg_param);
  qi_params.inv_param.preconditioner = mg_preconditioner;
}

static void closeMG(){
  destroyMultigridQuda(mg_preconditioner);
}

void initQI(char* params, int params_len){
  initCommsQuda(params,params_len);
  getArgs_QI_qcd(params,params_len);
  char gauge_name[1024];
  strcpy(gauge_name,getParam("<cfg_name>",params,params_len));
  for (int dir = 0; dir < 4; dir++)
    gauge[dir] = malloc(V*gaugeSiteSize*sizeof(double));
  readLimeGauge(gauge,gauge_name,&qi_params.gauge_param,&qi_params.inv_param,qi_geo.gridsize);
  applyBoundaryCondition(gauge, V/2 ,&qi_params.gauge_param);
  initQuda(-1); // the value -1 is for multi gpu code
  loadGaugeQuda((void*)gauge, &qi_params.gauge_param);

  if(qi_params.mg_param.smoother_solve_type[0] == QUDA_DIRECT_PC_SOLVE || qi_params.inv_param.solve_type == QUDA_DIRECT_PC_SOLVE ){ // in case we use MG
    QudaSolveType tmpSl = qi_params.inv_param.solve_type;
    qi_params.inv_param.solve_type == QUDA_DIRECT_PC_SOLVE;
    if (qi_params.inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH || qi_params.inv_param.dslash_type == QUDA_TWISTED_CLOVER_DSLASH)
      loadCloverQuda(NULL, NULL, &qi_params.inv_param);
    qi_params.inv_param.solve_type = tmpSl;
  }
  else{ // in case we use CG
    if (qi_params.inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH || qi_params.inv_param.dslash_type == QUDA_TWISTED_CLOVER_DSLASH)
      loadCloverQuda(NULL, NULL, &qi_params.inv_param);    
  }

  if(qi_params.inv_param.inv_type_precondition == QUDA_MG_INVERTER) initMG();
}


void closeQI(){
  if(qi_params.inv_param.inv_type_precondition == QUDA_MG_INVERTER) closeMG();
  freeGaugeQuda();
  for (int dir = 0; dir < 4; dir++)
    free(gauge[dir]);
  if (qi_params.inv_param.dslash_type == QUDA_CLOVER_WILSON_DSLASH || qi_params.inv_param.dslash_type == QUDA_TWISTED_CLOVER_DSLASH) freeCloverQuda();
  endQuda();
}

void checkInvert(){
  void *spinorIn = malloc(V*spinorSiteSize*sizeof(double));
  void *spinorOut = malloc(V*spinorSiteSize*sizeof(double));
  memset(spinorOut,0,V*spinorSiteSize*sizeof(double));
  for(int i = 0 ; i < V*spinorSiteSize ; i++)
    if(i==0 && comm_rank() == 0)
      ((double*)spinorIn)[i]=1.;
    else
    ((double*)spinorIn)[i]=0.;
  qi_params.inv_param.twist_flavor=QUDA_TWIST_MINUS;
  invertQuda(spinorOut,spinorIn,&qi_params.inv_param);
  mapEvenOddToNormal(spinorOut,qi_params.inv_param,qi_geo.xdim,qi_geo.ydim,qi_geo.zdim,qi_geo.tdim);
  printVector("/home/khadjiyiannakou/QUDA_interface_qcd/spinorOut",spinorOut);
  free(spinorIn);
  free(spinorOut);
}

EXTERN_C_END
