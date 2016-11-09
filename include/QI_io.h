#ifndef _QI_IO_H
#define _QI_IO_H
#include <quda.h>
#include <Cdefs.h>

EXTERN_C
void readLimeGauge(void **gauge, char *fname, QudaGaugeParam *param, QudaInvertParam *inv_param, int gridSize[4]);
void readLimeGaugeSmeared(void **gauge, char *fname, QudaGaugeParam *param, QudaInvertParam *inv_param, int gridSize[4]);
void readGaugeAndreas(void **gauge, char *fname, QudaGaugeParam *param, QudaInvertParam *inv_param, int gridSize[4]);
void applyBoundaryCondition(void **gauge, int Vh ,QudaGaugeParam *gauge_param);
int printVector(char file_pref[], void *vec);
void setSpinorSiteSize(int n);
void setDims(int *X);
EXTERN_C_END
#endif
