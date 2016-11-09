#ifndef _QI_MAPPING_PARITY_H
#define _QI_MAPPING_PARITY_H
#include <quda.h>
#include <Cdefs.h>

EXTERN_C
void mapNormalToEvenOdd(void *spinor, QudaInvertParam param, int nx , int ny , int nz, int nt);
void mapEvenOddToNormal(void *spinor, QudaInvertParam param, int nx , int ny , int nz, int nt);
void mapNormalToEvenOddGauge(void **gauge, QudaGaugeParam param, int nx , int ny , int nz, int nt);
void mapEvenOddToNormalGauge(void **gauge, QudaGaugeParam param, int nx , int ny , int nz, int nt);
EXTERN_C_END
#endif
