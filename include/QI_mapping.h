#ifndef _QI_MAPPING_PARITY_H
#define _QI_MAPPING_PARITY_H
#include <quda.h>
#include <Cdefs.h>

EXTERN_C
void mapNormalToEvenOdd(void *spinor, QudaInvertParam param, int nx , int ny , int nz, int nt);
void mapEvenOddToNormal(void *spinor, QudaInvertParam param, int nx , int ny , int nz, int nt);
void mapNormalToEvenOddGauge(void **gauge, QudaGaugeParam param, int nx , int ny , int nz, int nt);
void mapEvenOddToNormalGauge(void **gauge, QudaGaugeParam param, int nx , int ny , int nz, int nt);
void mapGauge_qcd_to_QUDA_eo(const void *gauge_qcd, void **gauge_Quda);
void mapVector_qcd_to_QUDA(const void *vec_qcd, void *vec_QUDA);
void mapVector_QUDA_to_qcd(void *vec_QUDA, void *vec_qcd);
EXTERN_C_END
#endif
