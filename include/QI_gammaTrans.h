#ifndef _QI_GAMMATRANS_H
#define _QI_GAMMATRANS_H
#include <quda.h>
#include <Cdefs.h>
EXTERN_C

/* 
this performs the transformation \gamma_4_dr \gamma_5_dr U_ch_dr^\dag
where dr = Degrand_Rossi and ch = chiral basis used in tmLQCD
 */
void gammaRotate_dg_ch(void *vec);

EXTERN_C_END
#endif
