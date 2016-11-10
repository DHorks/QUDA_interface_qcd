#ifndef _QI_BOUNDARYPHASES_H
#define _QI_BOUNDARYPHASES_H
#include <quda.h>
#include <Cdefs.h>
EXTERN_C

/* this function multiplies a vector with u=exp(-i\pi \theta t/Lt)
 For now is enough to twist the fields only in the time direction
 For theta=0 we will have periodic boundary conditions and when
 theta=1 we will have antiperiodic
*/
void boundary_phases(void *vec, bool dagger);

EXTERN_C_END
#endif
