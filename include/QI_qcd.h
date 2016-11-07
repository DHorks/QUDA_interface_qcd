#ifndef _QI_QCD_H
#define _QI_QCD_H

#include <stdio.h>
#include <Cdefs.h>
#include <QI_comms.h>

EXTERN_C

void getGridInfo(char* params, int params_len);
void getArgs_QI_qcd(char* params, int params_len);
char* getParams(char* fname,int *len);
EXTERN_C_END

#endif
