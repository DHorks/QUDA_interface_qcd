#ifndef _QI_QCD_H
#define _QI_QCD_H

#include <stdio.h>
#include <Cdefs.h>
#include <QI_comms.h>
#include <QI_io.h>
EXTERN_C

void getGridInfo(char* params, int params_len);
void getArgs_QI_qcd(char* params, int params_len);
char* getParams(char* fname,int *len);
char* getParam(const char token[],char* params,int len);

void initQI(char* params, int params_len);
void closeQI();
void checkInvert_Up();
void checkInvert_Down();
void initQI_qcd(void *gauge_qcd,char* params, int params_len);
void closeQI_qcd();
void invert_QI_qcd_up(const void *spinorIn, void *spinorOut);
void invert_QI_qcd_down(const void *spinorIn, void *spinorOut);
EXTERN_C_END

#endif
