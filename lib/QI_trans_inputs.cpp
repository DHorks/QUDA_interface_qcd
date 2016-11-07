#include <string.h>
#include <stdlib.h>
#include <QI_params.h>

EXTERN_C

QudaPrecision get_prec(char* s)
{
  QudaPrecision ret = QUDA_DOUBLE_PRECISION;
  if (strcmp(s, "double") == 0){
    ret = QUDA_DOUBLE_PRECISION;
  }else if (strcmp(s, "single") == 0){
    ret = QUDA_SINGLE_PRECISION;
  }else if (strcmp(s, "half") == 0){
    ret = QUDA_HALF_PRECISION;
  }else{
    fprintf(stderr, "Error: invalid precision type\n");
    exit(1);
  }
  return ret;
}

const char* get_prec_str(QudaPrecision prec)
{
  const char* ret;
  switch( prec){
  case QUDA_DOUBLE_PRECISION:
    ret=  "double";
    break;
  case QUDA_SINGLE_PRECISION:
    ret= "single";
    break;
  case QUDA_HALF_PRECISION:
    ret= "half";
    break;
  default:
    ret = "unknown";
    break;
  }
  return ret;
}


QudaReconstructType get_recon(char* s)
{
  QudaReconstructType  ret;
  if (strcmp(s, "8") == 0){
    ret =  QUDA_RECONSTRUCT_8;
  }else if (strcmp(s, "9") == 0){
    ret =  QUDA_RECONSTRUCT_9;
  }else if (strcmp(s, "12") == 0){
    ret =  QUDA_RECONSTRUCT_12;
  }else if (strcmp(s, "13") == 0){
    ret =  QUDA_RECONSTRUCT_13;
  }else if (strcmp(s, "18") == 0){
    ret =  QUDA_RECONSTRUCT_NO;
  }else{
    fprintf(stderr, "Error: invalid reconstruct type\n");
    exit(1);
  }
  return ret;
}

const char* get_recon_str(QudaReconstructType recon)
{
  const char* ret;
  switch(recon){
  case QUDA_RECONSTRUCT_13:
    ret="13";
    break;
  case QUDA_RECONSTRUCT_12:
    ret= "12";
    break;
  case QUDA_RECONSTRUCT_9:
    ret="9";
    break;
  case QUDA_RECONSTRUCT_8:
    ret = "8";
    break;
  case QUDA_RECONSTRUCT_NO:
    ret = "18";
    break;
  default:
    ret="unknown";
    break;
  }
  return ret;
}


 QudaDslashType get_dslash_type(char* s)
{
  QudaDslashType ret =  QUDA_INVALID_DSLASH;
  if (strcmp(s, "wilson") == 0){
    ret = QUDA_WILSON_DSLASH;
  }else if (strcmp(s, "clover") == 0){
    ret = QUDA_CLOVER_WILSON_DSLASH;
  }else if (strcmp(s, "twisted-mass") == 0){
    ret = QUDA_TWISTED_MASS_DSLASH;
  }else if (strcmp(s, "twisted-clover") == 0){
    ret = QUDA_TWISTED_CLOVER_DSLASH;
  }else{
    fprintf(stderr, "Error: invalid dslash type\n");
    exit(1);
  }
  return ret;
}

const char* get_dslash_str(QudaDslashType type)
{
  const char* ret;
  switch( type){
  case QUDA_WILSON_DSLASH:
    ret=  "wilson";
    break;
  case QUDA_CLOVER_WILSON_DSLASH:
    ret= "clover";
    break;
  case QUDA_TWISTED_MASS_DSLASH:
    ret= "twisted-mass";
    break;
  case QUDA_TWISTED_CLOVER_DSLASH:
    ret= "twisted-clover";
    break;
  default:
    ret = "unknown";
    break;
  }
  return ret;
}


 QudaInverterType get_solver_type(char* s)
{
  QudaInverterType ret =  QUDA_INVALID_INVERTER;
  if (strcmp(s, "cg") == 0){
    ret = QUDA_CG_INVERTER;
  } else if (strcmp(s, "bicgstab") == 0){
    ret = QUDA_BICGSTAB_INVERTER;
  } else if (strcmp(s, "gcr") == 0){
    ret = QUDA_GCR_INVERTER;
  } else if (strcmp(s, "pcg") == 0){
    ret = QUDA_PCG_INVERTER;
  } else if (strcmp(s, "mpcg") == 0){
    ret = QUDA_MPCG_INVERTER;
  } else if (strcmp(s, "mpbicgstab") == 0){
    ret = QUDA_MPBICGSTAB_INVERTER;
  } else if (strcmp(s, "mr") == 0){
    ret = QUDA_MR_INVERTER;
  } else if (strcmp(s, "sd") == 0){
    ret = QUDA_SD_INVERTER;
  } else if (strcmp(s, "eigcg") == 0){
    ret = QUDA_EIGCG_INVERTER;
  } else if (strcmp(s, "inc-eigcg") == 0){
    ret = QUDA_INC_EIGCG_INVERTER;
  } else if (strcmp(s, "gmresdr") == 0){
    ret = QUDA_GMRESDR_INVERTER;
  } else if (strcmp(s, "gmresdr-proj") == 0){
    ret = QUDA_GMRESDR_PROJ_INVERTER;
  } else if (strcmp(s, "gmresdr-sh") == 0){
    ret = QUDA_GMRESDR_SH_INVERTER;
  } else if (strcmp(s, "fgmresdr") == 0){
    ret = QUDA_FGMRESDR_INVERTER;
  } else if (strcmp(s, "mg") == 0){
    ret = QUDA_MG_INVERTER;
  } else {
    fprintf(stderr, "Error: invalid solver type\n");
    exit(1);
  }
  return ret;
}

const char* get_solver_str(QudaInverterType type)
{
  const char* ret;
  switch(type){
  case QUDA_CG_INVERTER:
    ret = "cg";
    break;
  case QUDA_BICGSTAB_INVERTER:
    ret = "bicgstab";
    break;
  case QUDA_GCR_INVERTER:
    ret = "gcr";
    break;
  case QUDA_PCG_INVERTER:
    ret = "pcg";
    break;
  case QUDA_MPCG_INVERTER:
    ret = "mpcg";
    break;
  case QUDA_MPBICGSTAB_INVERTER:
    ret = "mpbicgstab";
    break;
  case QUDA_MR_INVERTER:
    ret = "mr";
    break;
  case QUDA_SD_INVERTER:
    ret = "sd";
    break;
  case QUDA_EIGCG_INVERTER:
    ret = "eigcg";
    break;
  case QUDA_INC_EIGCG_INVERTER:
    ret = "inc-eigcg";
    break;
  case QUDA_GMRESDR_INVERTER:
    ret = "gmresdr";
    break;
  case QUDA_GMRESDR_PROJ_INVERTER:
    ret = "gmresdr-proj";
    break;
  case QUDA_GMRESDR_SH_INVERTER:
    ret = "gmresdr-sh";
    break;
  case QUDA_FGMRESDR_INVERTER:
    ret = "fgmresdr";
    break;
  case QUDA_MG_INVERTER:
    ret= "mg";
    break;
  default:
    ret = "unknown";
    fprintf(stderr, "Error: invalid solver type %d\n",type);
    exit(-1);
    break;
  }
  return ret;
}

EXTERN_C_END
