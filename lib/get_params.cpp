#include <iostream>
#include <QI_params.h>
#include <QI_qcd.h>
#include <comm_quda.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <cassert>

using namespace qi_qcd;

QI_params qi_params; // global variable to store the parameters for inversions, gauge and multigrid
QI_geo qi_geo; 
#define QI_PRINT_PARAMS
EXTERN_C

static std::string to_string( int x ) {
  int length = snprintf( NULL, 0, "%d", x );
  assert( length >= 0 );
  char* buf = new char[length + 1];
  snprintf( buf, length + 1, "%d", x );
  std::string str( buf );
  delete[] buf;
  return str;
}

static char* getParam(char token[],char* params,int len)
{
  bool check_match = false;
  int i,token_len=strlen(token);
  for(i=0;i<len-token_len;i++)
    {
      if(memcmp(token,params+i,token_len)==0)
	{
	  i+=token_len;
	  *(strchr(params+i,'<'))='\0';
	  check_match=true;
	  break;
	}
      check_match=false;
    }
  if(!check_match){
    fprintf(stderr,"Error: I am looking for token [%s] but I didnt find it please provide it\n",token);
    exit(-1);
  }
  return params+i;
}

char* getParams(char* fname,int *len)
{
  FILE *pfile;
  char *params;
  int i;
  if ((pfile=fopen(fname,"r"))==NULL)
    {
      fprintf(stderr,"Error, cannot open %s for reading\n",fname);
      return(NULL);
    }
  i=0;
  while(!feof(pfile))
    {
      fgetc(pfile);
      i++;
    }
  *(len)=i;
  rewind(pfile);
  params = (char*)malloc(*(len)*sizeof(char));
  fread(params,sizeof(char),*len,pfile);
  fclose(pfile);
  return params;
}


void getGridInfo(char* params, int params_len){
  sscanf(getParam("<processors_txyz>",params,params_len),"%hd %hd %hd %hd",&qi_geo.gridsize[3], &qi_geo.gridsize[0], &qi_geo.gridsize[1], &qi_geo.gridsize[2]);
  sscanf(getParam("<lattice_txyz>",params,params_len),"%hd %hd %hd %hd",&qi_geo.tdim, &qi_geo.xdim, &qi_geo.ydim, &qi_geo.zdim);
}



void getArgs_QI_qcd(char* params, int params_len){
  QudaPrecision prec = get_prec( getParam("<QUDA_prec>",params,params_len) ); // Precision in GPU <double/single/half>
  QudaPrecision prec_sloppy = get_prec( getParam("<QUDA_prec_sloppy>",params,params_len) ); // Sloppy precision in GPU <double/single/half>
  QudaPrecision prec_preco = get_prec( getParam("<QUDA_prec_preco>",params,params_len) ); // Preconditioner precision in GPU <double/single/half>
  QudaReconstructType recon = get_recon( getParam("<QUDA_recon>",params,params_len) ); // Link reconstruction type <8/9/12/13/18>
  QudaReconstructType recon_sloppy = get_recon( getParam("<QUDA_recon_sloppy>",params,params_len) ); // Link reconstruction type for sloppy <8/9/12/13/18>
  QudaReconstructType recon_preco = get_recon( getParam("<QUDA_recon_preco>",params,params_len) ); // Link reconstruction type for preconditioner <8/9/12/13/18>
  QudaDslashType dslash_type = get_dslash_type(getParam("<QUDA_dslash_type>",params,params_len) ); // The dslash type we want to use <wilson/clover/twisted-mass/twisted-clover>
  double kappa; sscanf(getParam("<QUDA_kappa>",params,params_len),"%lf",&kappa);
  double mu_val=0; if(dslash_type == QUDA_TWISTED_MASS_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) sscanf(getParam("<QUDA_mu>",params,params_len),"%lf",&mu_val); // mu value for TMF
  double csw=0; if(dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH)sscanf(getParam("<QUDA_csw>",params,params_len),"%lf",&csw); // clover coeff
  QudaInverterType inv_type = get_solver_type( getParam("<QUDA_solver_type>",params,params_len) ); // inverter type 
  if( (inv_type != QUDA_CG_INVERTER) && (inv_type != QUDA_MG_INVERTER) ){ // I will focus on these inverters now
    if(comm_rank() == 0) fprintf(stderr,"Error: The Quda interface supports only CG and MG inverters up to now\n");
    exit(-1);
  }
  int mg_nlvls=0;
  std::vector<int> mg_nNulls;
  int tmp_int;
  std::string str_nullVec;
  std::vector<int> block_xyzt(4,0);
  std::vector<std::vector<int> > all_block_xyzt;
  std::string str_blocks;
  double mg_mu_coarse=0.;
  int mg_nu_pre=0;
  int mg_nu_post=0;
  if(inv_type == QUDA_MG_INVERTER){
    sscanf(getParam("<QUDA_MG_nlvls>",params,params_len),"%d",&mg_nlvls); // number of multigrid levels
    assert(mg_nlvls > 0);
    for(int i = 0 ; i < mg_nlvls; i++){
      str_nullVec = "<QUDA_MG_lvl" + to_string(i) + "_nNulls>";
      sscanf(getParam(&str_nullVec[0],params,params_len),"%d",&tmp_int); // number of null vectors for each level
      mg_nNulls.push_back(tmp_int);
    }
    for(int i = 0 ; i < mg_nlvls; i++){
      str_blocks = "<QUDA_MG_lvl" + to_string(i) + "_blck_size_xyzt>";
      sscanf(getParam(&str_blocks[0],params,params_len),"%d %d %d %d",&block_xyzt[0],&block_xyzt[1],&block_xyzt[2],&block_xyzt[3]); // block size for each level
      all_block_xyzt.push_back(block_xyzt);
    }
    sscanf(getParam("<QUDA_MG_mu_coarse>",params,params_len),"%lf",&mg_mu_coarse); // the mu value for multigrid at coarse level
    sscanf(getParam("<QUDA_MG_nu_pre>",params,params_len),"%d",&mg_nu_pre); //The number of pre-smoother applications to do at each multigrid level
    sscanf(getParam("<QUDA_MG_nu_post>",params,params_len),"%d",&mg_nu_post); //The number of post-smoother applications to do at each multigrid level
  }

#ifdef QI_PRINT_PARAMS
  if(comm_rank() == 0){
    printf("<prec> = %s\n",get_prec_str(prec));
    printf("<prec_sloppy> = %s\n",get_prec_str(prec_sloppy));
    printf("<prec_preco> = %s\n",get_prec_str(prec_preco));
    printf("<recon> = %s\n",get_recon_str(recon));
    printf("<recon_sloppy> = %s\n",get_recon_str(recon_sloppy));
    printf("<recon_preco> = %s\n",get_recon_str(recon_preco));
    printf("<dslash_type> = %s\n",get_dslash_str(dslash_type));
    printf("<inv_type> = %s\n",get_solver_str(inv_type));
    printf("<mu_val> = %f\n",mu_val);
    printf("<kappa> = %f\n",kappa);
    printf("<csw = %f>\n",csw);
    if(inv_type == QUDA_MG_INVERTER){
      printf("<mg_nlvls> = %d\n",mg_nlvls);
      for(int i = 0 ; i < mg_nlvls; i++)
	printf("<mg_lvl%d> <Nnulls=%d>\n",i,mg_nNulls[i]);
      for(int i = 0 ; i < mg_nlvls; i++)
	printf("<mg_lvl%d> <block=(x=%d,y=%d,z=%d,t=%d)>\n",i,all_block_xyzt[i][0],all_block_xyzt[i][1],all_block_xyzt[i][2],all_block_xyzt[i][3]);
      printf("<mg_mu_coarse> = %f\n",mg_mu_coarse);
      printf("<mg_nu_pre> = %d\n",mg_nu_pre);
      printf("<mg_nu_post> = %d\n",mg_nu_post);
    }
  }
#endif

}


EXTERN_C_END
