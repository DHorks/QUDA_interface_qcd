#include <iostream>
#include <QI_params.h>
#include <QI_qcd.h>
#include <QI_io.h>
#include <comm_quda.h>
#include <face_quda.h>
#include <util_quda.h>
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

char* getParam(const char token[],char* params,int len)
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
  sscanf(getParam("<processors_txyz>",params,params_len),"%d %d %d %d",&qi_geo.gridsize[3], &qi_geo.gridsize[0], &qi_geo.gridsize[1], &qi_geo.gridsize[2]);
  sscanf(getParam("<lattice_txyz>",params,params_len),"%d %d %d %d",&qi_geo.L[3], &qi_geo.L[0], &qi_geo.L[1], &qi_geo.L[2]);
  assert(qi_geo.L[3]%qi_geo.gridsize[3] == 0); qi_geo.tdim = qi_geo.L[3] / qi_geo.gridsize[3];
  assert(qi_geo.L[0]%qi_geo.gridsize[0] == 0); qi_geo.xdim = qi_geo.L[0] / qi_geo.gridsize[0];
  assert(qi_geo.L[1]%qi_geo.gridsize[1] == 0); qi_geo.ydim = qi_geo.L[1] / qi_geo.gridsize[1];
  assert(qi_geo.L[2]%qi_geo.gridsize[2] == 0); qi_geo.zdim = qi_geo.L[2] / qi_geo.gridsize[2];
}



void getArgs_QI_qcd(char* params, int params_len){
  /*
  QudaPrecision prec = get_prec( getParam("<QUDA_prec>",params,params_len) ); // Precision in GPU <double/single/half>
  QudaPrecision prec_sloppy = get_prec( getParam("<QUDA_prec_sloppy>",params,params_len) ); // Sloppy precision in GPU <double/single/half>
  QudaPrecision prec_preco = get_prec( getParam("<QUDA_prec_preco>",params,params_len) ); // Preconditioner precision in GPU <double/single/half>
  QudaReconstructType recon = get_recon( getParam("<QUDA_recon>",params,params_len) ); // Link reconstruction type <8/9/12/13/18>
  QudaReconstructType recon_sloppy = get_recon( getParam("<QUDA_recon_sloppy>",params,params_len) ); // Link reconstruction type for sloppy <8/9/12/13/18>
  QudaReconstructType recon_preco = get_recon( getParam("<QUDA_recon_preco>",params,params_len) ); // Link reconstruction type for preconditioner <8/9/12/13/18>
  */
  QudaDslashType dslash_type = get_dslash_type(getParam("<QUDA_dslash_type>",params,params_len) ); // The dslash type we want to use <wilson/clover/twisted-mass/twisted-clover>
  QudaTboundary boundary_cond = get_boundary( getParam("<QUDA_boundary_cond>",params,params_len) ); // Boundary conditions that we want to use <periodic,antiperiodic>
  QudaVerbosity verbosity = get_verbosity( getParam("<QUDA_verbosity>",params,params_len) );
  double tol; sscanf(getParam("<QUDA_tolerance>",params,params_len),"%lf",&tol); // tolerance for the inverter
  double kappa; sscanf(getParam("<QUDA_kappa>",params,params_len),"%lf",&kappa); double mass = 0.5/kappa-4.;
  double mu_val=0; if(dslash_type == QUDA_TWISTED_MASS_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) sscanf(getParam("<QUDA_mu>",params,params_len),"%lf",&mu_val); // mu value for TMF
  double csw=0; if(dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH)sscanf(getParam("<QUDA_csw>",params,params_len),"%lf",&csw); // clover coeff
  QudaInverterType inv_type = get_solver_type( getParam("<QUDA_solver_type>",params,params_len) ); // inverter type 
  if( (inv_type != QUDA_CG_INVERTER) && (inv_type != QUDA_MG_INVERTER) ){ // I will focus on these inverters now
    if(comm_rank() == 0) fprintf(stderr,"Error: The Quda interface supports only CG and MG inverters up to now\n");
    exit(-1);
  }

  int mg_nlvls=2;
  int mg_nu_pre=0;
  int mg_nu_post=0;
  double mg_mu_factor=1;
#ifdef QMP_COMMS0
  char setup_filename[1024];
#else
  char *setup_filename = "";
#endif
  int read = 0;
  QudaInverterType smoother_type = QUDA_MR_INVERTER;
  QudaInverterType coarsest_inv_type = QUDA_GCR_INVERTER;
  if(inv_type == QUDA_MG_INVERTER){
    sscanf(getParam("<QUDA_MG_nlvls>",params,params_len),"%d",&mg_nlvls); // number of multigrid levels
    assert(mg_nlvls > 1);
    qi_params.mg_param.n_level = mg_nlvls;
    for(int i = 0 ; i < mg_nlvls-1; i++){
      // number of null vectors for each level
      std::string str_nullVec = "<QUDA_MG_lvl" + to_string(i) + "_nNulls>";
      sscanf(getParam(&str_nullVec[0],params,params_len),"%d",&qi_params.mg_param.n_vec[i]);
      // aggregation block size
      std::string str_blocks = "<QUDA_MG_lvl" + to_string(i) + "_blck_size_xyzt>";
      int * block_xyzt = qi_params.mg_param.geo_block_size[i];
      sscanf(getParam(&str_blocks[0],params,params_len),"%d %d %d %d",&block_xyzt[0],&block_xyzt[1],&block_xyzt[2],&block_xyzt[3]); // block size for each level
    }
    sscanf(getParam("<QUDA_MG_delta_mu>",params,params_len),"%lf",&mg_mu_factor);
#ifdef QMP_COMMS0
    sscanf(getParam("<QUDA_MG_setup_file>",params,params_len),"%s",setup_filename);
    sscanf(getParam("<QUDA_MG_read_setup>",params,params_len),"%d",&read);
    if(read && strcmp(setup_filename,"")==0) read = 0;
#endif
    sscanf(getParam("<QUDA_MG_nu_pre>",params,params_len),"%d",&mg_nu_pre); //The number of pre-smoother applications to do at each multigrid level
    sscanf(getParam("<QUDA_MG_nu_post>",params,params_len),"%d",&mg_nu_post); //The number of post-smoother applications to do at each multigrid level
    smoother_type = get_solver_type( getParam("<QUDA_MG_smoother_type>",params,params_len) ); // inverter type 
    coarsest_inv_type = get_solver_type( getParam("<QUDA_MG_solver_type>",params,params_len) ); // inverter type 
  }

#ifdef QI_PRINT_PARAMS
  if(comm_rank() == 0){
    printf("<dslash_type> = %s\n",get_dslash_str(dslash_type));
    printf("<boundary_conditions> = %s\n",get_boundary_str(boundary_cond));
    printf("<tolerance> = %e\n",tol);
    printf("<inv_type> = %s\n",get_solver_str(inv_type));
    printf("<mu_val> = %f\n",mu_val);
    printf("<kappa> = %f\n",kappa);
    printf("<csw = %f>\n",csw);
    if(inv_type == QUDA_MG_INVERTER){
      printf("<mg_nlvls> = %d\n",mg_nlvls);
      for(int i = 0 ; i < mg_nlvls-1; i++) {
        int * block_xyzt = qi_params.mg_param.geo_block_size[i];
        printf("<mg_lvl%d> <block=(x=%d,y=%d,z=%d,t=%d)>\n",i,block_xyzt[0],block_xyzt[1],block_xyzt[2],block_xyzt[3]);
	printf("<mg_lvl%d> <Nnulls=%d>\n",i,qi_params.mg_param.n_vec[i]);
      }
      printf("<mg_delta_mu> = %f\n",mg_mu_factor);
      printf("<mg_nu_pre> = %d\n",mg_nu_pre);
      printf("<mg_nu_post> = %d\n",mg_nu_post);
      printf("<mg_smoother_type> = %s\n",get_solver_str(smoother_type));
      printf("<mg_solver_type> = %s\n",get_solver_str(coarsest_inv_type));
    }
    fflush(stdout);
  }
#endif

  //============================== Gauge params =============================//
  qi_params.gauge_param.X[0]=qi_geo.xdim;
  qi_params.gauge_param.X[1]=qi_geo.ydim;
  qi_params.gauge_param.X[2]=qi_geo.zdim;
  qi_params.gauge_param.X[3]=qi_geo.tdim;
  setDims(qi_params.gauge_param.X);
  setSpinorSiteSize(24);
  qi_params.gauge_param.anisotropy=1.;
  qi_params.gauge_param.type = QUDA_WILSON_LINKS;
  qi_params.gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
  qi_params.gauge_param.t_boundary = boundary_cond;
  qi_params.gauge_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  qi_params.gauge_param.cuda_prec = QUDA_DOUBLE_PRECISION;
  qi_params.gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
  qi_params.gauge_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  qi_params.gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_NO;
  qi_params.gauge_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
  qi_params.gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_NO;
  qi_params.gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
  
#define MAX(a,b) ((a)>(b)?(a):(b))
  int x_face_size = qi_params.gauge_param.X[1]*qi_params.gauge_param.X[2]*qi_params.gauge_param.X[3]/2;
  int y_face_size = qi_params.gauge_param.X[0]*qi_params.gauge_param.X[2]*qi_params.gauge_param.X[3]/2;
  int z_face_size = qi_params.gauge_param.X[0]*qi_params.gauge_param.X[1]*qi_params.gauge_param.X[3]/2;
  int t_face_size = qi_params.gauge_param.X[0]*qi_params.gauge_param.X[1]*qi_params.gauge_param.X[2]/2;
  int pad_size =MAX(x_face_size, y_face_size);
  pad_size = MAX(pad_size, z_face_size);
  pad_size = MAX(pad_size, t_face_size);
  qi_params.gauge_param.ga_pad = pad_size;
#undef MAX
  //=========================== Invert params ===============================//
  qi_params.inv_param.dslash_type = dslash_type;
  qi_params.inv_param.mass = mass;
  qi_params.inv_param.mu = mu_val;
  qi_params.inv_param.twist_flavor = QUDA_TWIST_SINGLET;
  qi_params.inv_param.kappa = kappa;
  qi_params.inv_param.solution_type = QUDA_MAT_SOLUTION;
  qi_params.inv_param.inv_type = inv_type;
  qi_params.inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
  qi_params.inv_param.dagger = QUDA_DAG_NO;
  qi_params.inv_param.mass_normalization = QUDA_MASS_NORMALIZATION; // it will multiply the results with the factor 2*kappa;
  qi_params.inv_param.solver_normalization = QUDA_SOURCE_NORMALIZATION; // it will normalize the source and the rescale back to avoid underflow/overflow
  if(inv_type == QUDA_CG_INVERTER)
    qi_params.inv_param.solve_type = QUDA_NORMOP_PC_SOLVE;
  else
    qi_params.inv_param.solve_type = QUDA_DIRECT_PC_SOLVE;

  if(inv_type == QUDA_MG_INVERTER){
    qi_params.inv_param.inv_type = QUDA_GCR_INVERTER;
    qi_params.inv_param.inv_type_precondition = QUDA_MG_INVERTER;
  }
  qi_params.inv_param.pipeline = 0;
  qi_params.inv_param.Nsteps = 2;
  qi_params.inv_param.gcrNkrylov = 20;
  qi_params.inv_param.tol = tol;
  qi_params.inv_param.tol_restart = 1e-3; 
  double tol_hq = 0.;

  if(inv_type ==  QUDA_CG_INVERTER){
    qi_params.inv_param.residual_type = static_cast<QudaResidualType_s>(0);
    qi_params.inv_param.residual_type = (tol != 0) ? static_cast<QudaResidualType_s> ( qi_params.inv_param.residual_type | QUDA_L2_RELATIVE_RESIDUAL) : qi_params.inv_param.residual_type;
    qi_params.inv_param.residual_type = (tol_hq != 0) ? static_cast<QudaResidualType_s> (qi_params.inv_param.residual_type | QUDA_HEAVY_QUARK_RESIDUAL) : qi_params.inv_param.residual_type;
    qi_params.inv_param.tol_hq = 0.; // specify a tolerance for the residual for heavy quark residual
    qi_params.inv_param.maxiter = 50000;
  }
  else{
    // require both L2 relative and heavy quark residual to determine convergence
    qi_params.inv_param.residual_type = static_cast<QudaResidualType>(QUDA_L2_RELATIVE_RESIDUAL);
    qi_params.inv_param.tol_hq = tol_hq; // specify a tolerance for the residual for heavy quark residual           
    qi_params.inv_param.maxiter = 300;
  }

  for (int i=0; i<qi_params.inv_param.num_offset; i++) {
    qi_params.inv_param.tol_offset[i] = qi_params.inv_param.tol;
    qi_params.inv_param.tol_hq_offset[i] = qi_params.inv_param.tol_hq;
  }
  qi_params.inv_param.reliable_delta = 1e-4; // reliable delta is related to the mixed precision solver

  // domain decomposition preconditioner parameters
  qi_params.inv_param.schwarz_type = QUDA_ADDITIVE_SCHWARZ;
  qi_params.inv_param.precondition_cycle = 1;
  qi_params.inv_param.tol_precondition = 1e-1;
  qi_params.inv_param.maxiter_precondition = 1;
  qi_params.inv_param.omega = 1.0;

  qi_params.inv_param.use_sloppy_partial_accumulator = 1;
  qi_params.inv_param.max_res_increase = 1;
  qi_params.inv_param.verbosity_precondition = QUDA_SILENT;
  qi_params.inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;

  qi_params.inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
  qi_params.inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
  qi_params.inv_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
  qi_params.inv_param.preserve_source = QUDA_PRESERVE_SOURCE_YES;
  qi_params.inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;//QUDA_UKQCD_GAMMA_BASIS;
  qi_params.inv_param.dirac_order = QUDA_DIRAC_ORDER;

  qi_params.inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
  qi_params.inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
  qi_params.inv_param.tune = QUDA_TUNE_YES; 
  qi_params.inv_param.sp_pad = 0;
  qi_params.inv_param.cl_pad = 0;
  
  if (dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
    qi_params.inv_param.clover_cpu_prec = QUDA_DOUBLE_PRECISION;
    qi_params.inv_param.clover_cuda_prec = QUDA_DOUBLE_PRECISION;
    qi_params.inv_param.clover_cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    qi_params.inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
    qi_params.inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
    qi_params.inv_param.clover_coeff = csw*kappa;
    qi_params.inv_param.compute_clover = true;
    qi_params.inv_param.return_clover = 0;
    qi_params.inv_param.compute_clover_inverse = 1;
    qi_params.inv_param.return_clover_inverse = 0;
  }
  qi_params.inv_param.Ls = 1;
  qi_params.inv_param.verbosity = verbosity;
  qi_params.inv_param.verbosity_precondition = QUDA_SILENT;
  //=================================== Multigrid params =================//
  if(inv_type == QUDA_MG_INVERTER){
    qi_params.mg_inv_param.Ls = 1;
    
    qi_params.mg_inv_param.sp_pad = 0;
    qi_params.mg_inv_param.cl_pad = 0;
    
    qi_params.mg_inv_param.cpu_prec = QUDA_DOUBLE_PRECISION;
    qi_params.mg_inv_param.cuda_prec = QUDA_DOUBLE_PRECISION;
    qi_params.mg_inv_param.cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
    qi_params.mg_inv_param.cuda_prec_precondition = QUDA_SINGLE_PRECISION;
    qi_params.mg_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
    qi_params.mg_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
    qi_params.mg_inv_param.dirac_order = QUDA_DIRAC_ORDER;
    
    if (dslash_type == QUDA_CLOVER_WILSON_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
      qi_params.mg_inv_param.clover_cpu_prec = QUDA_DOUBLE_PRECISION;
      qi_params.mg_inv_param.clover_cuda_prec = QUDA_DOUBLE_PRECISION;
      qi_params.mg_inv_param.clover_cuda_prec_sloppy = QUDA_SINGLE_PRECISION;
      qi_params.mg_inv_param.clover_cuda_prec_precondition = QUDA_SINGLE_PRECISION;
      qi_params.mg_inv_param.clover_order = QUDA_PACKED_CLOVER_ORDER;
      qi_params.mg_inv_param.clover_coeff = csw*kappa;
    }
    
    qi_params.mg_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
    qi_params.mg_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
    
    qi_params.mg_inv_param.dslash_type = dslash_type;
    
    //Free field!
    qi_params.mg_inv_param.mass = mass;
    qi_params.mg_inv_param.kappa = kappa;
    
    if (dslash_type == QUDA_TWISTED_MASS_DSLASH || dslash_type == QUDA_TWISTED_CLOVER_DSLASH) {
      qi_params.mg_inv_param.mu = mu_val;
      qi_params.mg_inv_param.twist_flavor = QUDA_TWIST_SINGLET;
    }
    
    qi_params.mg_inv_param.dagger = QUDA_DAG_NO;
    qi_params.mg_inv_param.mass_normalization = QUDA_MASS_NORMALIZATION;
    
    qi_params.mg_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
    qi_params.mg_inv_param.solution_type = QUDA_MAT_SOLUTION;
    
    qi_params.mg_inv_param.solve_type = QUDA_DIRECT_SOLVE;

  // these need to tbe set for now but are actually ignored by the MG setup
  // needed to make it pass the initialization test
    qi_params.mg_inv_param.inv_type = QUDA_GCR_INVERTER;
    qi_params.mg_inv_param.tol = 1e-10;
    qi_params.mg_inv_param.maxiter = 1000;
    qi_params.mg_inv_param.reliable_delta = 1e-10;
    qi_params.mg_inv_param.gcrNkrylov = 10;

    qi_params.mg_inv_param.verbosity = QUDA_SUMMARIZE;
    qi_params.mg_inv_param.verbosity_precondition = QUDA_SUMMARIZE;

    qi_params.mg_param.invert_param = &qi_params.mg_inv_param;
    
    for(int i = 0 ; i < mg_nlvls; i++){
      qi_params.mg_param.verbosity[i] = QUDA_SILENT;
      qi_params.mg_param.setup_inv_type[i] = QUDA_CG_INVERTER;
      qi_params.mg_param.num_setup_iter[i] = 1;
      qi_params.mg_param.setup_tol[i] = 1e-4;
      qi_params.mg_param.mu_factor[i] = 1;
      qi_params.mg_param.spin_block_size[i] = 1;
      qi_params.mg_param.nu_pre[i] = mg_nu_pre;
      qi_params.mg_param.nu_post[i] = mg_nu_post;
      qi_params.mg_param.cycle_type[i] = QUDA_MG_CYCLE_RECURSIVE;
      qi_params.mg_param.smoother[i] = smoother_type;
      qi_params.mg_param.smoother_tol[i] = 0.2; // play with this tolerance later
      qi_params.mg_param.global_reduction[i] = QUDA_BOOLEAN_YES;
      qi_params.mg_param.smoother_solve_type[i] = QUDA_DIRECT_PC_SOLVE;
      qi_params.mg_param.coarse_grid_solution_type[i] = QUDA_MATPC_SOLUTION;
      qi_params.mg_param.omega[i] = 0.85; // over/under relaxation factor
      qi_params.mg_param.location[i] = QUDA_CUDA_FIELD_LOCATION;
    }
    qi_params.mg_param.setup_tol[0] = 5e-7;
    qi_params.mg_param.setup_type = QUDA_NULL_VECTOR_SETUP;
    qi_params.mg_param.pre_orthonormalize = QUDA_BOOLEAN_NO;
    qi_params.mg_param.post_orthonormalize = QUDA_BOOLEAN_YES;
    qi_params.mg_param.spin_block_size[0] = 2;
    qi_params.mg_param.smoother[mg_nlvls-1] = coarsest_inv_type;
    qi_params.mg_param.nu_pre[mg_nlvls-1] = 300;
    qi_params.mg_param.nu_post[mg_nlvls-1] = 0;
    if(read) {
      qi_params.mg_param.compute_null_vector = QUDA_COMPUTE_NULL_VECTOR_NO;
      strcpy(qi_params.mg_param.vec_outfile,"");
      strcpy(qi_params.mg_param.vec_infile,setup_filename);
    } else {
      qi_params.mg_param.compute_null_vector = QUDA_COMPUTE_NULL_VECTOR_YES;
      strcpy(qi_params.mg_param.vec_outfile,setup_filename);
      strcpy(qi_params.mg_param.vec_infile,"");
    }
    qi_params.mg_param.generate_all_levels = QUDA_BOOLEAN_YES;
    qi_params.mg_param.run_verify = QUDA_BOOLEAN_NO;
  }

#ifdef QI_PRINT_PARAMS
  if(verbosity == QUDA_VERBOSE) {
    printfQuda("\n\n qi_params.gauge_params\n");
    printQudaGaugeParam(&qi_params.gauge_param);
    
    printfQuda("\n\n qi_params.inv_params\n");
    printQudaInvertParam(&qi_params.inv_param);
    
    if(inv_type == QUDA_MG_INVERTER){
      printfQuda("\n\n qi_params.mg_params\n");
      printQudaMultigridParam(&qi_params.mg_param);
      
      printfQuda("\n\n qi_params.mg_param.invert_params\n");
      printQudaInvertParam(qi_params.mg_param.invert_param);
    }
  }
#endif
}


EXTERN_C_END
