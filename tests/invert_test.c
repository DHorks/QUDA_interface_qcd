#include <mpi.h>
#include <QI_qcd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc,char* argv[]){
  int numprocs;
  int myid;
  char*  params = NULL;
#ifdef QMP_COMMS
  QMP_thread_level_t tl;
  QMP_init_msg_passing(&argc, &argv, QMP_THREAD_MULTIPLE, &tl);
#else
  MPI_Init(&argc, &argv);
#endif
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);         // num. of processes taking part in the calculation                                                                                        
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);             // each process gets its ID                                         
  char   param_name[1024];
  int params_len;
  if(argc!=2)
    {
      if(myid==0) fprintf(stderr,"No input file specified\n");
      exit(EXIT_FAILURE);
    }

  strcpy(param_name,argv[1]);
  int ii=0;
  if(myid==0)
    {
      ii=0;
      printf("opening input file %s\n",param_name);
      params=getParams(param_name,&params_len);
      if(params==NULL)
	{
	  ii=1;
	}
    }
  MPI_Bcast(&ii,1,MPI_INT, 0, MPI_COMM_WORLD);
  if(ii==1) exit(EXIT_FAILURE);
  MPI_Bcast(&params_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(myid!=0) params = (char*) malloc(params_len*sizeof(char));
  MPI_Bcast(params, params_len, MPI_CHAR, 0, MPI_COMM_WORLD);


  initQI(params, params_len);
  checkInvert_Up();
  closeQI();


  /* initCommsQuda(params,params_len); */
  /* getArgs_QI_qcd(params,params_len); */

  /* char gauge_name[257]; */
  /* strcpy(gauge_name,getParam("<cfg_name>",params,params_len)); */
  /* printf("%s\n",gauge_name); */

#ifdef QMP_COMMS
  QMP_finalize_msg_passing();
#else
  MPI_Finalize();
#endif
  return 0;
}
