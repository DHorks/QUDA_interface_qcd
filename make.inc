############ make.inc included by the Makefile

CUDA_INSTALL_PATH = /opt/nvidia/cudatoolkit8.0/8.0.44_GA_2.2.7_g4a6c213-2.1
QUDA_INSTALL_PATH = /users/krikitos/build/build_quda-mg_latest
MAGMA_INSTALL_PATH = /users/krikitos/build/magma-1.7.0_intel/
MPI_HOME = /opt/cray/pe/mpt/7.5.0/gni/mpich-intel/16.0/
MY_HOME = /users/krikitos/build/QUDA_interface_qcd
LIME_HOME = /users/krikitos/install/intel
ARPACK_HOME= /users/krikitos/install/intel

FECC = cc
FECXX = CC

CC  = $(FECC)
CXX = $(FECXX)

INC = -I$(MY_HOME)/include
INC += -I${MPI_HOME}/include64
INC += -I$(QUDA_INSTALL_PATH)/include
INC += -I$(QUDA_INSTALL_PATH)/lib
INC += -I$(LIME_HOME)/include

LIB = -L$(QUDA_INSTALL_PATH)/lib -lquda
LIB += -L$(MAGMA_INSTALL_PATH)/lib -lmagma -lcublas -lgsl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lstdc++ -liomp5 -lhdf5
LIB += -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcufft
LIB += -L$(LIME_HOME)/lib -llime
LIB += -L$(ARPACK_HOME)/lib -lparpack -larpack

CFLAGS = -Wall -O3 -std=c99 $(INC)
CXXFLAGS = -Wall -O3 $(INC)
LDFLAGS = -fPIC $(LIB)
