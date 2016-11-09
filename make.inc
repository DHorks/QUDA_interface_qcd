############ make.inc included by the Makefile

CUDA_INSTALL_PATH = /home/buildsets/eb141118/software/CUDA/7.0.28
QUDA_INSTALL_PATH = /home/khadjiyiannakou/quda-develop-multi-grid_CLEAN
MAGMA_INSTALL_PATH = /home/khadjiyiannakou/magma-1.7.0
MPI_HOME = /home/buildsets/eb141118/software/impi/5.0.1.035-iccifort-2015.0.090/
MY_HOME = /home/khadjiyiannakou/QUDA_interface_qcd
LIME_HOME = /home/khadjiyiannakou/lime-1.3.2_intel

FECC = mpiicc
FECXX = mpiicpc

CC  = $(FECC)
CXX = $(FECXX)

INC = -I$(MY_HOME)/include
INC += -I${MPI_HOME}/include64
INC += -I$(QUDA_INSTALL_PATH)/include
INC += -I$(QUDA_INSTALL_PATH)/lib
INC += -I$(LIME_HOME)/include

LIB = -L$(QUDA_INSTALL_PATH)/lib -lquda
LIB += -L$(MAGMA_INSTALL_PATH)/lib -lmagma
LIB += -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcufft
LIB += -L$(LIME_HOME)/lib -llime

CFLAGS = -Wall -O3 -std=c99 $(INC)
CXXFLAGS = -Wall -O3 $(INC)
LDFLAGS = -fPIC $(LIB)