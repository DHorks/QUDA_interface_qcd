############ make.inc included by the Makefile

MY_HOME = /ccs/home/sbacchio
LIME_HOME = $(MY_HOME)
QMP_HOME = $(MY_HOME)
QIO_HOME = $(MY_HOME)
QUDA_INTERFACE_HOME = $(MY_HOME)/src/QUDA_interface_qcd
QUDA_INSTALL_PATH = $(MY_HOME)/src/quda/build
CUDA_INSTALL_PATH = /opt/nvidia/cudatoolkit7.5/7.5.18-1.0502.10743.2.1
MAGMA_INSTALL_PATH = /sw/xk6/magma/1.6.2/cle5.2up04_gnu4.9.0
MPI_HOME = /opt/cray/mpt/7.5.2/gni/mpich-gnu/4.9

FECC = cc
FECXX = CC

CC  = $(FECC)
CXX = $(FECXX)

INC = -I$(MY_HOME)/include
INC += -I$(QUDA_INTERFACE_HOME)/include
INC += -I$(MPI_HOME)/include64
INC += -I$(QUDA_INSTALL_PATH)/include
INC += -I$(QUDA_INSTALL_PATH)/lib
INC += -I$(LIME_HOME)/include
# Comment in case QMP not used in QUDA
INC += -I$(QMP_HOME)/include -DQMP_COMMS -DMULTI_GPU

LIB = -L$(QUDA_INSTALL_PATH)/lib -lquda
LIB += -L$(MAGMA_INSTALL_PATH)/lib -lmagma -lcublas -L/opt/acml/5.3.1/gfortran64/lib -lacml
LIB += -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcufft
# Comment in case QIO not used in QUDA
LIB += -L$(QIO_HOME)/lib -lqio
# Comment in case QMP not used in QUDA
LIB += -L$(QMP_HOME)/lib -lqmp
LIB += -L$(LIME_HOME)/lib -llime

CFLAGS = -Wall -O3 -std=c99 $(INC)
CXXFLAGS = -Wall -O3 -std=c++11 $(INC)
LDFLAGS = -fPIC $(LIB)
