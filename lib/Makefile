include ../make.inc
QUDA_INTFC = libquda_intfc.a
QI_OBJS= QI_get_params.o QI_comms.o QI_trans_inputs.o QI_io.o QI_mapping.o QI_interface.o QI_gammaTrans.o QI_boundaryPhases.o
QI_HDRS= Cdefs.h QI_comms.h QI_params.h QI_qcd.h QI_io.h QI_mapping.h QI_gammaTrans.h QI_boundaryPhases.h
INL_HDRS=
INC += -I./

HDRS = $(QI_HDRS:%=../include/%)
HDRS += $(INL_HDRS)

all: $(QUDA_INTFC)

$(QUDA_INTFC): $(QI_OBJS) ../make.inc
	ar cru $@ $(QI_OBJS)

clean:
	-rm -f *.o $(QUDA_INTFC)

%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $< -c -o $@

.PHONY: all clean
