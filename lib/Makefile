include ../make.inc
QUDA_INTFC = libquda_intfc.a
QI_OBJS=
QI_HDRS=
INL_HDRS=
INC += ./

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
