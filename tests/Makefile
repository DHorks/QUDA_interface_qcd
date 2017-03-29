include ../make.inc
QUDA_INTFC=../lib/libquda_intfc.a

TESTS= invert_test.exe
INC += -I../include -I../lib -I./

all:$(TESTS)

test_get_params.exe: test_get_params.o $(QUDA_INTFC)
	$(CC) $^ -o $@ $(LDFLAGS)

invert_test.exe: invert_test.o $(QUDA_INTFC)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c $(HDRS)
	$(CC) $(INC) $(CFLAGS) $< -c -o $@ $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) $(INC) $(CXXFLAGS) $< -c -o $@ $(LDFLAGS)

clean:
	-rm -f *.o $(TESTS)

.PHONY: all clean
