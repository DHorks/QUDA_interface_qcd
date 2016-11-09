include ../make.inc
QUDA_INTFC=../lib/libquda_intfc.a

TESTS=test_get_params.exe invert_test.exe
INC += -I../include -I../lib -I./

all:$(TESTS)

test_get_params.exe: test_get_params.o $(QUDA_INTFC)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDFLAGS)

invert_test.exe: invert_test.o $(QUDA_INTFC)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $< -c -o $@

%.o: %.cpp $(HDRS)
	$(CXX) $(INC) $(CXXFLAGS) $< -c -o $@

clean:
	-rm -f *.o $(TESTS)

.PHONY: all clean