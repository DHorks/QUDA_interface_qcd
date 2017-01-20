#!/bin/bash
export QUDA_RESOURCE_PATH=/users/krikitos/build/build_quda-mg_latest
srun -n 1 ./invert_test.exe test_mg.ini
