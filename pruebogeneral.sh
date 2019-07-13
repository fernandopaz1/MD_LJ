#!/bin/bash
rm md.e
reset
make clean
make all
./md.e
ipython3 -i problema1.py
#ipython3 -i problema2.py
#ipython3 -i problema3.py
