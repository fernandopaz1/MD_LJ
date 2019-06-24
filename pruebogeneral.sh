#!/bin/bash
reset
make clean
make all
./md.e
ipython3 -i pruebogeneral.py
