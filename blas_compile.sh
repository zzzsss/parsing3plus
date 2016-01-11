#!/usr/bin/bash

g++ -g -Wall  */*.cpp main.cpp -DBLAS_ATLAS -lblas -lboost_regex -o nnpgdp
