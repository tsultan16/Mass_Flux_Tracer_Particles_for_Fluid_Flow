#!/bin/sh
ifort xorshift32.f90 -o a1
ifort taus88.f90 -o a2
ifort genericrand.f90 -o a3

./a3
./a2
./a1
