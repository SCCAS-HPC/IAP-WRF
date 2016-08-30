#!/bin/bash

#gfortran -convert big_endian precision.F90 paramodel.F90 restart.F90
#gfortran -fconvert=big-endian precision.F90 paramodel.F90 restart_digger.F90
 gfortran -fconvert=big-endian precision.F90 paramodel.F90 snow_restart.F90
