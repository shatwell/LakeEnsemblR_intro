#!/bin/bash

gfortran data_parameters.f90 src_flake_derivedtypes.f90 src_flake_parameters.f90 src_flake_configure.f90 src_flake_albedo_ref.f90  src_flake_paramoptic_ref.f90 src_flake.f90 src_flake_radflux.f90 src_SfcFlx.f90 Qinflow.f90 src_flake_interface_1D.f90 FLake1D.f90 -o nixflake
