#!/bin/bash

PROG="../test/test_stokes_poiseuille_flow"

INPUT="../test/input_files/stokes_poiseuille_flow_sensitivity_input.in"

PETSC_OPTIONS="-ksp_monitor"

${LIBMESH_RUN:-} $PROG $INPUT $PETSC_OPTIONS
