#!/bin/bash

PROG="../test/test_stokes_poiseuille_flow"

INPUT="../test/input_files/stokes_poiseuille_flow_sensitivity_input.in"

${LIBMESH_RUN:-} $PROG $INPUT 