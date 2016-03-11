#!/bin/bash
set -o pipefail
R --vanilla --args seed 1003 n_rep 250 r 1 n_family 1000 p_dis 0.3 risk.variant 2 < mainSim.R > mainSim_250_1_1000_0.3_2.Rout1003 2>&1
exit $?
