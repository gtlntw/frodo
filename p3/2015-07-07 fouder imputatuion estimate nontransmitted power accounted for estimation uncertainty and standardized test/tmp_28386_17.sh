#!/bin/bash
set -o pipefail
R --vanilla --args seed 1016 n_rep 250 r 1.2 n_family 1000 p_dis 0.1 risk.variant 2 < mainSim.R > mainSim_250_1.2_1000_0.1_2.Rout1016 2>&1
exit $?
