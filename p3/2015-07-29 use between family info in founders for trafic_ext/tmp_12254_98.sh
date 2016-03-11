#!/bin/bash
set -o pipefail
R --vanilla --args seed 1097 n_rep 500 r 1.4 n_family 1000 p_dis 0.3 risk.variant 39 < mainSim.R > mainSim_500_1.4_1000_0.3_39.Rout1097 2>&1
exit $?
