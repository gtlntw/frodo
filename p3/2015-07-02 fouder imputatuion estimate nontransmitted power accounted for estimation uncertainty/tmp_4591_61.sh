#!/bin/bash
set -o pipefail
R --vanilla --args seed 1060 n_rep 500 r 1.6 n_family 1000 p_dis 0.3 risk.variant 7 < mainSim.R > mainSim_500_1.6_1000_0.3_7.Rout1060 2>&1
exit $?
