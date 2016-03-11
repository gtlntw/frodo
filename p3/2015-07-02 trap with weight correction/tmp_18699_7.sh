#!/bin/bash
set -o pipefail
R --vanilla --args seed 1006 n_rep 500 r 1.1 n_family 1000 p_dis 0.3 risk.variant 2 < mainSim.R > mainSim_500_1.1_1000_0.3_2.Rout1006 2>&1
exit $?
