#!/bin/bash
set -o pipefail
R --vanilla --args seed 1004 n_rep 300 r 1 n_family 2000 p_dis 0.3 risk.variant 2 < mainSim.R > mainSim_300_1_2000_0.3_2.Rout1004 2>&1
exit $?
