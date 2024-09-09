#!/bin/bash

pth=$(pwd)

## run Matlab step 1
matlab_run eegstats_condor_resubmit.m

## submit job to condor pool
matlab_submit eegstats_encode_LMM_run.sub

## run Matlab step 2
matlab_run eegstats_condor_monitor_resubmit.m
matlab_run eegstats_condor_move_resubmit.m
matlab_run eegstats_encode_compile_samples.m