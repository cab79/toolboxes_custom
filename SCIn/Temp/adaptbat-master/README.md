# Adaptive Beat Alignment Test (BAT)

## Summary of adaptive BAT procedure
Subjects are presented musical stimuli superimposed by a metronome. The metronome is either aligned with the beat or perturbed in phase or tempo with varying magnitude from trial to trial. Subjects respond by pressing the keyboard button Q if off-beat or P if on-beat.

The tempo and phase tests are presented in random order. In each test, 80% of trials are test trials, and 20% of trials are catch trials (i.e. no deviation). The task always starts with a catch trial.

Trial-to-trial deviation values and final thresholds are calculated using a Bayesian framework, specifically the Zippy Estimation by Sequential Testing (ZEST) threshold procedure. Default ZEST parameters should work fine for implementing the adaptive BAT test. If you wish to change these parameters, you can do so in bat_params.m

Stimulus audio is generated dynamically for each trial, with deviations governed by the subject’s performance on preceding trials and the current state of threshold estimation. Trials are presented until the threshold algorithm converges, at which time the subject’s data for that task are output to a CSV file and the subject progresses to the next task version.

## Requirements
- MATLAB
- [Psychtoolbox-3](http://psychtoolbox.org/)

## Usage
- Before running the experiment on your machine for the first time, you must open `bat_params.m` and change the `stim_fpath`, `data_fpath`, and `ITI_fpath` to the appropriate file paths. Most users can then run the experiment with all other parameters set at default values. If, after testing, you wish to modify some aspects of the experiment, check the `bat_params.m` file first before changing code in the other functions, as the change you wish to make is likely controllable from the params file.
- The experiment is executed with the command `BAT('subject_id_string')`, where `'subject_id_string'` is a string of the name or ID of the subject. `BAT.m` is a wrapper function and is the entry point for all functions used in the adaptive BAT.

### References 
For more details on BAT:
- Iversen, JR, Patel, AD. (2008) The Beat Alignment Test (BAT): Surveying beat processing abilities in the general population. Proceedings of the 10th International Conference on Music Perception and Cognition. 465-468.

For more details on ZEST parameters and algorithm:
- King-Smith, P. E., Grigsby, S. S., Vingrys, A. J., Benes, S. C., & Supowit, A. (1994). Efficient and unbiased modifications of the QUEST threshold method: Theory, simulations, experimental evaluation and practical implementation. Vision Research, 34(7), 885–912.
- Marvit, et al. (2003). A comparison of psychophysical procedures
for level-discrimination thresholds. Journal of the Acoustical Society of America. 113(6), 3348-3361.
