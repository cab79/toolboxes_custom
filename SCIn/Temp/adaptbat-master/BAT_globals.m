function globals = BAT_globals

% Get wav file names and categorize them by augmentation interval
data_loc =  '/data0/stimuli/audio/attmap/bat/stimuli_adaptive_short/';
globals.on0 = dir([data_loc,'*_B0_v2.wav']);dirlist_on0 = {globals.on0.name};
globals.offn1 = dir([data_loc,'*_B-1_v2.wav']);dirlist_offn1 = {globals.offn1.name};
globals.offn2 = dir([data_loc,'*_B-2_v2.wav']);dirlist_offn2 = {globals.offn2.name};
globals.offn3 = dir([data_loc,'*_B-3_v2.wav']);dirlist_offn3 = {globals.offn3.name};
globals.offn4 = dir([data_loc,'*_B-4_v2.wav']);dirlist_offn4 = {globals.offn4.name};
globals.offn5 = dir([data_loc,'*_B-5_v2.wav']);dirlist_offn5 = {globals.offn5.name};
globals.offn6 = dir([data_loc,'*_B-6_v2.wav']);dirlist_offn6 = {globals.offn6.name};
globals.offn7 = dir([data_loc,'*_B-7_v2.wav']);dirlist_offn7 = {globals.offn7.name};
globals.offn8 = dir([data_loc,'*_B-8_v2.wav']);dirlist_offn8 = {globals.offn8.name};
globals.offn9 = dir([data_loc,'*_B-9_v2.wav']);dirlist_offn9 = {globals.offn9.name};
globals.offn10 = dir([data_loc,'*_B-10_v2.wav']);dirlist_offn10 = {globals.offn10.name};
globals.offp1 = dir([data_loc,'*_B1_v2.wav']);dirlist_offp1 = {globals.offp1.name};
globals.offp2 = dir([data_loc,'*_B2_v2.wav']);dirlist_offp2 = {globals.offp2.name};
globals.offp3 = dir([data_loc,'*_B3_v2.wav']);dirlist_offp3 = {globals.offp3.name};
globals.offp4 = dir([data_loc,'*_B4_v2.wav']);dirlist_offp4 = {globals.offp4.name};
globals.offp5 = dir([data_loc,'*_B5_v2.wav']);dirlist_offp5 = {globals.offp5.name};
globals.offp6 = dir([data_loc,'*_B6_v2.wav']);dirlist_offp6 = {globals.offp6.name};
globals.offp7 = dir([data_loc,'*_B7_v2.wav']);dirlist_offp7 = {globals.offp7.name};
globals.offp8 = dir([data_loc,'*_B8_v2.wav']);dirlist_offp8 = {globals.offp8.name};
globals.offp9 = dir([data_loc,'*_B9_v2.wav']);dirlist_offp9 = {globals.offp9.name};
globals.offp10 = dir([data_loc,'*_B10_v2.wav']);dirlist_offp10 = {globals.offp10.name};



end