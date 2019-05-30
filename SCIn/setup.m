
addpath(genpath(pwd));

a=['root=''' pwd ''';'];
clear rootdir.m;
edit rootdir.m;
FID=fopen('rootdir.m','w');
fprintf(FID, '%s', a);
fclose(FID);