clear all
close all
clc

%input directory
nrrd_dir='/home/cc14/Desktop/wallthicknessMaps/wallthicknessMaps/901001';

%input files:

%segmentation
nrrd_file='wallthickness3SD.nrrd';
%File for ahortic valv
nrrd_valv='mvlayer2vox.nrrd';

npix=3; %number of voxels to take for delimiting PVs

%output dir and file name
out_dir='/home/cc14/Desktop/wallthicknessMaps/wallthicknessMaps/901001';
out_file='LA.inr';



[thickMap, meta] = nrrdread([nrrd_dir,'/',nrrd_file]);
[Valv, ~] = nrrdread([nrrd_dir,'/',nrrd_valv]);

clear('nrrd_dir','nrrd_file','nrrd_valv');


new_img=labeling(thickMap,Valv,npix);
spacing = arrayfun(@str2double,strsplit(meta.spacings,' '));
convertMat2Inr(new_img,spacing,[out_dir,'/',out_file]);

clear('out_dir','out_file');

