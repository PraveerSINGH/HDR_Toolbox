%
%      This HDR Toolbox demo creates virtual photographs from any HDR image:
%	   1) Read a stack of LDR images
%	   2) Read exposure values from the EXIF
%	   3) Estimate the Camera Response Function (CRF)
%	   4) Read the targetted HDR image
%	   5) Compute ambient luminance
%	   6) Calibrate luminance of image using ambient luminance
%      7) Build the Virtual Photographs using delta-T and CRF
%
%       Author: Praveer SINGH
%       Copyright 2016 (c)
%

clear all;
close all;

name_folder = 'stack';
format = 'jpg';

disp('1) Read a stack of LDR images');
[stack, norm_value] = ReadLDRStack(name_folder, format, 1);

disp('2) Read exposure values from the exif');
stack_exposure = ReadLDRStackInfo(name_folder, format);

disp('3) Estimate the Camera Response Function (CRF)');
[lin_fun, ~] = DebevecCRF(stack, stack_exposure);

disp('4) Reading HDR image');
img = hdrimread('Bottles_Small.hdr');

disp('5) Reading stack of deltaT');
deltaT = [1/4,1/2,1,2,4];

disp('6) Building stack of virtual Photographs using array of deltaT');
stackVP = BuildVirtualPhotographs(img, deltaT, stack, stack_exposure, lin_fun);

for i = 1:length(deltaT)    
    
    figure
    imshow(stackVP(:,:,:,i));
    imwrite(stackVP(:,:,:,i),['virtual_photograph_',num2str(i),'.png']);
end
