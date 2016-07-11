%
%      This HDR Toolbox demo creates an HDR radiance map:
%	   1) Read a stack of LDR images
%	   2) Read exposure values from the EXIF
%	   3) Estimate the Camera Response Function (CRF)
%	   4) Build the radiance map using the stack and stack_exposure
%	   5) Save the radiance map in .hdr format
%	   6) Show the tone mapped version of the radiance map
%
%       Author: Francesco Banterle
%       Copyright 2015 (c)
%

clear all;

name_folder = 'stack';
format = 'jpg';

disp('1) Read a stack of LDR images');
[stack, norm_value] = ReadLDRStack(name_folder, format, 1);

disp('2) Read exposure values from the exif');
stack_exposure = ReadLDRStackInfo(name_folder, format);

disp('3) Estimate the Camera Response Function (CRF)');
[lin_fun, ~] = DebevecCRF(stack, stack_exposure);
h = figure(1);
set(h, 'Name', 'The Camera Response Function (CRF)');
plot(0:255, lin_fun(:,1), 'r', 0:255, lin_fun(:,2),'g', 0:255, lin_fun(:,3), 'b');

disp('4) Build the radiance map using the stack and stack_exposure');
imgHDR = BuildHDR(stack, stack_exposure, 'LUT', lin_fun, 'Deb97', 'log');

disp('5) Save the radiance map in the .hdr format');
hdrimwrite(imgHDR, 'stack_hdr_image.exr');

disp('6) Show the tone mapped version of the radiance map with gamma encoding');
h = figure(2);
set(h, 'Name', 'Tone mapped version of the built HDR image');
GammaTMO(ReinhardTMO(imgHDR, 0.18), 2.2, 0, 1);


disp('7) Estimate the polynomial Camera Response Function (CRF)');
[lin_fun_poly, pp] = MitsunagaNayarCRF(stack, stack_exposure, -3, 1000, true);
h = figure(3);
set(h, 'Name', 'The polynomial Camera Response Function (CRF)');
plot(0:255, lin_fun_poly(:,1), 'r', 0:255, lin_fun_poly(:,2),'g', 0:255, lin_fun_poly(:,3), 'b');

disp('8) Build the radiance map using the stack and stack_exposure');
imgHDRpoly = BuildHDR(stack, stack_exposure, 'poly', pp, 'Deb97', 'log');

disp('9) Save the radiance map in the .hdr format');
hdrimwrite(imgHDRpoly, 'stack_hdr_image_poly.exr');

disp('10) Show the tone mapped version of the radiance map with gamma encoding');
h = figure(4);
set(h, 'Name', 'Tone mapped version of the built polynomial HDR image');
GammaTMO(ReinhardTMO(imgHDRpoly, 0.18), 2.2, 0, 1);
