function [stackVP] = BuildVirtualPhotographs(img, deltaT, stack, stack_exposure, lin_fun)
%
%
%      [imgOut, deltaT] = BuildVirtualPhotograph(img, deltaT, stack, stack_exposure, lin_fun)
%
%
%       Input:
%           -img: input HDR image
%
%           -deltaT: array of exposure time for each virtual photograph
%
%           -stack: an input stack of LDR images. This has to be set if we
%           the stack is already in memory and we do not want to load it
%           from the disk using the tuple (dir_name, format).
%           If the stack is a single or dobule, values are assumed to be in
%           the range [0,1].
%
%           -stack_exposure: an array containg the exposure time of each
%           image in the stack. Time is expressed in second (s).
%
%           -lin_fun: it is the camera response function of the camera that
%           took the pictures in the stack. If it is empty, [], and 
%           type is 'LUT' it will be estimated using Debevec and Malik's
%           method.
%
%       Output:
%
%           -stackVP: stack of Virtual Photographs corresponding to array of deltaT
%
%     Copyright (C) 2016  Praveer SINGH
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%     The paper describing this technique is:
%     "Virtual Photograph Based Saliency Analysis of High Dynamic Range Images"
% 	  by Gao, Xihe and Brooksy, Stephen and Arnold, Dirk V.
%     in Computational Aesthetics in Graphics, Visualization, and Imaging
%     2013
%

%do we have the inverse camera response function?
if(~exist('lin_fun', 'var'))
    lin_fun = [];
end

%deltaT checks
if (isempty(deltaT))
    error('deltaT (exposure time) is empty')
end

%stack exposure checks
if(isempty(stack) || isempty(stack_exposure))
    error('The stack is set empty!');
end

stack_exposure_check = unique(stack_exposure);

if(length(stack_exposure) ~= length(stack_exposure_check))
    error('The stack contains images with the same exposure value. Please remove these duplicated images!');
end

%recomute the inverse Camera response function if found empty
if(isempty(lin_fun))
    [lin_fun, ~] = DebevecCRF(single(stack) / scale, stack_exposure);        
end

%calibrating luminance using logarithmic mean
Lav = logMean(img);
Lcalib = img/Lav;
y = 0:1/255:1;
[m,n,o] = size(img);
stackVP = zeros(m,n,o,length(deltaT));
for i = 1:length(deltaT)
    
    for j = 1:3
    
        Lcalib_perchannel = Lcalib(:,:,j);
        xq = log(Lcalib_perchannel(:)*deltaT(i)); %xq are target HDR query points to be interpolated using CRF 
        
        % g = ln(f-1) where f-1 (lin_fun) is the inverse CRF function
        %Display value computed from the inverse of this lookup table 
        %Z(x,y) = g-1(ln(L(x,y)*deltaT))
    
        g = log(lin_fun(:,1)); %we chose R-channel CRf since it has largest range for this stack
        vq1 = interp1(g,y,xq,'nearest','extrap');
        vq1 = reshape(vq1,m,n);
        stackVP(:,:,j,i) = vq1;
    end
end

end

