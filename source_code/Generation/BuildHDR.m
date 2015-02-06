function [imgHDR, lin_fun] = BuildHDR(stack, stack_exposure, lin_type, lin_fun, weightFun, bRobertson)
%
%       [imgHDR, lin_fun] = BuildHDR(stack, stack_exposure, lin_type, lin_fun, weightFun, bRobertson)
%
%       This function builds an HDR image from a stack of LDR images.
%
%        Input:
%           -stack: an input stack of LDR images. This has to be set if we
%           the stack is already in memory and we do not want to load it
%           from the disk using the tuple (dir_name, format).
%           -stack_exposure: an array containg the exposure time of each
%           image. Time is expressed in second (s).
%           -lin_type: the linearization function:
%                      - 'linear': images are already linear
%                      - 'gamma2.2': gamma function 2.2 is used for
%                                    linearisation;
%                      - 'sRGB': images are encoded using sRGB
%                      - 'tabledDeb97': a tabled RGB function is used for
%                                       linearisation passed as input in
%                                       lin_fun using Debevec and Malik 97
%                                       method
%           -lin_fun: it is the camera response function of the camera that
%           took the pictures in the stack. If it is empty, [], this
%           function will be estimated when lin_type == tabledDeb97.
%           -weight_type:
%               - 'all':   weight is set to 1
%               - 'hat':   hat function 1-(2x-1)^12
%               - 'Deb97': Debevec and Malik 97 weight function
%               - 'Gauss': Gaussian function as weight function.
%                          This function produces good results when some 
%                          under-exposed or over-exposed images are present
%                          in the stack.
%            -bRobertson: if it is set to 1 it enables the Robertson's
%             modification for assembling exposures for reducing noise.
%
%        Output:
%           -imgHDR: the final HDR image
%           -lin_fun:
%
%        Example:
%           This example line shows how to load a stack from disk:
%
%               stack = ReadLDRStack('stack_alignment', 'jpg');               
%               stack_exposure = ReadLDRExif('stack_alignment', 'jpg');
%               BuildHDR(stack, stack_exposure,'tabledDeb97',[],'Deb97');
%
%           In the case we previously loaded LDR images, in stack, and
%           their EXIF information, in stack_exposure, we have to use
%           the following line:
%               BuildHDR('','','tabledDeb97','Deb97',stack,stack_exposure);
%
%
%     Copyright (C) 2011-15  Francesco Banterle
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

if(~exist('bRobertson', 'var'))
    bRobertson = 0;
end

%is a weight function defined?
if(~exist('weightFun', 'var'))
    weightFun = 'all';
end

%is the linearization type of the images defined?
if(~exist('lin_type', 'var'))
    lin_type = 'gamma2.2';
end

%do we have the inverse camera response function?
if(~exist('lin_fun', 'var'))
    lin_fun = [];
end

if(isempty(stack) || isempty(stack_exposure))
    error('The stack is set empty!');
end

%do we have a camera response function?
bFun = 0;
if(~isempty(lin_fun))
    bFun = (length(lin_fun) == 256);
end

%is the inverse camera function ok? Do we need to recompute it?
if((strcmp(lin_type, 'tabledDeb97') == 1) && bFun)
    lin_fun = ComputeCRF(stack, stack_exposure);        
end

%combining the LDR images
imgHDR = double(CombineLDR(stack, stack_exposure, lin_type, lin_fun, weightFun, bRobertson));

end