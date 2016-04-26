function [motionMap, uv] = MotionEstimationHDR(img1, img2, blockSize, maxSearchRadius, lambda_reg, bVisualize)
%
%        motionMap = MotionEstimationHDR(img1, img2, blockSize)
%
%       This computes motion estimation between HDR frames
%
%       input:
%         - img1: source
%         - img2: target
%         - blockSize: size of the block
%         - maxSearchRadius: search size in blocks
%         - lambda_reg: regularization coefficient
%         - bVisualize: if it is set to 1 it visualizes the motion field 
%
%       output:
%         - motionMap: motion map for each pixel
%         - uv: motion map for each pixel
%
%     Copyright (C) 2013-16  Francesco Banterle
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

[r, c, ~] = size(img1);

if(~exist('bVisualize', 'var'))
    bVisualize = 0;
end

bAuto = 0;
if(~exist('blockSize', 'var'))
    bAuto = 1;
else
    bAuto = strcmp('blockSize', 'auto');
end

if(bAuto)
    nPixels = r * c;
    blockSize = max([2^ceil(log10(nPixels)), 4]);
end

if(~exist('maxSearchRadius', 'var'))
    maxSearchRadius = 2; %size in blocks
end

if(~exist('lambda_reg', 'var'))
    lambda_reg = 0;
end

[motionMap, uv] = MotionEstimation(log10(img1 + 1e-6), log10(img2 + 1e-6), blockSize, maxSearchRadius, lambda_reg, bVisualize);

end