function [motionMap, uv] = MotionEstimationBiDiHDR(img_prev, img_next, blockSize, maxSearchRadius, lambda_reg, bVisualize)
%
%        motionMap = MotionEstimationBiDiHDR(img_prev, img_next, blockSize)
%
%       This computes bi-directional motion estimation between HDR frames
%
%       input:
%         - img_prev: previous frame
%         - img_next: next frame
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
    [r, c, ~] = size(img_prev);
    nPixels = r * c;
    blockSize = max([2^ceil(log10(nPixels)), 4]);
end

if(~exist('maxSearchRadius', 'var'))
    maxSearchRadius = 2; %size in blocks
end

if(~exist('lambda_reg', 'var'))
    lambda_reg = 0;
end

[motionMap, uv] = MotionEstimation(log10(img_prev + 1e-6), log10(img_next + 1e-6), blockSize, maxSearchRadius, lambda_reg, bVisualize);

end