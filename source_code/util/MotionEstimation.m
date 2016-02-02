function [motionMap, uv] = MotionEstimation(img1, img2, blockSize, maxSearchRadius, lambda_reg, bVisualize)
%
%       [motionMap, uv] = MotionEstimation(img1, img2, blockSize, maxSearchRadius, lambda_reg, bVisualize)
%
%       This computes motion estimation between frames
%
%       input:
%         - img1: source
%         - img2: target
%         - blockSize: size of the block
%         - maxSearchRadius: search size in blocks
%         - lambda_reg: regularization coefficient
%         - bVisualize: 
%
%       output:
%         - motionMap: motion map for each pixel
%
%     Copyright (C) 2013-15  Francesco Banterle
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

if(~exist('blockSize', 'var'))
    nPixels = r * c;
    blockSize = max([2^ceil(log10(nPixels)), 4]);
end

if(~exist('maxSearchRadius', 'var'))
    maxSearchRadius = 2; %size in blocks
end

if(~exist('maxSearchRadius', 'var'))
    maxSearchRadius = 2; %size in blocks
end

if(~exist('lambda_reg', 'var'))
    lambda_reg = 0;
end

shift = round(blockSize * maxSearchRadius);

block_r = ceil(r / blockSize);
block_c = ceil(c / blockSize);

motionMap = zeros(r, c, 3);

uv = zeros(block_r, block_c, 4);

k_vec = [];
l_vec = [];
n_vec = [];
for k=(-shift):shift
	for l=(-shift):shift
        k_vec = [k_vec, k];
        l_vec = [l_vec, l];
        n_vec = [n_vec, abs(k) + abs(l)];
    end
end

n_vec = n_vec / (2 * shift);

vec_n = length(k_vec);

for i=1:block_r   
    for j=1:block_c     
        dx = 0;
        dy = 0;
        err = 1e20;
                
        i_b = (i - 1) * blockSize + 1;
        i_e = min([i_b + blockSize - 1, r]);
        j_b = (j - 1) * blockSize + 1;
        j_e = min([j_b + blockSize - 1, c]);
        
        for p=1:vec_n
            i_b1 = i_b - k_vec(p);
            i_e1 = i_e - k_vec(p);            
            j_b1 = j_b - l_vec(p);
            j_e1 = j_e - l_vec(p);

            i_b2 = i_b + k_vec(p);
            i_e2 = i_e + k_vec(p);            
            j_b2 = j_b + l_vec(p);
            j_e2 = j_e + l_vec(p);

            if( (i_b1 > 0) && (j_b1 > 0) && (i_e1 <= r) && (j_e1 <= c) &&...
                (i_b2 > 0) && (j_b2 > 0) && (i_e2 <= r) && (j_e2 <= c))
                tmp_err = abs(img1(i_b1:i_e1, j_b1:j_e1, :) - img2(i_b2:i_e2, j_b2:j_e2, :));
                tmp_err = mean(tmp_err(:)) + lambda_reg * n_vec(p);
                
                if(tmp_err < err)
                    err = tmp_err;
                    dx = l_vec(p);
                    dy = k_vec(p);
                end 
            end
        end
        
        motionMap(i_b:i_e,j_b:j_e,1) = dx;
        motionMap(i_b:i_e,j_b:j_e,2) = dy;
        motionMap(i_b:i_e,j_b:j_e,3) = err;
        
        uv(i, j, 1) = (j_b + j_e) / 2.0;
        uv(i, j, 2) = (i_b + i_e) / 2.0;
        uv(i, j, 3) = dx;
        uv(i, j, 4) = dy;
    end
end

if(bVisualize > 0)
    figure(bVisualize)    
    quiver(uv(:, :, 1), r - uv(:, :, 2) + 1, uv(:, :, 3), -uv(:, :, 4));
end

end