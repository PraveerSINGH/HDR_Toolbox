function lin_fun = RobertsonCRF(stack, stack_exposure, max_iterations, err_threshold, bNormalize, bPolyFit)
%
%       lin_fun = RobertsonCRF(stack, stack_exposure, max_iterations, err_threshold, bNormalize, bPolyFit)
%
%       This function computes camera response function using Mitsunaga and
%       Nayar method.
%
%        Input:
%           -stack: a stack of LDR images. If the stack is a single or
%           double values are assumed to be in [0,1].
%           -stack_exposure: an array containg the exposure time of each
%           image. Time is expressed in second (s)
%           -max_iterations: max number of iterations
%           -err_threshold: threshold after which the function can stop
%           -bNormalize: if 1 it enables function normalization
%           -bPolyFit: if it is set to 1, a polynomial fit will be computed;
%           otherwise it will not (default)
%
%        Output:
%           -lin_fun: tabled function
%           -max_lin_fun: maximum value of the inverse CRF
%
%     Copyright (C) 2016  Francesco Banterle
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

if(isempty(stack))
    error('RobertsonCRF: a stack cannot be empty!');
end

if(isempty(stack_exposure))
    error('RobertsonCRF: a stack_exposure cannot be empty!');
end

if(~exist('err_threshold', 'var'))
    err_threshold = 1e-7;
end

if(~exist('bNormalize', 'var'))
    bNormalize = 1;
end

if(~exist('max_iterations', 'var'))
    max_iterations = 5;
end

if(~exist('bPolyFit', 'var'))
    bPolyFit = 0;
end

if(isa(stack, 'double'))
    stack = uint8(round(stack * 255));
end

if(isa(stack, 'single'))
    stack = uint8(round(stack * 255));
end

if(isa(stack, 'uint16'))
    stack = uint8(stack / 255);
end

col =  size(stack, 3);
lin_fun = zeros(256, col);
 
for i=1:col
    lin_fun(:, i) = (0:255) / 255;
end

for i=1:max_iterations
    lin_fun_prev = lin_fun;
    
    %iCRF normalization step
    for j=1:col
        lf = lin_fun(:, j);
        
        [~, index_min] = min(lf(lf > 0.0));
        [~, index_max] = max(lf(lf > 0.0));
        
        index = index_min + round((index_max - index_min) / 2);
        
        mid = lf(index);
        
        if(mid == 0.0)
            while(k < 256 && lf(index) == 0.0)
                index = index + 1;
            end
            
            mid = lf(index);
        end
        
        if(mid > 0.0)
            lin_fun(:,j) = lin_fun(:,j) / mid;
        end
    end

    %update X
    x_tilde = Update_X(stack, stack_exposure, lin_fun, scale);
    
    %update the iCRF
    lin_fun = Update_lin_fun(x_tilde, stack, stack_exposure, lin_fun); 
            
    %compute error
    delta = (lin_fun_prev - lin_fun).^2;
    err = mean(delta(:));
    
    if(err < err_threshold)
        break;
    end
    
    disp([i, err]);
end

% %poly-fit (rational)
% if(bPolyFit)
%     x = (0:255) / 255;
%     ft = fittype( 'rat33' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = ones(7, 1);
% 
%     for i=1:col    
%         [xData, yData] = prepareCurveData(x', lin_fun(:,i));
%         [fit_ret, ~] = fit( xData, yData, ft, opts );    
%         lin_fun(:,i) = feval(fit_ret, x');
%     end
% end

end

function imgOut = Update_X(stack, stack_exposure, lin_fun)
    [r, c, col, n] = size(stack);

    %for each LDR image...
    imgOut    = zeros(r, c, col, 'single');
    totWeight = zeros(r, c, col, 'single');
    
    for i=1:n
        tmpStack = ClampImg(single(stack(:,:,:,i)) / 255.0, 0.0, 1.0);

        %computing the weight function    
        weight  = WeightFunction(tmpStack, 'Robertson', 0);
        
        tmpStack = tabledFunction(tmpStack, lin_fun); 

        %summing things up...
        t = stack_exposure(i);    
        if(t > 0.0)                
            imgOut = imgOut + (weight .* tmpStack) * t;
            totWeight = totWeight + weight * t * t;
        end
    end

    saturation = 1e-4;

    imgOut = imgOut ./ totWeight;
    imgOut(totWeight < saturation) = -1.0;
end

function f_out = Update_lin_fun(x_tilde, stack, stack_exposure, lin_fun)
    col = size(x_tilde, 3);

    n = length(stack_exposure);
    f_out = zeros(size(lin_fun));

    for i=1:col
        tmp_x_tilde = x_tilde(:,:,i);
                       
        for j=0:255
            jp1 = j + 1;

            card_m = 0;
            for k=1:n
                t = stack_exposure(k);
                
                tmp = stack(:,:,i,k);

                tmp_x_tilde_ind = tmp_x_tilde(tmp == j & tmp_x_tilde > 0.0);

                ind = find(tmp_x_tilde_ind >= 0.0);  
                f_out(jp1, i) = f_out(jp1, i) + t * sum(tmp_x_tilde_ind(ind));

                card_m = card_m + length(ind);
            end

            if(card_m > 0.0)
                f_out(jp1, i) = f_out(jp1, i) / card_m;
            else
                f_out(jp1, i) = 0.0;
            end    
        end
    end
end