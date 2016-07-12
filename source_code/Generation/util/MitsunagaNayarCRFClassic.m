function [pp, err] = MitsunagaNayarCRFClassic(stack_samples, stack_exposure, N)
%
%       [pp, err] = MitsunagaNayarCRFClassic(stack_samples, stack_exposure, N)
%
%       This function computes camera response function using Mitsunaga and
%       Nayar method.
%
%        Input:
%           -stack_samples: a stack of samples from LDR images
%           -nSamples: number of samples for computing the CRF
%           -N: polynomial degree of the inverse CRF
%
%        Output:
%           -pp: a polynomial encoding the inverse CRF
%           -err: the error
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

col = size(stack_samples, 3);

Q = length(stack_exposure);

Mmax = 1.0;

%recovering the CRF
function d = MN_d(c, q, n)
    q_p = q + 1;

    M_q   = stack_samples(:, q  , c) ;
    M_q_p = stack_samples(:, q_p, c) ;
    
    indx = find(M_q > 0.0 & M_q_p > 0.0);
        
    d = M_q(indx).^n - R(q) * (M_q_p(indx).^n);
end

pp = zeros(N + 1, col);

err = 0.0;

R0 = zeros(Q - 1, 1);
for q=1:(Q - 1)
    R0(q) = stack_exposure(q) / stack_exposure(q + 1);
end

for c=1:col
    
    R = R0;

    A = zeros(N, N);
    b = zeros(N, 1);

    for i=1:N
        %init A
        for j=1:N
            for q=1:(Q - 1)
                delta  = MN_d(c, q, j - 1) - MN_d(c, q, N);
                A(i,j) = A(i,j) + sum(MN_d(c, q, i - 1) .* delta);
            end

        end

        %init b
        for q=1:(Q - 1)
            b(i) = b(i) - sum(Mmax * MN_d(c, q, i - 1) .* MN_d(c, q, N));
        end
    end  

    coeff = A \ b;    
    coeff_n = Mmax - sum(coeff);

    pp(:,c) = flip([coeff; coeff_n]);
      
    %compute err
    for q=1:(Q - 1)
        s1 = stack_samples(:, q    , c);
        s2 = stack_samples(:, q + 1, c);   
        
        indx = find(s1 > 0.0 & s2 > 0.0);
        
        e1 = polyval(pp(:,c), s1(indx) );
        e2 = polyval(pp(:,c), s2(indx) );
        
        err = err + sum((e1 - R(q) * e2).^2);
    end
end

end