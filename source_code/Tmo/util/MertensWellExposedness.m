function We = MertensWellExposedness(img)
%
%
%        We = MertensWellExposedness(img)
%
% 
%     Copyright (C) 2010-15 Francesco Banterle
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

%sigma for the Well-exposedness weights.
sigma  = 0.2; %as in the original paper
sigma2 = 2.0 * sigma^2;

[r, c, col] = size(img);
We = ones(r,c );

for i=1:col
    We = We .* exp(-(img(:,:,i) - 0.5).^2 / sigma2);
end

end