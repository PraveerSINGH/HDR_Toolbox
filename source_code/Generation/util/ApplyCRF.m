function imgOut = ApplyCRF(img, table)
%
%
%        img = ApplyCRF(img, table)
%
%
%        Input:
%           -img: an LDR image with values in [0,1]
%           -table: three functions for remapping image pixels values
%
%        Output:
%           -imgOut: an LDR image with with values in [0,2^nBit - 1]
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

col = size(img, 3);

total_values = size(table, 1);
x = 0 : (total_values - 1);
imgOut = zeros(size(img));

for i=1:col
    imgOut(:,:,i) = interp1(table(:,i), x, img(:,:,i), 'nearest', 'extrap');
end

end