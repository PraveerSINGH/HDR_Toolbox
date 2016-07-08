function image_target = reExpose(img_source, source_exposure, target_exposure, lin_type, lin_fun)
%
%       image_target = reExpose(img_source, source_exposure, target_exposure, lin_type, lin_fun)
%       
%       This function re-exposes an LDR image set        
%       with respect to the target image.
%
%        Input:
%           -img_source:
%           -source_exposure:
%           -target_exposure:
%           -lin_type:
%           -lin_fun:
%
%        Output:
%           -image_target:
%
%     Copyright (C) 2015  Damla Ezgi Akcora
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

Z = RemoveCRF(img_source, lin_type, lin_fun);

image_target = (Z * target_exposure) / source_exposure;
image_target = ClampImg(image_target.^(1.0 / 2.2), 0.0, 1.0);

end