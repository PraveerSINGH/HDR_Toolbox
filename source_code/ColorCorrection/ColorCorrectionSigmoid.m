function S = ColorCorrectionSigmoid(L_cone, sr_n, sigma_cone, B_cone)
%
%         S = ColorCorrectionSigmoid(L_cone, sr_n, sigma_cone, B_cone)
%
%        This function computes color correction for sigmoid.
%
%        Input:
%           -L_cone: the luminance value for cones
%           -sr_n: sigmoid parameter (typically set to 0.73)
%           -sigma_cone: saturation parameter for cones
%           -B_cone: bleaching parameter for cones
%
%        Output:
%           -S: the saturation value
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
%
%     The paper describing this technique is:
%     "Time-Dependent Visual Adaptation For Fast Realistic Image Display"
% 	  by Sumanta N. Pattanaik, Jack Tumblin, Hector Yee, Donald P. Greenberg
%     in SIGGRAPH 2000
%

S = sr_n .* B_cone .* (L_cone.^sr_n) .* (sigma_cone.^sr_n);
S = S ./ (L_cone.^sr_n + sigma_cone.^sr_n).^2;

end

