function output = circmask(r)

% Binary mask with 1 for any pixel with distance from center < r+sqrt(0.5). 
% Roughly matches ImageJ circmask definition.
%
% __________________________________
%
% Copyright 2024 Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel
% see license.txt
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%

s = ceil(r);

output = zeros(2*s+1,2*s+1);
for i = 1:2*s+1
    for j = 1:2*s+1
        if sqrt((s+1-i)^2 + (s+1-j)^2) < r+sqrt(0.5)
            output(i,j) = 1;
        end
    end
end
