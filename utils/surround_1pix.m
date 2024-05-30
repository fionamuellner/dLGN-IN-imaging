function [sur,suroi]=surround_1pix(roi)

% Calculates the 1-pixel surround of a roi mask.
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

f=find(roi>0);
suroi=0*roi;
[y,x]=ind2sub(size(roi),f);
i=1;
while i<=p
    ny=vertcat(y,y-1,y+1,y,y,y-1,y-1,y+1,y+1);
    nx=vertcat(x,x,x,x-1,x+1,x-1,x+1,x-1,x+1);
    y=ny;
    x=nx;
    i=i+1;
end
ny=max(ny,1);
ny=min(ny,size(roi,1));
nx=max(nx,1);
nx=min(nx,size(roi,2));
sur=sub2ind(size(roi),ny,nx);
sur=unique(sur);
sur=setdiff(sur,f);
suroi(sur)=1;
