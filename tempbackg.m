function [bgt_init] = tempbackg(proj,rad,n)

% Determine darkest pixels within a given radius for temporal background estimation;
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% proj: projected image
% rad: radius
% n: number of background pixels
% __________________________________
% Outputs:
%
% bgt_init: matrix containing background pixel coordinates; dimension 1: coordinates, 2/3: xy position in projected image
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

if nargin<4
    rad=20;
end
if nargin<5
    n=5;
end

d2=size(proj,1);
d3=size(proj,2);
bgt_init=zeros(n,d2,d3);
cm=circmask_part(rad,rad+1,rad+1,d2,d3);
[iy,ix]=find(cm);
for i=1:d2
    for j=1:d3
        iy_sur=iy+i-rad-1;
        ix_sur=ix+j-rad-1;
        f=find(iy_sur<1 | iy_sur>d2 | ix_sur<1 | ix_sur>d3);
        iy_sur(f)=[];
        ix_sur(f)=[];

        f=sub2ind([d2,d3],iy_sur,ix_sur);
        [~,si]=sort(proj(f),'ascend');
        bgt_init(:,i,j)=f(si(1:n));
    end
end

