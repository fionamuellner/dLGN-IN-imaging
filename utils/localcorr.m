function [Cn] = localcorr(data,chunks,rem,filters,d2,d3)

% Calculate local correlation matrix across time intervals for ROI detection (see ROIdetect.m),
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% data: raw data, 1st dim: time, 2nd dim: space
% chunks: data intervals to calculate correlations (default: per 90 seconds)
% rem: binary vector indicating datapoints to treat as NaN
% filters: filter length for running average (default: 0.5 seconds)
% d2,d3: raw data dimensions in space
% __________________________________
% Outputs:
%
% Cn: matrix of local correlations; dimension 1,2: x,y coordinates, 3: time interval
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

Cn=0;
d1=size(data,1);
if nargin<5 
    if length(size(data))>2
        d2=size(data,2);
        d3=size(data,3);
        data=reshape(data,[d1,d2*d3]);
    else
        disp('Too few inputs provided.');
        Cn=[];
    end
elseif size(data,2)~=d2*d3
    disp('Dimensions mismatch.');
    Cn=[];
end

if ~isempty(Cn)
    Cn=zeros(d2,d3,floor(d1/chunks));
    for si=1:floor(d1/chunks)
        image=single(data((si-1)*chunks+1:min(si*chunks,d1),:));
        image=image(rem((si-1)*chunks+1:min(si*chunks,d1))==0,:);
        le=size(image,1);
        t1=tic;
        if le>chunks/2
            [Cn(:,:,si)]=neighborcorr(image,le,d2,d3,filters,1,3);
            t2=toc(t1);
            fprintf('Time for local correlations: %3.1f seconds\n',t2);
        end
    end
end
