function output = medianfilter2(input,r,mode)

% Medianfilter for 2-dimensional matrices.
% __________________________________
% Inputs:
%
% r: radius
% mode: 'fast' or 'accurate'
%
% __________________________________
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


% The following two methods give equivalent results, 
% except for reduced edge effects:
if nargin<3
    mode='fast';
end

if strcmp(mode,'fast')
    % Using the ordfilt2 function (much faster):
    if r~=0
        cmat = circmask(r);
        msize = (sum(sum(cmat))+1)/2;

        % To prevent edge effects:
        input = vertcat(repmat(input(1,:),msize,1),input,repmat(input(end,:),msize,1));
        input = horzcat(repmat(input(:,1),1,msize),input,repmat(input(:,end),1,msize));
        temp1 = ordfilt2(input,msize,cmat);
        temp2 = temp1(:,msize+1:end-msize);
        output = temp2(msize+1:end-msize,:);
    else
        output = input;
    end
else

    % 2-dimensional medianfilter from scratch (accurate but slow):
    k = size(input,1);
    l = size(input,2);
    cmat = circmask(r);
    s = (size(cmat,1)-1)/2; % radius of the circular mask

    % for the mean filter response, convolution is faster:
    % conv2(output,cmat);

    extmat = zeros(k+4*s,l+4*s);
    extmat(2*s+1:k+2*s,2*s+1:l+2*s) = input;
    
    % to prevent edge effects:
    extmat(1:2*s,2*s+1:l+2*s) = repmat(input(1,:),2*s,1);
    extmat(end-2*s+1:end,2*s+1:l+2*s) = repmat(input(end,:),2*s,1);
    extmat(2*s+1:k+2*s,1:2*s) = repmat(input(:,1),1,2*s);
    extmat(2*s+1:k+2*s,end-2*s+1:end) = repmat(input(:,end),1,2*s);
    extmat(1:2*s,1:2*s) = input(1,1);
    extmat(end-2*s+1:end,end-2*s+1:end) = input(end,end);
    extmat(1:2*s,end-2*s+1:end) = input(1,end);
    extmat(end-2*s+1:end,1:2*s) = input(end,1);
    for i = 1:k+2*s
        for j = 1:l+2*s
            tempmat=extmat(i:i+2*s,j:j+2*s);
            tempvect = tempmat(logical(cmat));
            extmat(i,j) = my_nanmedian(single(tempvect));
        end   
    end
    output = extmat(s+1:k+s,s+1:l+s);
end

