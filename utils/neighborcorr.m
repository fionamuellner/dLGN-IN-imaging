function [C]=neighborcorr(data,d1,d2,d3,filters,downs,ordfilts)

% Calculates for each pixel its correlation with the average of its four direct neighbor pixels; 
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% data: raw data, 1st dim: time, 2nd dim: space
% d2,d3: raw data dimensions in space
% filters: filter length for running average (default: 0.5 seconds) 
% downs: factor for downsampling in time
% ordfilts: filter size of optional pre-processing medianfilter in the time-domain 
% __________________________________
% Outputs:
%
% C: correlation matrix 
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt


if nargin<5
    filters=1;
end
if nargin<6
    downs=1;
end
if nargin<7
    ordfilts=1;
end
if downs>1
    data=downsamp(data',downs);
    data=data';
end
if ordfilts>1
   if ceil(ordfilts/2)-(ordfilts+1)/2~=0
       disp('ordfilts was rounded to the nearest smaller odd number:');
       ordfilts=(floor(ordfilts/2))*2+1;
       disp(ordfilts)
   end
   data=ordfilt2(data,(ordfilts+1)/2,ones(ordfilts,1)); 
   data=data((ordfilts+1)/2:end-(ordfilts+1)/2+1,:);
   d1=d1-ordfilts+1;
end
C=zeros(d2,d3);
for i=2:d2-1
    temp=squeeze(single(data(:,i+[0:d3-1]*d2)));
    t1=filter(ones(filters,1),1,temp);
    t1=t1(filters:end,:);
    t2=zeros(d1-filters+1,d3);
    for j=2:d3-1
        sur=[i+[-1:1]+(j-2)*d2,i+[-1,1]+(j-1)*d2,i+[-1:1]+j*d2];
        temp=mean(data(:,sur),2);
        temp=filter(ones(filters,1),1,temp);
        t2(:,j)=temp(filters:end);
    end
    temp=corr(t1,t2);
    C(i,:)=diag(temp);
end
