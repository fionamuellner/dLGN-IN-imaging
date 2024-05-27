function [mi1,mi2,newrois]=findLocalMax2d(data,prcthr,maxnum)

% Finds local maxima in 2D; 
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% data: raw 2D data
% prcthr: percentile cutoff for detection
% maxnum: maximum number of detected peaks
% __________________________________
% Outputs:
%
% mi1, mi2: coordinates of local maxima
% newrois: local environments of the maxima
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt


if nargin<3
    maxnum=10^4;
end
data=data-min(data(:));
count=1;
mi1=[];
mi2=[];
if nargin<2
   prcthr=99.9;
end
thr=prctile(data(:),prcthr);
data(isnan(data))=thr;
data(data<thr)=thr;
data=horzcat(data(:,1)*NaN,data,data(:,1)*NaN);
data=vertcat(data(1,:)*NaN,data,data(1,:)*NaN);
data1=(circshift(data,[-1,0])>=data)==1;
data2=(circshift(data,[1,0])>=data)==1;
data3=(circshift(data,[0,1])>=data)==1;
data4=(circshift(data,[0,-1])>=data)==1;
data5=(circshift(data,[1,1])>=data)==1;
data6=(circshift(data,[-1,1])>=data)==1;
data7=(circshift(data,[1,-1])>=data)==1;
data8=(circshift(data,[-1,-1])>=data)==1;
data=data(2:end-1,2:end-1);
data1=data1(2:end-1,2:end-1);
data2=data2(2:end-1,2:end-1);
data3=data3(2:end-1,2:end-1);
data4=data4(2:end-1,2:end-1);
data5=data5(2:end-1,2:end-1);
data6=data6(2:end-1,2:end-1);
data7=data7(2:end-1,2:end-1);
data8=data8(2:end-1,2:end-1);
[sortdata,si]=sort(data(:),'descend');
[~,invs]=sort(si);
mi=1;
mi1=zeros(maxnum,1);
mi2=zeros(maxnum,1);
newrois=cell(maxnum,1);
while ~isempty(mi) && sortdata(mi)>=thr && count<=maxnum
    [mi1(count,1),mi2(count,1)]=ind2sub(size(data),si(mi));
    roi=data1*0;
    roi(mi1(count),mi2(count))=1;
    [~,sur]=surround_1pix(roi);
    sur=sur.*(data<=data(mi1(count),mi2(count))).*(data>thr);
    while any(sur(:)) 
        roi=sur|roi;
        lm=max(max(data(sur==1)));
        [~,sur]=surround_1pix(roi);
        sur=sur.*(data<=lm).*(data>thr);
    end

    f=find(roi);
    roi(data<=0.3*(data(mi1(count,1),mi2(count,1)))+0.7*thr)=0;
    data(f)=thr;
    newrois{count,1}=roi;
    count=count+1;
    sortdata(invs(f))=thr;
    mi=find(sortdata>thr,1,'first');
end
mi1=mi1(1:count-1);
mi2=mi2(1:count-1);
newrois=newrois(1:count-1,1);