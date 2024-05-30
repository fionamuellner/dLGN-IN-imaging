function [roi,sur,orig,C,Cfilt,Cbg,M,nroi,maxnroi,copyInd,Ind,mCn,corrthr] = flexROIdetect(Cn,maxdata,data,chunks,rem,bgt,maxradius,maxnroi,filters,file_name,root)

% ROI detection routine based on correlation-based seeding, adaptive thresholding and expansion;
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% Cn: local correlation matrix to seed ROIs (see localcorr.m); optional 3rd dimension provides corr per time intervals (default: per 90 seconds)
% maxdata: projected data to add seeds for optimized detection of sparsely responding cells (default: 95th percentile projection, see myprctile2D.m)
% data: motion corrected data; dimension 1: time, 2/3: space
% chunks: data intervals to calculate correlations (default: per 90 seconds)
% rem: binary vector indicating datapoints to treat as NaN
% bgt: matrix for temporal background calculation (see tempbackg.m)
% maxradius: maximum ROI radius (default: 15 um)
% maxnroi: maximum number of ROIs detected (default: 10000)
% filters: filter length for running average (default: 0.5 seconds)
% file_name: filename to save code
% root: folder to save code 
% __________________________________
% Outputs:
%
% roi: spatial ROIs 
% sur: ROI surrounds
% orig: ROI origins
% C: temporal footprints
% Cfilt: filtered temporal footprints
% Cbg: temporal background
% M: origin coordinates
% nroi: number of detected ROIs
% copyInd: priority matrix for initialization of ROIs
% Ind: remaining priority matrix (with intialized ROIs removed)
% mCn: matrix of local correlations (projected across time intervals)
% corrthr: adaptive correlation threshold per ROI
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



warning('off','images:imshow:magnificationMustBeFitForDockedFigure');

% Save current file version together with the analysis results:
FileNameAndLocation=[mfilename('fullpath')];
FileName=mfilename();
newbackup=sprintf('%s\\%s_%s.m',root,file_name,FileName);
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

t1=tic;

d1=size(data,1);
d2=size(data,2);
d3=size(data,3);

% Define matrix indicating correlated pixels:
mCn=prctile(Cn(:,:,1:end),95,3);
Ind=mCn;

% Remove pixels not significantly above baseline correlations:
stdest=prctile(movingstd(Ind(~isnan(Ind)),20),2);
indthr=prctile(Ind(~isnan(Ind)),20)+5*stdest;
Ind(Ind<indthr)=0;

% Remove singular pixels:
Ind2=ordfilt2(Ind,8,ones(3,3))>0;
Ind=Ind.*Ind2;

% Add maxima of maxdata to include structures with low correlations:
Ind3=maxdata;
Ind3(maxdata<prctile(maxdata(:),99))=0;
[mi1,mi2]=findLocalMax2d(ordfilt2(Ind3,25,ones(7,7)),50);
s=find(mi1>10 & mi1<d2-9 & mi2>10 & mi2<d3-9);
mi1=mi1(s);
mi2=mi2(s);
cut=zeros(21,21,length(mi1));
for i=1:length(mi2)
    cut(:,:,i)=maxdata(mi1(i)+[-10:10],mi2(i)+[-10:10]);
end
temp=mean(cut,3);


if any(temp(:))
    c=normxcorr2(temp,maxdata);
    c=c(11:end-10,11:end-10);
    c(c<0.5)=0;
    [mi1,mi2]=findLocalMax2d(c,99.5);
    for i=1:length(mi2)
        Ind(mi1(i),mi2(i))=max(mCn(mi1(i),mi2(i)),indthr);
    end
end

copyInd=Ind;
nind=my_nansum(Ind(:));

data=reshape(data,[d1,d2*d3]);

% Maximum number of ROIs to be detected:
nroi=maxnroi;
roi=cell(nroi,1);
sur=cell(nroi,1);
orig=cell(nroi,1);
C=NaN*zeros(d1,nroi);
Cfilt=NaN*zeros(d1,nroi);
Cbg=NaN*zeros(d1,nroi);
M=zeros(nroi,2);
data=reshape(data,[d1,d2,d3]);
corrthr=zeros(nroi,1);
nchunks=max(1,floor(d1/chunks));

% Define the dataintervals for each chunk:
int=cell(nchunks,1);
L=zeros(nchunks,1);
for sj=1:nchunks-1
    int{sj,1}=(sj-1)*chunks+1:min(sj*chunks,d1);
end
sj=nchunks;
int{sj,1}=(sj-1)*chunks+1:d1;
for sj=1:nchunks
    int{sj,1}=int{sj,1}(rem(int{sj,1})==0);
    L(sj)=length(int{sj,1});
end
int=int(L>chunks/2,:);
int{1,1}=int{1,1}((filters-1)/2+1:end);
nchunks=size(int,1);
int{nchunks,1}=int{nchunks,1}(1:end-filters);




ri=1;
fprintf('%1.4f',ceil(10^4*(1-my_nansum(Ind(:))/nind))/10^4);
while any(Ind(:)) && ri<=nroi
    roi{ri,1}=zeros(d2,d3,'logical');
    [~,mi]=max(Ind(:));
    [mi1,mi2]=ind2sub(size(Ind),mi);
    roi{ri,1}(mi1,mi2)=1;
    Ind(mi1,mi2)=NaN;
    
    % initialize surround for each ROI
    sur{ri,1}=[];
    for j=-1:1
        for k=-1:1
            if mi1+j>0 && mi1+j<=d2 && mi2+k>0 && mi2+k<=d3 && ...
                    roi{ri,1}(mi1+j,mi2+k)==0 && ~isnan(Ind(mi1+j,mi2+k))
                sur{ri,1}=vertcat(sur{ri,1},[mi1+j,mi2+k]);
            end
        end
    end
    temp=single(data(:,bgt(:,mi1,mi2)));
    temp(rem==1,:)=NaN;
    Cbg(:,ri)=mean(temp,2);
    C(rem==0,ri)=double(data(rem==0,mi1,mi2));
    Cfilt(rem==0,ri)=circshift(filter(ones(filters,1)/filters,1,C(rem==0,ri)),[-floor(filters/2),0]);
    orig{ri,1}=roi{ri,1};

    count=1;
    M(ri,1)=mi1;
    M(ri,2)=mi2;
    fi=surround_1pix(roi{ri,1});
    [sur_d,ix]=ind2sub([d2,d3],fi);
    sur_d(:,2)=ix;
    ai=Cfilt(:,ri);
    chunkcorr=zeros(length(fi),nchunks);
    chunkcorr_bg=zeros(1,nchunks);
    bg=NaN*C(:,1);
    bg(rem==0)=circshift(filter(ones(filters,1)/filters,1,Cbg(rem==0,ri)),[-floor(filters/2),0]);
    
    % Calculate maximum correlation across datachunks for each surround
    % pixel with the center:
    aj=NaN*C(:,1);
    for si=1:length(fi)
        aj(rem==0)=circshift(filter(ones(filters,1)/filters,1,double(squeeze(data(rem==0,sur_d(si,1),sur_d(si,2))))),[-floor(filters/2),0]);
        for sj=1:nchunks
            chunkcorr(si,sj)=corr(ai(int{sj,1}),aj(int{sj,1}));
        end
    end
    corrmat=max(chunkcorr,[],2);
    
    % Correlate maximum correlation with the background:
    for sj=1:nchunks
        chunkcorr_bg(1,sj)=corr(ai(int{sj,1}),bg(int{sj,1}));
    end
    bgcorr=max(chunkcorr_bg);
    
    % Define local threshold for the correlations:
    corrthr_local=max([(max(corrmat)*0.5+min(corrmat)*0.5),0.2,(0.5*max(corrmat)+0.5*bgcorr)]);

    while(~isempty(sur{ri,1})) && ~isempty(intersect(sur{ri,1},sur_d)) 
        new1=sur{ri,1}(1,1);
        new2=sur{ri,1}(1,2);
        
        % add pixels to the ROI if the correlation exceeds the local
        % threshold:
        if 1
            ai=Cfilt(:,ri);
            aj(rem==0)=circshift(filter(ones(filters,1)/filters,1,double(squeeze(data(rem==0,new1,new2)))),[-floor(filters/2),0]);
            chunkcorr=zeros(nchunks,1);
            for sj=1:nchunks
                chunkcorr(sj,1)=corr(ai(int{sj,1}),aj(int{sj,1}));
            end
            cval=max(chunkcorr);
            if cval>=corrthr_local
                bg=mean(data(:,bgt(:,new1,new2)),2);
                temp_bg=Cbg(:,ri)*count+bg;
                temp=C(:,ri)*count+double(squeeze(data(:,new1,new2)));
                count=count+1;
                C(:,ri)=temp/count;
                Cbg(:,ri)=temp_bg/count;
                Cfilt(rem==0,ri)=circshift(filter(ones(filters,1)/filters,1,C(rem==0,ri)),[-floor(filters/2),0]);
                orig{ri,1}(new1,new2)=true;
                roi{ri,1}(new1,new2)=1;
                Ind(new1,new2)=NaN;
                bg(rem==0,:)=circshift(filter(ones(filters,1)/filters,1,Cbg(rem==0,ri)),[-floor(filters/2),0]);
                for sj=1:nchunks
                    chunkcorr_bg(1,sj)=corr(Cfilt(int{sj,1},ri),bg(int{sj,1}));
                end
                bgcorr=max(chunkcorr_bg);
                
                % Update local threshold for the correlations
                corrthr_local=max([corrthr_local,(0.5*max(corrmat)+0.5*bgcorr)]);
                
                % Update the surround:
                for j=-1:1
                    for k=-1:1
                        if new1+j>0 && new1+j<=d2 && new2+k>0 && new2+k<=d3 ...
                                && roi{ri,1}(new1+j,new2+k)==0 ...
                                && ~isnan(Ind(new1+j,new2+k)) ...
                                && sqrt((new1+j-mi1)^2+(new2+k-mi2)^2) < maxradius ...
                                && ~ismember([new1+j,new2+k],sur{ri,1},'rows') 
                            sur{ri,1}=vertcat(sur{ri,1},[new1+j,new2+k]);
                        end
                    end
                end
            end
        end
        sur{ri,1}=sur{ri,1}(2:end,1:2);       
    end
    
    if ri<size(roi,1)
        fprintf(repmat('\b',1,6))
        fprintf('%1.4f',ceil(10^4*(1-my_nansum(Ind(:))/nind))/10^4);
    end
    corrthr(ri)=corrthr_local;
    if count>0 
        ri=ri+1;
    end
    % Remove singular pixels:
    Ind2=ordfilt2(Ind,8,ones(3,3))>0;
    Ind=Ind.*Ind2;
end
nroi=ri-1;

t2=toc(t1);
fprintf('\nTime for ROI initialization: %3.1f minutes\n',t2/60); 
%%
% Remove single-pixel ROIs:
keep=zeros(nroi,1);
for j=1:nroi
    keep(j)=sum(roi{j,1}(:))>1;
end
keep=find(keep);
roi=roi(keep,:);
sur=sur(keep,:);
orig=orig(keep,:);
C=C(:,keep);
Cfilt=Cfilt(1:end-floor(filters/2),keep);
Cbg=Cbg(:,keep);
corrthr=corrthr(keep);




