function [corrx,corry,corrv]=rigidmotioncorr(image,options)

%
% Rigid motion correction for calcium imaging data, 
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
%%
% NOTE: this function can be considerably sped up by restricting the
% calculation of the normalized cross correlation to the area of interest.
% To achieve this, e.g. add the following two lines to the Matlab function
% normxcorr2.m at line #84, following the definition of variable 'mn':
%     [ma, na] = size(A);
%     xcorr_TA = xcorr_TA(m:ma,n:na);
% If modified, uncomment the lines following normxcorr2 below.
%%
% __________________________________
% Inputs:
% 
% image: raw data; 1st/2nd dimension space, 3rd dimension time
% options.downs: factor for downsampling in time
% options.ups: factor for upsampling in space
% options.maxoff: maximally allowed offset in x or y, default 30
% __________________________________
% Outputs:
% 
% corrx: shift in 1st dimension
% corry: shift in 2nd dimension
% corrv: correlation coefficient with reference image
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


if nargin<2
    options.downs=1;
    options.ups=1;
    options.maxoff=30;
end
downs=options.downs;
ups=options.ups;
maxoff=options.maxoff;


disp(datetime)

% Optionally, you can use parallel processing to potentially speed up:
% LASTN=maxNumCompThreads(1);
% pool=gcp;

% Without parallel processing:
pool.NumWorkers=1;

d1=size(image,1); 
d2=size(image,2);
lx=size(image,3);

tic
if downs>1
    newimage=zeros(d1,d2,floor(lx/downs));
    for k=downs+1:downs:lx
        newimage(:,:,(k-1)/downs)=mean(image(:,:,k-downs:k),3);
    end
    toc
    image=newimage;
    lx=size(image,3);
    clear newimage
    disp('Downsampled data in time');
end

d1=size(image,1);
d2=size(image,2);


% tile size in x and y
sizey=d1;
sizex=d2;

% reference image size: s1*2
% template image size: (s1-maxoff)*2
s1=min(round((sizey)*0.4),floor(sizey/2)-maxoff);
s2=min(round((sizex)*0.4),floor(sizex/2)-maxoff);
corrx=zeros(lx,1);
corry=zeros(lx,1);
corrv=zeros(lx,1);

firstreg=floor(size(image,3)/2)-24:floor(size(image,3)/2)+25;

% take brightest image region (after filtering) as alignment center 
avg=mean(image(:,:,firstreg),3);
[~,mi]=max(medianfilter2(avg(maxoff+1:end-maxoff,maxoff+1:end-maxoff),5),[],'all','linear');
[c1,c2]=ind2sub(size(avg(maxoff+1:end-maxoff,maxoff+1:end-maxoff)),mi);
c1=c1+maxoff;
c2=c2+maxoff;

%         % take image center as aligment center:
%         c1=round(d1/2);
%         c2=round(d2/2);

% define template borders
sel=[min(max(c1-s1,1)+maxoff,d1-s1*2+maxoff) min(max(c2-s2,1)+maxoff,d2-s2*2+maxoff) (s1-maxoff)*2 (s2-maxoff)*2];
% define reference image borders
cut=[min(max(c1-s1,1),d1-s1*2) min(max(c2-s2,1),d2-s2*2) s1*2 s2*2];

% template image
im4=mean(reshape(image(:,:,firstreg),[d1,d2,50]),3);
im4=im4(sel(1):sel(1)+sel(3),sel(2):sel(2)+sel(4));
if ups>1
    [X,Y] = meshgrid(1:size(im4,2),1:size(im4,1));
    [nX,nY] = meshgrid(1:1/ups:size(im4,2),1:1/ups:size(im4,1));
    im4=interp2(X,Y,im4,nX,nY);
end

% align first 50 images to their overall average
tim1=tic;

image_t=image(:,:,firstreg);
for j=firstreg

   if ups>1
        [X,Y] = meshgrid(1:d2,1:d1);
        [nX,nY] = meshgrid(1:1/ups:d2,1:1/ups:d1);
        image_i=interp2(X,Y,single(image(:,:,j)),nX,nY);
    else
        image_i=image(:,:,j);
    end
    im3=image_i(1+(cut(1)-1)*ups:1+(cut(1)+cut(3)-1)*ups,1+(cut(2)-1)*ups:1+(cut(2)+cut(4)-1)*ups);    
    im3=medianfilter2(im3,2);
    
    c=normxcorr2(im4,im3); % calculates only size(im4):size(im3) 
    c=c(size(im4,1):end-size(im4,1)+1,size(im4,2):end-size(im4,2)+1); % uncomment if normxcorr2 is modified
    [~,mi]=max(c(:));
    [y,x]=ind2sub(size(c),mi);

    % expected correct position in full newimage: 
    % (sel(1)-cut(1))*ups
    y=y-((sel(1)-cut(1))*ups+1);
    x=x-((sel(2)-cut(2))*ups+1);

    corrx(j,1)=-x/ups;
    corry(j,1)=-y/ups;

    % motion correct the first 50 images to obtain new reference;
    % always use interpolation, in order to sharpen the estimate
    image_t(:,:,j-firstreg(1)+1)=circshift(image_t(:,:,j-firstreg(1)+1),[-y,-x]);

end
tim2=toc(tim1);
fprintf('Expected time: %d seconds\n',ceil(tim2/20*lx));

% calculate new reference image
im4=mean(reshape(image_t,[d1,d2,50]),3);
im4=im4(sel(1):sel(1)+sel(3),sel(2):sel(2)+sel(4));

if ups>1
    [X,Y] = meshgrid(1:size(im4,2),1:size(im4,1));
    [nX,nY] = meshgrid(1:1/ups:size(im4,2),1:1/ups:size(im4,1));
    im4=interp2(X,Y,im4,nX,nY);
end

% split image for parallel processing
% length of each partition
stretch=ceil(lx/pool.NumWorkers);
% indices of each partition
parstr=cell(pool.NumWorkers,1);

% partitioned images and correction vectors
% corrnewimage_str=cell(pool.NumWorkers,1);
% newimage_str=cell(pool.NumWorkers,1);
corrv_str=cell(pool.NumWorkers,1);
corrx_str=cell(pool.NumWorkers,1);
corry_str=cell(pool.NumWorkers,1);

for k=1:pool.NumWorkers
    parstr{k,1}=intersect(1:lx,(1:stretch)+(k-1)*stretch);
end
% for k=1:pool.NumWorkers
%     newimage_str{k,1}=image(:,:,parstr{k,1});
% end
for k=1:pool.NumWorkers
%     corrnewimage_str{k,1}=newimage_str{k,1};
    corrv_str{k,1}=zeros(length(parstr{k,1}),1);
    corrx_str{k,1}=zeros(length(parstr{k,1}),1);
    corry_str{k,1}=zeros(length(parstr{k,1}),1);

    for j=parstr{k,1}(1:end)
        jk=j-parstr{k,1}(1)+1;

        if ups>1
            [X,Y] = meshgrid(1:d2,1:d1);
            [nX,nY] = meshgrid(1:1/ups:d2,1:1/ups:d1);
            image_i=interp2(X,Y,single(image(:,:,j)),nX,nY);
        else
            image_i=image(:,:,j);
        end
        im3=image_i(1+(cut(1)-1)*ups:1+(cut(1)+cut(3)-1)*ups,1+(cut(2)-1)*ups:1+(cut(2)+cut(4)-1)*ups);
        im3=medianfilter2(im3,2);
        c=normxcorr2(im4,im3); % calculates only size(im4):size(im3) 
        c=c(size(im4,1):end-size(im4,1)+1,size(im4,2):end-size(im4,2)+1); % uncomment if normxcorr2 is modified
        [mv,mi]=max(c(:));
        [y,x]=ind2sub(size(c),mi);

        % expected correct position in full newimage: 
        % (sel(1)-cut(1))*ups
        y=y-((sel(1)-cut(1))*ups+1);
        x=x-((sel(2)-cut(2))*ups+1);
        
        corrx_str{k,1}(jk,1)=-x/ups;
        corry_str{k,1}(jk,1)=-y/ups;
        corrv_str{k,1}(jk,1)=mv;
        
    end
end

% bring back together images from parallel processing
for k=pool.NumWorkers:-1:1
    corrx(parstr{k,1},1)=corrx_str{k,1};
    corry(parstr{k,1},1)=corry_str{k,1};
    corrv(parstr{k,1},1)=corrv_str{k,1};
end

% clear corrnewimage_str
tim2=toc(tim1);
fprintf('Elapsed time: %d seconds\n',ceil(tim2));


disp(datetime)
