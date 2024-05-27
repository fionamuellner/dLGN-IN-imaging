function [y,x,m,c]=imcorr_centr(im1,im2,cut,showim,ups)

% Performs normalized cross correlation of two images im1 and im2;
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% cut: The lowest cut % of the images is set to 0 (in order to avoid correlations in this range e.g. due to uneven illumination)
% showim: plot results if true
% ups: upsample images by integer factor ups
% __________________________________
% Outputs:
%
% y,x: offset
% To correct offset: im1=circshift(im1,[y,x]);
% m: maximum correlation coefficient
% c: cross correlation matrix
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt


if nargin<5
    ups=1;
end
if nargin<4
    showim=1;
end
if ups>1
    [X1,Y1] = meshgrid(1:size(im1,2),1:size(im1,1));
    [nX1,nY1] = meshgrid(1:1/ups:size(im1,2),1:1/ups:size(im1,1));
    [X2,Y2] = meshgrid(1:size(im2,2),1:size(im2,1));
    [nX2,nY2] = meshgrid(1:1/ups:size(im2,2),1:1/ups:size(im2,1));
    im1=interp2(X1,Y1,single(im1),nX1,nY1);
    im2=interp2(X2,Y2,single(im2),nX2,nY2);
end

im3=im1;
im4=im2;
if length(cut)==2
    p3=prctile(im3(:),cut(1));
    p4=prctile(im4(:),cut(1));
    im3(im3<p3)=p3;
    im4(im4<p4)=p4;
    p3=prctile(im3(:),cut(2));
    p4=prctile(im4(:),cut(2));
    im3(im3>p3)=p3;
    im4(im4>p4)=p4;
elseif cut>0
    p3=prctile(im3(:),cut);
    p4=prctile(im4(:),cut);
    im3(im3<p3)=p3;
    im4(im4<p4)=p4;
end
c=normxcorr2_part(im3,im4,[]);
if showim==1
    figure
    imshow(c,[]);
end
[m,i]=max(c(:));
[y,x]=ind2sub(size(c),i);
y=(y-1)/ups;
x=(x-1)/ups;
