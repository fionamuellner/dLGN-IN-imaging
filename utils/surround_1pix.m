function [sur,suroi]=surround_1pix(roi)

% Calculates the 1-pixel surround of a roi mask.
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt

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
