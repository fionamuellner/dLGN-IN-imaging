function c=normxcorr2_part(temp,image,csiz)

% Calculates the accurate (non fft based) normalized cross-correlation for overlapping parts of temp and image.
% csiz: [csiz1 csiz2 csiz3 csiz4] will only calculate the rectangle [csiz1:csiz2,csiz3:csiz4]
% instead of the full [1:size(image,1)-size(temp,1)+1,1:size(image,2)-size(temp,2)+1] matrix;
% can speed up computation if csiz is much smaller than the full matrix.
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt


ty=size(temp,1);
tx=size(temp,2);
    
if nargin<3 || isempty(csiz)
    csiz=[1,size(image,1)-ty+1,1,size(image,2)-tx+1];
end
   
c=zeros(csiz(2)-csiz(1)+1,csiz(4)-csiz(3)+1);
V1=sum(sum((temp-mean(temp(:))).^2));
ntemp=temp-mean(temp(:));
for offy=csiz(1)-1:csiz(2)-1
    for offx=csiz(3)-1:csiz(4)-1
        temp2=image([1:ty]+offy,[1:tx]+offx);
        ntemp2=temp2-mean(temp2(:));
        S=sum(sum(ntemp.*ntemp2));
        V=V1*sum(sum(ntemp2.^2));
        c(offy+1-csiz(1)+1,offx+1-csiz(3)+1)=S/sqrt(V);
    end
end

end