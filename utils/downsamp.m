function newimage = downsamp(Y,downs,type)

% Downsample data in its last (2nd or 3rd) dimension
%
% __________________________________
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024

if nargin<4
    type=class(Y);
end

if strcmp(type,'uint16')
    if length(size(Y))==3
        newimage=zeros(size(Y,1),size(Y,2),floor(size(Y,3)/downs),'uint16');
        for j=downs+1:downs:size(Y,3)+1
            newimage(:,:,(j-1)/downs)=uint16(mean(Y(:,:,j-downs:j-1),3));
        end
    elseif length(size(Y))==2
        newimage=zeros(size(Y,1),floor(size(Y,2)/downs),'uint16');
        for j=downs+1:downs:size(Y,2)+1
            newimage(:,(j-1)/downs)=uint16(mean(Y(:,j-downs:j-1),2));
        end
    end
else
    if length(size(Y))==3
        newimage=zeros(size(Y,1),size(Y,2),floor(size(Y,3)/downs));
        for j=downs+1:downs:size(Y,3)+1
            newimage(:,:,(j-1)/downs)=uint16(mean(Y(:,:,j-downs:j-1),3));
        end
    elseif length(size(Y))==2
        newimage=zeros(size(Y,1),floor(size(Y,2)/downs));
        for j=downs+1:downs:size(Y,2)+1
            newimage(:,(j-1)/downs)=uint16(mean(Y(:,j-downs:j-1),2));
        end
    end
end