function proj=myprctile2D(data,p,dim)

% Projects data at the p'th percentile along dimension dim
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024

if dim==3
    proj=zeros(size(data,1),size(data,2));
    i=0;
    while i<size(data,1)
       I=i+1:min(i+50,size(data,1));
       proj(I,:)=prctile(data(I,:,:),p,3);
       i=I(end);
    end
elseif dim==2
    proj=zeros(size(data,1),size(data,3));
    i=0;
    while i<size(data,1)
       I=i+1:min(i+50,size(data,1));
       proj(I,:)=prctile(data(I,:,:),p,2);
       i=I(end);
    end
elseif dim==1
    proj=zeros(size(data,2),size(data,3));
    i=0;
    while i<size(data,2)
       I=i+1:min(i+50,size(data,2));
       proj(I,:)=prctile(data(:,I,:),p,1);
       i=I(end);
    end
else
    disp('Dimension not recognized')
    proj=[];
end