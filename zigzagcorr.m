function [bdshift,corrx1,corrx2]=zigzagcorr(image,zdim,rdim)

% Performs zigzag correction for calcium imaging data;
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% image: raw data
% zdim: dimension containing zigzag artifacts
% rdim: channel used for correction, if image has more than one channel
% __________________________________
% Outputs:
%
% bdshift: average shift per column (to double-check that average across columns is adequate)
% corrx1: shift in forward lines
% corrx2: shift in backward lines
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt


if zdim==1
    m11=mean(image(1:2:end-1,:,1:2:end,rdim),3); % forth
    m21=mean(image(2:2:end,:,1:2:end,rdim),3);   % back
    m12=mean(image(1:2:end-1,:,2:2:end,rdim),3); % back
    m22=mean(image(2:2:end,:,2:2:end,rdim),3);   % forth
    
    cut=0;
    lx=size(image,2);
    x1=zeros(1,lx);
    y1=zeros(1,lx);
    x2=zeros(1,lx);
    y2=zeros(1,lx);
    m1=zeros(1,lx);
    m2=zeros(1,lx);

    int=10;
    maxoff=10;
    for j=1:lx
        % correlate [-int +int] pieces across the back and forth lines:
        [y1(1,j),x1(1,j),m1(1,j)]=imcorr_centr(m11(:,max(1,j-int):min(j+int,lx)),m21(:,max(1,j-int-maxoff):min(j+int+maxoff,lx)),cut,0,3);
        [y2(1,j),x2(1,j),m2(1,j)]=imcorr_centr(m22(:,max(1,j-int):min(j+int,lx)),m12(:,max(1,j-int-maxoff):min(j+int+maxoff,lx)),cut,0,3);
        x1(1,j)=x1(1,j)-max(1,j-int)+max(1,j-int-maxoff);
        x2(1,j)=x2(1,j)-max(1,j-int)+max(1,j-int-maxoff);
    end

    % expect y1=-0.5 and y2=0.5

    f=find(abs(y1(1,:))<=1 & abs(y2(1,:))<=1 & m1>prctile(m1,10) & m2>prctile(m2,10));
    f=setdiff(f,[1,lx]);
    x1(1,:)=interp1([1,f,lx],[x1(f(1)),x1(1,f),x1(f(end))],1:lx);
    x2(1,:)=interp1([1,f,lx],[x2(f(1)),x2(1,f),x2(f(end))],1:lx);
    
    if x1(1,int+maxoff+1)<0
        for j=1:int+maxoff
            x1(1,j)=min(x1(1,j),x1(1,int+maxoff+1));
        end
    end
    if x1(1,end-int-maxoff)>0
        for j=size(x1,2)-int-maxoff+1:size(x1,2)
            x1(1,j)=max(x1(1,j),x1(1,end-int-maxoff));
        end
    end
    
    if x2(1,int+maxoff+1)<0
        for j=1:int+maxoff
            x2(1,j)=min(x2(1,j),x2(1,int+maxoff+1));
        end
    end
    if x2(1,end-int-maxoff)>0
        for j=size(x2,2)-int-maxoff+1:size(x2,2)
            x2(1,j)=max(x2(1,j),x2(1,end-int-maxoff));
        end
    end
    
    shift1=x1;
    shift2=x2;
    shift=(shift1+shift2)/2;
    
    % average both directions, since problem is symmetric, and round to full pixels:
    bdshift=(shift(end:-1:1)+shift)/2;

    % take overall average, because variation is usually within 0.5 pixels,
    % and rounding errors will introduce artifacts if not upsampled.
    % This simplifies correction significantly.
    bdfshift=round(mean(bdshift(50:end-49)));
    corrx1=ceil((bdfshift+0.5)/2);
    corrx2=floor((bdfshift-0.5)/2);

elseif zdim==2
    
    m11=mean(image(:,1:2:end-1,1:2:end,rdim),3); % forth
    m21=mean(image(:,2:2:end,1:2:end,rdim),3);   % back
    m12=mean(image(:,1:2:end-1,2:2:end,rdim),3); % back
    m22=mean(image(:,2:2:end,2:2:end,rdim),3);   % forth
    
    cut=0;
    lx=size(image,1);
    x1=zeros(1,lx);
    y1=zeros(1,lx);
    x2=zeros(1,lx);
    y2=zeros(1,lx);
    m1=zeros(1,lx);
    m2=zeros(1,lx);

    int=10;
    maxoff=10;
    for j=1:lx
        % correlate [-int +int] pieces across the back and forth lines:
        [y1(1,j),x1(1,j),m1(1,j)]=imcorr_centr(m11(max(1,j-int):min(j+int,lx),:),m21(max(1,j-int-maxoff):min(j+int+maxoff,lx),:),cut,0,3,[]);
        [y2(1,j),x2(1,j),m2(1,j)]=imcorr_centr(m22(max(1,j-int):min(j+int,lx),:),m12(max(1,j-int-maxoff):min(j+int+maxoff,lx),:),cut,0,3,[]);
        y1(1,j)=y1(1,j)-max(1,j-int)+max(1,j-int-maxoff);
        y2(1,j)=y2(1,j)-max(1,j-int)+max(1,j-int-maxoff);
    end

    % expect y1=-0.5 and y2=0.5

    f=find(abs(x1(1,:))<=1 & abs(x2(1,:))<=1 & m1>prctile(m1,10) & m2>prctile(m2,10));
    f=setdiff(f,[1,lx]);
    y1(1,:)=interp1([1,f,lx],[y1(f(1)),y1(1,f),y1(f(end))],1:lx);
    y2(1,:)=interp1([1,f,lx],[y2(f(1)),y2(1,f),y2(f(end))],1:lx);
    
    if y1(1,int+maxoff+1)<0
        for j=1:int+maxoff
            y1(1,j)=min(y1(1,j),y1(1,int+maxoff+1));
        end
    end
    if y1(1,end-int-maxoff)>0
        for j=size(y1,2)-int-maxoff+1:size(y1,2)
            y1(1,j)=max(y1(1,j),y1(1,end-int-maxoff));
        end
    end
    
    if y2(1,int+maxoff+1)<0
        for j=1:int+maxoff
            y2(1,j)=min(y2(1,j),y2(1,int+maxoff+1));
        end
    end
    if y2(1,end-int-maxoff)>0
        for j=size(y2,2)-int-maxoff+1:size(y2,2)
            y2(1,j)=max(y2(1,j),y2(1,end-int-maxoff));
        end
    end
    
    shift1=y1;
    shift2=y2;
    shift=(shift1+shift2)/2;
    
    % average both directions, since problem is symmetric, and round to full pixels:
    bdshift=(shift(end:-1:1)+shift)/2;
    
    % take overall average, because variation is usually within 0.5 pixels,
    % and rounding errors will introduce artifacts if not upsampled.
    % This simplifies correction significantly.
    bdfshift=round(mean(bdshift(50:end-49)));
    corrx1=ceil((bdfshift+0.5)/2);
    corrx2=floor((bdfshift-0.5)/2);
    
else
    
    disp('Dimension mismatch')

end