function [roip,proi_new]=GetROIOutline(proi)

% Draws the outline of a region of interest (ROI);
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% proi: logical map of the ROI
% __________________________________
% Outputs:
%
% roip: outline of the ROI as vector of x and y coordinates; 
%       if the ROI contains several unconnected subregions,
%       roip will return a structure with coordinates for each component.
% 
% proi_new: structure containing a separate binary mask for each unconnected region
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


if ~any(proi)
    roip=zeros(0,2);
    proi_new=[];
else
    t1=circshift(proi,[0,1])-proi<0;
    f1=find(t1);
    [y1,x1]=ind2sub(size(proi),f1);
    x1=[x1-0.5;x1-0.5;x1-0.5];
    y1=[y1-0.5;y1;y1+0.5];

    t2=circshift(proi,[0,-1])-proi<0;
    f2=find(t2);
    [y2,x2]=ind2sub(size(proi),f2);
    x2=[x2+0.5;x2+0.5;x2+0.5];
    y2=[y2-0.5;y2;y2+0.5];

    t3=circshift(proi,[1,0])-proi<0;
    f3=find(t3);
    [y3,x3]=ind2sub(size(proi),f3);
    x3=[x3-0.5;x3;x3+0.5];
    y3=[y3-0.5;y3-0.5;y3-0.5];

    t4=circshift(proi,[-1,0])-proi<0;
    f4=find(t4);
    [y4,x4]=ind2sub(size(proi),f4);
    x4=[x4-0.5;x4;x4+0.5];
    y4=[y4+0.5;y4+0.5;y4+0.5];
    
    t5=find(proi(1,:).*proi(end,:));
    if isempty(t5)
        x5=[];
        y5=[];
    else
        x5=[repmat(t5,[2,1]);repmat(t5-0.5,[2,1]);repmat(t5+0.5,[2,1])];
        y5=repmat([0.5*ones(length(t5),1);(size(proi,2)+0.5)*ones(length(t5),1)],[3,1]);
    end
    
    t6=find(proi(:,1).*proi(:,end));
    if isempty(t6)
        x6=[];
        y6=[];
    else
        x6=repmat([0.5*ones(length(t6),1);(size(proi,2)+0.5)*ones(length(t6),1)],[3,1]);
        y6=[repmat(t6,[2,1]);repmat(t6-0.5,[2,1]);repmat(t6+0.5,[2,1])];
    end
    
    x=[x1;x2;x3;x4;x5;x6];
    y=[y1;y2;y3;y4;y5;y6];
    
    temp=[x,y];
    temp=unique(temp,'rows');
    x=temp(:,1);
    y=temp(:,2);
    
    mi=[x(1),y(1)];
    ms=mi;
    roip{1,1}(1,:)=mi;
    count=1;
    countrois=1;
    f=find(sqrt((roip{countrois,1}(count,1)-x).^2+(roip{countrois,1}(count,2)-y).^2)==0.5);
    mp=mi;
    mi=[x(f(1)),y(f(1))];
    count=count+1;
    roip{countrois,1}(count,:)=mi;
    while ~isempty(x) && ~isempty(f)
        f=find(sqrt((roip{countrois,1}(count,1)-x).^2+(roip{countrois,1}(count,2)-y).^2)==0.5);
        if length(f)>2
           if sqrt((mp(1)-x(f(1))).^2+(mp(2)-y(f(1))).^2)==1
                mp=mi;
                mi=[x(f(2)),y(f(2))];
                count=count+1;
                roip{countrois,1}(count,:)=mi;
           else
                mp=mi;
                mi=[x(f(1)),y(f(1))];
                count=count+1;
                roip{countrois,1}(count,:)=mi;
           end
        elseif ~isempty(f)
            if isequal([x(f(1)),y(f(1))],mp)
                mp=mi;
                mi=[x(f(2)),y(f(2))];
                count=count+1;
                roip{countrois,1}(count,:)=mi;
                f=find(x==mp(1)&y==mp(2));
                x(f)=[];
                y(f)=[];
            else
                mp=mi;
                mi=[x(f(1)),y(f(1))];
                count=count+1;
                roip{countrois,1}(count,:)=mi;
                f=find(x==mp(1)&y==mp(2));
                x(f)=[];
                y(f)=[];
            end
        end
        if isequal(mi,ms)
            f=find(x==mi(1)&y==mi(2));
            x(f)=[];
            y(f)=[];
            if ~isempty(x)
                mi=[x(1),y(1)];
                ms=mi;
                countrois=countrois+1;
                count=1;
                roip{countrois,1}(count,:)=mi;
                f=find(sqrt((roip{countrois,1}(count,1)-x).^2+(roip{countrois,1}(count,2)-y).^2)==0.5);
                mp=mi;
                mi=[x(f(1)),y(f(1))];
                count=count+1;
                roip{countrois,1}(count,:)=mi;    
            end
        end
    end
    count=count+1;
    roip{countrois,1}(count,:)=roip{countrois,1}(1,:);

    for c=1:countrois
       proi_new{c,1}=false(size(proi));
    end
    f=find(proi);
    [fi,fj]=ind2sub(size(proi),f);
    for i=1:length(fi)
        for c=1:countrois
            g1=find(roip{c,1}(:,1)==fj(i));
            g2=find(roip{c,1}(:,2)==fi(i));
            if any(roip{c,1}(g1,2)<fi(i)) && any(roip{c,1}(g1,2)>fi(i)) && ...
                any(roip{c,1}(g2,1)<fj(i)) && any(roip{c,1}(g2,1)>fj(i))
                
                proi_new{c,1}(fi(i),fj(i))=true;

            end
        end
    end

    if size(roip,1)==1
        roip=roip{1,1};
        proi_new=proi_new{1,1};
    end
    
end

        
