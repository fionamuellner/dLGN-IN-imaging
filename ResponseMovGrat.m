function ResponseMovGrat(data,proi,plotdata,screenrot,time_first,time_last,orientation_val,speed_val,data_frame_start_ms,save_dir,file_name,filters,plotsw,roin,savesw)

% Analyzes calcium responses to moving grating stimuli;
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% data: temporal footprint of the ROI
% proi: binary mask of the ROI
% plotdata: e.g. maximum projection of the full dataset
% screenrot: screen rotation with respect to horizontal axis
% time_first: onset of stimuli
% time_last: offset of stimuli
% orientation_val: stimuli orientation in degree 
% speed_val: stimuli speed 
% data_frame_start_ms: start of each imaging frame in ms
% save_dir: directory to save Results
% file_name: filename to save Results
% filters: filter length
% plotsw: 1 = plot results, 0 = do not plot results 
% roin: ROI number to save Results
% savesw: 1 = save results, 0 = do not save results 
% __________________________________
% Outputs:
%
% saved as 'Results' structure
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


if nargin<12
    filters=10;
end
if nargin<13
    plotsw=1;
end
if nargin<14
    roin=1;
end
if nargin<15
    savesw=1;
end
cd(save_dir);
scrsz=get(0,'ScreenSize');

select=1;

data_frame_start_ms=data_frame_start_ms(end-size(data,1)+1:end); 

unidir=unique(orientation_val);
unidir=mod(unidir-screenrot,360);
univel=setdiff(unique(speed_val),0);

Results.DSI=zeros(length(univel),length(select),2);
Results.OSI=zeros(length(univel),length(select),2);
Results.DSIavg=zeros(length(univel),length(select),2);
Results.OSIavg=zeros(length(univel),length(select),2);
Results.SbC=zeros(length(univel),length(unidir),length(select));
Results.SbCInd=zeros(length(univel),length(select));
Results.RespInd=zeros(length(univel),length(select),2);
Results.PrefDir=zeros(length(univel),length(select),2);
Results.PrefOri=zeros(length(univel),length(select),2);
Results.PrefDiravg=zeros(length(univel),length(select),2);
Results.PrefOriavg=zeros(length(univel),length(select),2);
Results.DirResp=zeros(length(unidir),length(univel),length(select));
Results.AbsDirResp=zeros(length(unidir),length(univel),length(select));
Results.PosPeak=zeros(length(unidir),length(univel),length(select));
Results.NegPeak=zeros(length(unidir),length(univel),length(select));
Results.AvgPeak=zeros(length(unidir),length(univel),length(select));
Results.PosResp=zeros(length(unidir),length(univel),length(select));
Results.NegResp=zeros(length(unidir),length(univel),length(select));
Results.AvgResp=zeros(length(unidir),length(univel),length(select));
Results.PosPeakTime=zeros(length(unidir),length(univel),length(select));
Results.NegPeakTime=zeros(length(unidir),length(univel),length(select));
Results.OffResp=zeros(length(unidir),length(univel),length(select));
Results.Selection=select;
Results.OrigSelection=origsel;
Results.Para.screenrot=screenrot;
Results.Para.time_first=time_first;
Results.Para.time_last=time_last;
Results.Para.orientation_val=orientation_val;
Results.Para.speed_val=speed_val;
Results.Para.data_frame_start_ms=data_frame_start_ms;
Results.Para.filters=filters;
Results.Plotdata=plotdata;
Results.Script=mfilename('fullpath');
Results.F0=zeros(length(select),1);

Results.DSI_shuf=zeros(length(univel),length(select),2,1000);
Results.OSI_shuf=zeros(length(univel),length(select),2,1000);
Results.DSIavg_shuf=zeros(length(univel),length(select),2,1000);
Results.OSIavg_shuf=zeros(length(univel),length(select),2,1000);
Results.SbC_shuf=zeros(length(univel),length(unidir),length(select),1000);
Results.SbCInd_shuf=zeros(length(univel),length(select),1000);
Results.RespInd_shuf=zeros(length(univel),length(select),2,1000);


point_first=time_first*0;
point_last=time_last*0;
point_middle=time_last*0;

dt=median(diff(data_frame_start_ms));

for i=1:length(time_first)
    if speed_val(i)>0
        point_first(i)=find(data_frame_start_ms>=time_first(i),1,'first');
        point_last(i)=find(data_frame_start_ms<=time_last(i),1,'last');
        point_middle(i)=find(data_frame_start_ms>=(time_first(i)+time_last(i))/2,1,'first');
    else 
        point_first(i)=NaN;
    end
end
f=find(~isnan(point_first));
plotint=min(point_first(f(2:end))-point_last(f(1:end-1)));
baselen=round(1500/mean(diff(data_frame_start_ms)));
    
Results.Para.PlotInt=plotint;
Results.Para.BaseLen=baselen;

data_orig=NaN*zeros(round(length(data_frame_start_ms)/4),length(unidir),length(univel),20,length(select));

for ri=1:length(select)
    r=select(ri);
    timec=data(:,r);
    
    if plotsw
        figure
        h1=gcf;

        set(h1,'Position',[scrsz(1:3),scrsz(4)*0.8])
        set(h1,'Visible','off');
    end
    
    fresponse=filter([1],[1],timec);
    yl=[nanmin(fresponse)*0.95,nanmax(fresponse)*1.05];
    
    data=NaN*zeros(round(size(data_frame_start_ms,1)/4),length(unidir),length(univel),20);
    countavg=zeros(length(unidir),length(univel));
    allstimlen=zeros(length(unidir),length(univel),10);
    count=zeros(length(setdiff(univel,0))*(length(unidir)+3),1);
    colors=[0 0.5 1.0; 0 0.75 1.0; 0 0.85 1.0; 0.25 0.85 1.0; 0.5 0.85 1.0; 0.75 0.85 1.0;  ];
    colors=colors./repmat(sum(colors,2),[1,size(colors,2)]);
    colors=min(colors*2,1);
    
    for si=1:length(f)
        i=f(si);
        j=find(unidir==orientation_val(i));
        k=find(univel==speed_val(i));
        countavg(j,k)=countavg(j,k)+1;
        data(plotint-point_first(i)+1+max(point_first(i)-plotint,1):plotint-point_first(i)+1+min(point_last(i)+plotint,length(fresponse)),j,k,countavg(j,k))=...
            fresponse(max(point_first(i)-plotint,1):min(point_last(i)+plotint,length(fresponse)));
        allstimlen(j,k,countavg(j,k))=(point_last(i)-point_first(i));
        if plotsw
            pn=j+2+(length(unidir)+3)*(k-1);
            count(pn)=count(pn)+1;
            ah=subplot(length(setdiff(univel,0)),length(unidir)+3,pn,'parent',h1);
            plot(ah,10^(-3)*data_frame_start_ms(max(point_first(i)-plotint,1):min(point_last(i)+plotint,length(fresponse)))-10^(-3)*data_frame_start_ms(max(point_first(i),1)),...
                fresponse(max(point_first(i)-plotint,1):min(point_last(i)+plotint,length(fresponse))),'color',colors(count(pn),:));
            hold(ah,'all')
            line(10^(-3)*(time_first(i)-data_frame_start_ms(max(point_first(i)-plotint,1)+plotint))*[1 1],yl,'color','r','parent',ah)
            line(10^(-3)*(time_last(i)-data_frame_start_ms(max(point_first(i)-plotint,1)+plotint))*[1 1],yl,'color','b','parent',ah)
            if i>1
                line(10^(-3)*(time_last(i-1)-data_frame_start_ms(max(point_first(i)-plotint,1)+plotint))*[1 1],yl,'color','k','parent',ah)
            end
            if i<length(time_first(f))
                line(10^(-3)*(time_first(i+1)-data_frame_start_ms(max(point_first(i)-plotint,1)+plotint))*[1 1],yl,'color','k','parent',ah)
            end
            xlim(ah,10^(-3)*[data_frame_start_ms(max(point_first(i)-plotint,1))-data_frame_start_ms(max(point_first(i),1))-200,...
                data_frame_start_ms(min(point_last(i)+plotint,length(data_frame_start_ms)))-data_frame_start_ms(max(point_first(i),1))-20]);
            ylim(ah,yl);
        end
    end
    dir_resp=zeros(length(unidir),length(univel));
    abs_dir_resp=zeros(length(unidir),length(univel));
    DSI=zeros(length(univel),2);
    OSI=zeros(length(univel),2);
    DSIavg=zeros(length(univel),2);
    OSIavg=zeros(length(univel),2);
    SbC=zeros(length(univel),length(unidir));
    SbCInd=zeros(length(univel),1);
    RespInd=zeros(length(univel),2);
    PrefDir=zeros(length(univel),2);
    PrefOri=zeros(length(univel),2);
    PrefDiravg=zeros(length(univel),2);
    PrefOriavg=zeros(length(univel),2);
    filtavgdata=NaN*data;
    
    basel=zeros(length(unidir),length(univel));
    stimlen=zeros(length(unidir),length(univel));
    minstimlen=zeros(length(unidir),length(univel));
    maxstimlen=zeros(length(unidir),length(univel));
    stdbasel=zeros(length(unidir),length(univel));
    allbasel=zeros(size(data,2),size(data,3),size(data,4));
    mediandata=zeros(size(data,1),size(data,2),size(data,3));
    avgdata=zeros(size(data,1),size(data,2),size(data,3));
    filtdata=NaN*data;
    
    for k=1:length(univel)
        for j=1:length(unidir)
            for i=1:countavg(j,k)
                allbasel(j,k,i)=prctile(data(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k,i),50);
                filtdata(:,j,k,i)=filter(ones(filters,1),1,data(:,j,k,i));
            end
            mediandata(:,j,k)=nanmedian(data(:,j,k,1:countavg(j,k))-repmat(permute(allbasel(j,k,1:countavg(j,k)),[4,1,2,3]),[size(data,1),1,1,1]),4);
            avgdata(:,j,k)=nanmean(data(:,j,k,1:countavg(j,k)),4);
            stimlen(j,k)=round(mean(allstimlen(j,k,1:countavg(j,k)),3));
            minstimlen(j,k)=floor(min(allstimlen(j,k,1:countavg(j,k)),[],3));
            maxstimlen(j,k)=ceil(max(allstimlen(j,k,1:countavg(j,k)),[],3));
            
            basel(j,k)=prctile(avgdata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k),50);
            stdbasel(j,k)=std(avgdata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k));
            mediandata(:,j,k)=mediandata(:,j,k)-prctile(mediandata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k),50)+basel(j,k);
            
            filtavgdata(:,j,k)=filter(ones(filters,1),1,mediandata(:,j,k,1));
        end
    end
    if any(basel(:)>0)
        basel_global=min(basel(basel(:)>0));
    else 
        basel_global=1;
    end
    pos_peak=zeros(length(unidir),length(univel))*NaN;
    neg_peak=zeros(length(unidir),length(univel))*NaN;
    all_pos_peak=zeros(length(unidir),length(univel),max(countavg(:)))*NaN;
    all_neg_peak=zeros(length(unidir),length(univel),max(countavg(:)))*NaN;
    
    pos_peak_time=zeros(length(unidir),length(univel))*NaN;
    neg_peak_time=zeros(length(unidir),length(univel))*NaN;
    pos_resp=zeros(length(unidir),length(univel))*NaN;
    avg_peak=zeros(length(unidir),length(univel))*NaN;
    neg_resp=zeros(length(unidir),length(univel))*NaN;
    avg_resp=zeros(length(unidir),length(univel))*NaN;
    abs_pos_resp=zeros(length(unidir),length(univel))*NaN;
    abs_neg_resp=zeros(length(unidir),length(univel))*NaN;
    abs_avg_resp=zeros(length(unidir),length(univel))*NaN;
    offresp=zeros(length(unidir),length(univel))*NaN;
    for k=1:length(univel) 
        for j=1:length(unidir)
            SbC(k,j)=-(mean(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k))-basel(j,k))/stdbasel(j,k);
            [mv,mp]=min(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
            neg_peak(j,k)=mv;
            neg_peak_time(j,k)=mp*dt;
            for i=1:countavg(j,k)
                all_neg_peak(j,k,i)=-filtdata(plotint+mp+floor(filters/2),j,k,i)+allbasel(j,k,i);
            end
            [mv,mp]=max(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
            pos_peak(j,k)=mv;
            pos_peak_time(j,k)=mp*dt;
            for i=1:countavg(j,k)
                all_pos_peak(j,k,i)=filtdata(plotint+mp+floor(filters/2),j,k,i)-allbasel(j,k,i);
            end
            offresp(j,k)=(max(filtavgdata(plotint+1+maxstimlen(j,k)+floor(filters/2):plotint+1+maxstimlen(j,k)+2*baselen+floor(filters/2),j,k))-max(filtavgdata(plotint+1+minstimlen(j,k)-baselen-floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k)))/max(filtavgdata(plotint+1+minstimlen(j,k)-baselen-floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
            avg_peak(j,k)=mean(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
        end
        SbCInd(k)=mean(SbC(k,:)>0);
        pos_resp(:,k)=(pos_peak(:,k)-basel(:,k))/basel_global;
        neg_resp(:,k)=-(neg_peak(:,k)-basel(:,k))/basel_global;
        avg_resp(:,k)=(avg_peak(:,k)-basel(:,k))/basel_global;
        if SbCInd(k)>0.5
            dir_resp(:,k)=neg_resp(:,k);
            abs_dir_resp(:,k)=max(dir_resp(:,k),0);
            if plotsw
                for j=1:length(unidir)
                    pn=j+2+(length(unidir)+3)*(k-1);
                    ah=subplot(length(setdiff(univel,0)),length(unidir)+3,pn,'parent',h1);
                    plot(ah,10^(-3)*median(diff(data_frame_start_ms))*([0:length(avgdata(:,j,k))-1]-plotint),filtavgdata(:,j,k),'Linewidth',2,'color','k')
                    line(10^(-3)*median(diff(data_frame_start_ms))*([-baselen 0]-floor(filters/2)),[1 1]*basel(j,k),'color','g','linewidth',2,'parent',ah);
                    line(10^(-3)*median(diff(data_frame_start_ms))*[0+floor(filters/2) stimlen(j,k)-floor(filters/2)],[1 1]*neg_peak(j,k),'color','r','linewidth',2,'parent',ah);
                end
            end
        else
            dir_resp(:,k)=pos_resp(:,k);
            abs_dir_resp(:,k)=max(dir_resp(:,k),0);
            if plotsw
                for j=1:length(unidir)
                    pn=j+2+(length(unidir)+3)*(k-1);
                    ah=subplot(length(setdiff(univel,0)),length(unidir)+3,pn,'parent',h1);
                    plot(ah,10^(-3)*median(diff(data_frame_start_ms))*([0:length(avgdata(:,j,k))-1]-plotint),filtavgdata(:,j,k),'Linewidth',2,'color','k')
                    line(10^(-3)*median(diff(data_frame_start_ms))*([-baselen 0]-floor(filters/2)),[1 1]*basel(j,k),'color','g','linewidth',2,'parent',ah);
                    line(10^(-3)*median(diff(data_frame_start_ms))*[0+floor(filters/2) stimlen(j,k)-floor(filters/2)],[1 1]*pos_peak(j,k),'color','r','linewidth',2,'parent',ah);
                end
            end
        end
        
        abs_pos_resp(:,k)=max(pos_resp(:,k),0);
        [ma,mi]=max(abs_pos_resp(:,k));
        DSI(k,1)=sqrt((sum(abs_pos_resp(:,k).*sind(unidir)))^2+(sum(abs_pos_resp(:,k).*cosd(unidir)))^2)/sum(abs_pos_resp(:,k));
        OSI(k,1)=sqrt((sum(abs_pos_resp(:,k).*sind(2*unidir)))^2+(sum(abs_pos_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_pos_resp(:,k));
        
        RespInd(k,1)=ma*basel_global/stdbasel(mi,k);
        
        abs_avg_resp(:,k)=max(avg_resp(:,k),0);
        DSIavg(k,1)=sqrt((sum(abs_avg_resp(:,k).*sind(unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(unidir)))^2)/sum(abs_avg_resp(:,k));
        OSIavg(k,1)=sqrt((sum(abs_avg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_avg_resp(:,k));
                
        a=(sum(abs_pos_resp(:,k).*sind(unidir))*1i+sum(abs_pos_resp(:,k).*cosd(unidir)))/sum(abs_pos_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefDir(k,1)=acosd(real(a));
        elseif imag(a)<0
            PrefDir(k,1)=(360-acosd(real(a)));
        end
        
        a=(sum(abs_pos_resp(:,k).*sind(unidir*2))*1i+sum(abs_pos_resp(:,k).*cosd(unidir*2)))/sum(abs_pos_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefOri(k,1)=acosd(real(a))/2;
        elseif imag(a)<0
            PrefOri(k,1)=(360-acosd(real(a)))/2;
        end
        
        
        a=(sum(abs_avg_resp(:,k).*sind(unidir))*1i+sum(abs_avg_resp(:,k).*cosd(unidir)))/sum(abs_avg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefDiravg(k,1)=acosd(real(a));
        elseif imag(a)<0
            PrefDiravg(k,1)=(360-acosd(real(a)));
        end
        
        a=(sum(abs_avg_resp(:,k).*sind(unidir*2))*1i+sum(abs_avg_resp(:,k).*cosd(unidir*2)))/sum(abs_avg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefOriavg(k,1)=acosd(real(a))/2;
        elseif imag(a)<0
            PrefOriavg(k,1)=(360-acosd(real(a)))/2;
        end
        
        
        abs_neg_resp(:,k)=max(neg_resp(:,k),0);
        [ma,mi]=max(abs_neg_resp(:,k));
        DSI(k,2)=sqrt((sum(abs_neg_resp(:,k).*sind(unidir)))^2+(sum(abs_neg_resp(:,k).*cosd(unidir)))^2)/sum(abs_neg_resp(:,k));
        OSI(k,2)=sqrt((sum(abs_neg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_neg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_neg_resp(:,k));
        RespInd(k,2)=ma*basel_global/stdbasel(mi,k);
        
        abs_avg_resp(:,k)=max(-avg_resp(:,k),0);
        DSIavg(k,2)=sqrt((sum(abs_avg_resp(:,k).*sind(unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(unidir)))^2)/sum(abs_avg_resp(:,k));
        OSIavg(k,2)=sqrt((sum(abs_avg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_avg_resp(:,k));
        
        
        a=(sum(abs_neg_resp(:,k).*sind(unidir))*1i+sum(abs_neg_resp(:,k).*cosd(unidir)))/sum(abs_neg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefDir(k,2)=acosd(real(a));
        elseif imag(a)<0
            PrefDir(k,2)=(360-acosd(real(a)));
        end
        
        a=(sum(abs_neg_resp(:,k).*sind(unidir*2))*1i+sum(abs_neg_resp(:,k).*cosd(unidir*2)))/sum(abs_neg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefOri(k,2)=acosd(real(a))/2;
        elseif imag(a)<0
            PrefOri(k,2)=(360-acosd(real(a)))/2;
        end
        
         
        a=(sum(abs_avg_resp(:,k).*sind(unidir))*1i+sum(abs_avg_resp(:,k).*cosd(unidir)))/sum(abs_avg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefDiravg(k,2)=acosd(real(a));
        elseif imag(a)<0
            PrefDiravg(k,2)=(360-acosd(real(a)));
        end
        
        a=(sum(abs_avg_resp(:,k).*sind(unidir*2))*1i+sum(abs_avg_resp(:,k).*cosd(unidir*2)))/sum(abs_avg_resp(:,k));
        a=a/abs(a);
        if imag(a)>=0
            PrefOriavg(k,2)=acosd(real(a))/2;
        elseif imag(a)<0
            PrefOriavg(k,2)=(360-acosd(real(a)))/2;
        end
        
        
        if plotsw
            ah=subplot(length(setdiff(univel,0)),length(unidir)+3,3+(length(unidir)+3)*(k-1),'parent',h1);
            title(ah,num2str(univel(k)));
            ah=subplot(length(setdiff(univel,0)),length(unidir)+3,3);
            ah=subplot(length(setdiff(univel,0)),length(unidir)+3,(length(unidir)+3)*(k-1)+length(unidir)+3,'parent',h1);
            px=cosd([unidir;unidir(1)]).*[abs_pos_resp(:,k);abs_pos_resp(1,k)];
            py=sind([unidir;unidir(1)]).*[abs_pos_resp(:,k);abs_pos_resp(1,k)];
            plot(ah,px,py,'-k','Linewidth',2)
            sc1=max([abs(px);abs(py)]);
            hold(ah,'on')
            px=cosd([unidir;unidir(1)]).*[abs_neg_resp(:,k);abs_neg_resp(1,k)];
            py=sind([unidir;unidir(1)]).*[abs_neg_resp(:,k);abs_neg_resp(1,k)];
            plot(ah,px,py,'-','color',[1 1 1]*0.8,'Linewidth',2)
            sc2=max([abs(px);abs(py)]);
            sc=max(sc1,sc2);
            if sc==0
                sc=1;
            end                
            line([-1,1]*sc,[0,0],'color','k','parent',ah');
            line([0,0],[-1,1]*sc,'color','k','parent',ah');
            pos=get(ah,'position');
            set(ah,'position',[pos(1) pos(2)-pos(4)/3 pos(3)*1.3 pos(4)*1.3]);
            title(ah,sprintf('DSI: %1.2f , %1.2f \nOSI: %1.2f , %1.2f \nResp. Ind.: %1.2f , %1.2f \nSbC: %1.1f',...
                DSI(k,1),DSI(k,2),OSI(k,1),OSI(k,2),RespInd(k,1),RespInd(k,2),SbCInd(k)));
            set(ah,'fontsize',8)
            set(ah,'xlim',[-1 1]*sc,'ylim',[-1 1]*sc)
            axis square
        end
    end
    
    roip=GetROIOutline(proi{r,1});
            
    if plotsw
        ah=subplot(length(setdiff(univel,0)),length(unidir)+3,[1,2,length(unidir)+4,length(unidir)+5],'parent',h1);
        if size(plotdata,3)==1
            imshow(plotdata,[prctile(plotdata(:),5),prctile(plotdata(:),99.5)],'parent',ah);
        else
            imshow(plotdata,[],'parent',ah);
        end
        hold(ah,'all')

        if iscell(roip)
            for k=1:size(roip,1)
                plot(ah,[roip{k,1}(1:end,1);roip{k,1}(1,1)],[roip{k,1}(1:end,2);roip{k,1}(1,2)],'color',[0 0 1],'linewidth',2);
            end
        else
            plot(ah,[roip(1:end,1);roip(1,1)],[roip(1:end,2);roip(1,2)],'color',[0 0 1],'linewidth',2);
        end
        pos=get(gca,'position');
        set(gca,'position',[0.01 0   pos(3)*1.8   pos(4)*1.8])
        if savesw
            saveas(h1,sprintf('%s_roi%d',file_name,roin),'png')
        end

    end
    
    Results.DSI(:,ri,:)=DSI;
    Results.OSI(:,ri,:)=OSI;
    Results.DSIavg(:,ri,:)=DSIavg;
    Results.OSIavg(:,ri,:)=OSIavg;
    Results.SbC(:,:,ri)=SbC;
    Results.SbCInd(:,ri)=SbCInd;
    Results.RespInd(:,ri,:)=RespInd;
    Results.PrefDir(:,ri,:)=PrefDir;
    Results.PrefOri(:,ri,:)=PrefOri;
    Results.PrefDiravg(:,ri,:)=PrefDiravg;
    Results.PrefOriavg(:,ri,:)=PrefOriavg;
    Results.DirResp(:,:,ri)=dir_resp;
    Results.AllPosPeak(:,:,:,ri)=all_pos_peak;
    Results.AllNegPeak(:,:,:,ri)=all_neg_peak;
    Results.PosPeakTime(:,:,ri)=pos_peak_time;
    Results.NegPeakTime(:,:,ri)=neg_peak_time;
    Results.PosPeak(:,:,ri)=pos_peak;
    Results.NegPeak(:,:,ri)=neg_peak;
    Results.AvgPeak(:,:,ri)=avg_peak;
    Results.PosResp(:,:,ri)=pos_resp;
    Results.NegResp(:,:,ri)=neg_resp;
    Results.AvgResp(:,:,ri)=avg_resp;
    Results.OffResp(:,:,ri)=offresp;
    Results.Basel(:,:,ri)=basel;
    Results.StdBasel(:,:,ri)=stdbasel;
    Results.F0(ri)=basel_global;
    Results.AbsDirResp(:,:,ri)=abs_dir_resp;
    Results.roi{ri,1}=proi{r,1};
    Results.Para.AllStimLen=allstimlen;
    Results.Para.StimLen=stimlen;
    for k=1:length(univel)
        for j=1:length(unidir)
            Results.AvgData{k,1}(:,j,ri)=avgdata(plotint-baselen:max(stimlen(:,k))+plotint,j,k);
            Results.MedianData{k,1}(:,j,ri)=mediandata(plotint-baselen:max(stimlen(:,k))+plotint,j,k);
        end
    end
    data_orig(:,:,:,:,ri)=data;
end

data_orig=data_orig(:,:,:,1:max(countavg(:)),:);


% Repeat for shuffled data:
numtrials=sum(max(countavg,[],2));
ref_mat=zeros(numtrials,2);
count=0;
for j=1:length(unidir)
    for i=1:max(countavg(j,:))
        count=count+1;
        ref_mat(count,1)=j;
        ref_mat(count,2)=i;
    end
end
for ri=1:length(select)
    for rep=1:1000
        shuf=randperm(numtrials);
        for k=1:length(univel)
            count=0;
            for j=1:length(unidir)
                for i=1:countavg(j,k)
                    count=count+1;
                    data(:,j,k,i)=data_orig(:,ref_mat(shuf(count),1),k,ref_mat(shuf(count),2),ri);
                end
            end
        end
        for k=1:length(univel)
            for j=1:length(unidir)

                for i=1:countavg(j,k)
                    allbasel(j,k,i)=prctile(data(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k,i),50);
                end
                mediandata(:,j,k)=nanmedian(data(:,j,k,1:countavg(j,k))-repmat(permute(allbasel(j,k,1:countavg(j,k)),[4,1,2,3]),[size(data,1),1,1,1]),4);
                avgdata(:,j,k)=nanmean(data(:,j,k,1:countavg(j,k)),4);

                basel(j,k)=prctile(avgdata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k),50);
                stdbasel(j,k)=std(avgdata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k));
                mediandata(:,j,k)=mediandata(:,j,k)-prctile(mediandata(plotint-baselen-floor(filters/2):plotint-floor(filters/2),j,k),50)+basel(j,k);

                filtavgdata(:,j,k)=filter(ones(filters,1),1,mediandata(:,j,k,1));

                SbC(k,j)=-(mean(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k))-basel(j,k))/stdbasel(j,k);
                [mv,mp]=min(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
                neg_peak(j,k)=mv;
                neg_peak_time(j,k)=mp*dt;
                [mv,mp]=max(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
                pos_peak(j,k)=mv;
                pos_peak_time(j,k)=mp*dt;
                offresp(j,k)=(max(filtavgdata(plotint+1+maxstimlen(j,k)+floor(filters/2):plotint+1+maxstimlen(j,k)+2*baselen+floor(filters/2),j,k))-max(filtavgdata(plotint+1+minstimlen(j,k)-baselen-floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k)))/max(filtavgdata(plotint+1+minstimlen(j,k)-baselen-floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
                avg_peak(j,k)=mean(filtavgdata(plotint+1+floor(filters/2):plotint+1+minstimlen(j,k)-floor(filters/2),j,k));
            end

            SbCInd(k)=mean(SbC(k,:)>0);
            pos_resp(:,k)=(pos_peak(:,k)-basel(:,k))/basel_global;
            neg_resp(:,k)=-(neg_peak(:,k)-basel(:,k))/basel_global;
            avg_resp(:,k)=(avg_peak(:,k)-basel(:,k))/basel_global;
            if SbCInd(k)>0.5
                dir_resp(:,k)=neg_resp(:,k);
                abs_dir_resp(:,k)=max(dir_resp(:,k),0);
            else
                dir_resp(:,k)=pos_resp(:,k);
                abs_dir_resp(:,k)=max(dir_resp(:,k),0);
            end

            abs_pos_resp(:,k)=max(pos_resp(:,k),0);
            [ma,mi]=max(abs_pos_resp(:,k));
            DSI(k,1)=sqrt((sum(abs_pos_resp(:,k).*sind(unidir)))^2+(sum(abs_pos_resp(:,k).*cosd(unidir)))^2)/sum(abs_pos_resp(:,k));
            OSI(k,1)=sqrt((sum(abs_pos_resp(:,k).*sind(2*unidir)))^2+(sum(abs_pos_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_pos_resp(:,k));
            
            RespInd(k,1)=ma*basel_global/stdbasel(mi,k);

            abs_avg_resp(:,k)=max(avg_resp(:,k),0);
            DSIavg(k,1)=sqrt((sum(abs_avg_resp(:,k).*sind(unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(unidir)))^2)/sum(abs_avg_resp(:,k));
            OSIavg(k,1)=sqrt((sum(abs_avg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_avg_resp(:,k));
            

            abs_neg_resp(:,k)=max(neg_resp(:,k),0);
            [ma,mi]=max(abs_neg_resp(:,k));
            DSI(k,2)=sqrt((sum(abs_neg_resp(:,k).*sind(unidir)))^2+(sum(abs_neg_resp(:,k).*cosd(unidir)))^2)/sum(abs_neg_resp(:,k));
            OSI(k,2)=sqrt((sum(abs_neg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_neg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_neg_resp(:,k));
            
            abs_avg_resp(:,k)=max(-avg_resp(:,k),0);
            DSIavg(k,2)=sqrt((sum(abs_avg_resp(:,k).*sind(unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(unidir)))^2)/sum(abs_avg_resp(:,k));
            OSIavg(k,2)=sqrt((sum(abs_avg_resp(:,k).*sind(2*unidir)))^2+(sum(abs_avg_resp(:,k).*cosd(2*unidir)))^2)/sum(abs_avg_resp(:,k));
            
            RespInd(k,2)=ma*basel_global/stdbasel(mi,k);


            Results.DSI_shuf(:,ri,:,rep)=DSI;
            Results.OSI_shuf(:,ri,:,rep)=OSI;
            Results.DSIavg_shuf(:,ri,:,rep)=DSIavg;
            Results.OSIavg_shuf(:,ri,:,rep)=OSIavg;
            Results.SbCInd_shuf(:,ri,rep)=SbCInd;
            Results.RespInd_shuf(:,ri,:,rep)=RespInd;
        end
    end
end
if savesw
    save(sprintf('%s_results_roi%d',file_name,roin),'Results')
end
    
