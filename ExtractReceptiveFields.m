function [onmat,offmat,onmat_shuf,offmat_shuf,on_est,off_est,onoff_est,x_on,y_on,x_off,y_off,x_onoff,y_onoff]=ExtractReceptiveFields(timec,angle,reftime,time_first,time_last,pre,ontime,offtime,filters,nrep)

% Analyzes calcium responses to receptive field sparse noise stimuli (white flashes on black background);
% see Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751 
% __________________________________
% Inputs:
%
% timec: temporal footprint, 1st dim: time, 2nd dim: trials
% angle: visual angle, 1st dim: trials, 2nd dim: space axes
% reftime: start of each imaging frame in ms
% time_first: onset of stimuli
% time_last: offset of stimuli
% pre: interval considered baseline in ms
% ontime: flash length in ms
% offtime: time between flashes in ms
% filters: filter length
% nrep: shuffling repetitions
% __________________________________
% Outputs:
%
% onmat: matrix of ON response
% offmat: matrix of OFF responses
% onmat_shuf: shuffling control matrices for ON responses
% offmat_shuf: shuffling control matrices for OFF responses 
% on_est: 3x3 average filtered ON responses
% off_est: 3x3 average filtered OFF responses 
% onoff_est: sumed on_est+off_est
% x_on,y_on: ON recpective field center
% x_off,y_off: OFF receptive field center
% x_onoff,y_onoff: ON-OFF receptive field center
% __________________________________
%
% copyright: Fiona Muellner, Institute of Molecular and Clinical Ophthalmology Basel, 24.5.2024
% license: BSD (use/copy/modify at own risk), see license.txt

ny=length(unique(angle(:,1,1)));
nx=length(unique(angle(:,2,1)));
vecy=sort(unique(angle(:,1,1)));
vecx=sort(unique(angle(:,2,1)));

mit=min(timec(:));
mat=max(timec(:));
onmat=zeros(nx,ny);
offmat=zeros(nx,ny);

h=mysubplot(ny,nx);
counter=length(time_first);
eps=median(diff(reftime));
F0=prctile(timec(:),10);
filt_timec=filter(ones(filters,1),1,timec);
avg_timec=mean(timec,2);
filt_avg_timec=filter(ones(filters,1),1,avg_timec);
cut=floor(filters/2); % avoid edge effects due to filtering


j0=zeros(length(time_first));
j1=zeros(length(time_first));
j2=zeros(length(time_first));
j3=zeros(length(time_first));
j4=zeros(length(time_first));

for i=1:length(time_first)
    j0(i)=find(reftime>time_first(i)-pre,1,'first');
    j1(i)=find(reftime<time_first(i),1,'last');
    j2(i)=find(reftime>time_last(i)-pre,1,'first');
    j3(i)=find(reftime<time_last(i),1,'last');
    j4(i)=find(reftime<time_last(i)+offtime-eps,1,'last');
end

for i=1:counter
    ix=find(vecx==angle(i,2));
    iy=find(vecy==angle(i,1));
    hold(h(iy+(ix-1)*ny),'all');

    yl=[mit, mat];

    
    line([-pre 0],[1 1]*mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)),'color','g','Linewidth',2,'parent',h(iy+(ix-1)*ny));
    line([-pre 0]+time_last(i)-time_first(i),[1 1]*mean(filt_avg_timec(j2(i):j3(i)-cut)),'color','c','Linewidth',2,'parent',h(iy+(ix-1)*ny));

    if (max(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut))) >= -(min(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))
        onmat(ix,ny-iy+1)=(max(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))/F0*100;
        line([0 ontime],[1 1]*max(filt_avg_timec(j1(i)+1+cut:j3(i)-cut)),'color','r','Linewidth',2,'parent',h(iy+(ix-1)*ny))
    else
        onmat(ix,ny-iy+1)=(min(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))/F0*100;
        line([0 ontime],[1 1]*min(filt_avg_timec(j1(i)+1+cut:j3(i)-cut)),'color','r','Linewidth',2,'parent',h(iy+(ix-1)*ny))
    end
    if (max(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut))) >= -(min(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))
        offmat(ix,ny-iy+1)=(max(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))/F0*100;
        line([0 offtime]+ontime-eps,[1 1]*max(filt_avg_timec(j3(i)+1+cut:j4(i)-cut)),'color','b','Linewidth',2,'parent',h(iy+(ix-1)*ny))
    else
        offmat(ix,ny-iy+1)=(min(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))/F0*100;
        line([0 offtime]+ontime-eps,[1 1]*min(filt_avg_timec(j3(i)+1+cut:j4(i)-cut)),'color','b','Linewidth',2,'parent',h(iy+(ix-1)*ny))
    end

    line([1 1]*0,yl,'color','r','parent',h(iy+(ix-1)*ny))
    line([1 1]*(time_last(i)-time_first(i)),yl,'color','b','parent',h(iy+(ix-1)*ny))
    for j=1:size(timec,2)
        plot(h(iy+(ix-1)*ny),reftime-time_first(i),filt_timec(:,j),'Color',[1 1 1]*0.5);
    end
    plot(h(iy+(ix-1)*ny),reftime-time_first(i),filt_avg_timec,'k-','LineWidth',2);
    xlim(h(iy+(ix-1)*ny),[-pre,ontime+offtime-eps]);
    set(h(iy+(ix-1)*ny),'ylim',yl);
    axis(h(iy+(ix-1)*ny),'off')

end  
line([-pre 0],[1 1]*F0,yl,'color','c','LineWidth',2,'parent',h(1))


timec=cat(1,timec,timec(end,:));
sniplen=(size(timec,1))/counter;
if round(sniplen)~=sniplen
    keyboard
end
ntrial=size(timec,2);
timec=reshape(timec,[sniplen,counter,ntrial]);
shuf_timec=timec;

onmat_shuf=zeros(nx,ny,nrep);
offmat_shuf=zeros(nx,ny,nrep);

for r=1:nrep
    for j=1:ntrial
        s=randperm(counter,counter);
        shuf_timec(:,:,j)=timec(:,s,j);
    end
    shuf_timec=reshape(shuf_timec,[sniplen*counter,ntrial]);

    avg_timec=mean(shuf_timec,2);
    filt_avg_timec=filter(avg_timec,filters);
    
    for i=1:counter
        ix=find(vecx==angle(i,2));
        iy=find(vecy==angle(i,1));
        
        if (max(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut))) >= -(min(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))
            onmat_shuf(ix,ny-iy+1,r)=(max(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))/F0*100;
        else
            onmat_shuf(ix,ny-iy+1,r)=(min(filt_avg_timec(j1(i)+1+cut:j3(i)-cut))-mean(filt_avg_timec(j0(i)+cut:j1(i)-cut)))/F0*100;
        end
        if (max(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut))) >= -(min(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))
            offmat_shuf(ix,ny-iy+1,r)=(max(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))/F0*100;
        else
            offmat_shuf(ix,ny-iy+1,r)=(min(filt_avg_timec(j3(i)+1+cut:j4(i)-cut))-mean(filt_avg_timec(j2(i)+cut:j3(i)-cut)))/F0*100;
        end
    
    end 
    shuf_timec=reshape(shuf_timec,[sniplen,counter,ntrial]);
end

if plotsw

    cvec=zeros(nrep,1);
    
    for rep=1:nrep
        temp1=onmat_shuf(2:end-1,2:end-1,rep);
        temp2=circshift(onmat_shuf(:,:,rep),[1 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=corr(temp1(:),temp2(:));
        temp2=circshift(onmat_shuf(:,:,rep),[-1 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=max(cv,corr(temp1(:),temp2(:)));
        temp2=circshift(onmat_shuf(:,:,rep),[0 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=max(cv,corr(temp1(:),temp2(:)));
        temp2=circshift(onmat_shuf(:,:,rep),[1 0]);
        temp2=temp2(2:end-1,2:end-1);
        cvec(rep)=max(cv,corr(temp1(:),temp2(:)));
    end

    temp1=onmat(2:end-1,2:end-1);
    temp2=circshift(onmat,[1 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=corr(temp1(:),temp2(:));
    temp2=circshift(onmat,[-1 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    temp2=circshift(onmat,[0 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    temp2=circshift(onmat,[1 0]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    pval_on=sum(cvec>cv)/nrep;

    for rep=1:nrep
        temp1=offmat_shuf(2:end-1,2:end-1,rep);
        temp2=circshift(offmat_shuf(:,:,rep),[1 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=corr(temp1(:),temp2(:));
        temp2=circshift(offmat_shuf(:,:,rep),[-1 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=max(cv,corr(temp1(:),temp2(:)));
        temp2=circshift(offmat_shuf(:,:,rep),[0 1]);
        temp2=temp2(2:end-1,2:end-1);
        cv=max(cv,corr(temp1(:),temp2(:)));
        temp2=circshift(offmat_shuf(:,:,rep),[1 0]);
        temp2=temp2(2:end-1,2:end-1);
        cvec(rep)=max(cv,corr(temp1(:),temp2(:)));
    end

    temp1=offmat(2:end-1,2:end-1);
    temp2=circshift(offmat,[1 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=corr(temp1(:),temp2(:));
    temp2=circshift(offmat,[-1 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    temp2=circshift(offmat,[0 1]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    temp2=circshift(offmat,[1 0]);
    temp2=temp2(2:end-1,2:end-1);
    cv=max(cv,corr(temp1(:),temp2(:)));
    pval_off=sum(cvec>cv)/nrep;


    [Y,X]=meshgrid(1:size(onmat,2),1:size(onmat,1));
                            
    figure
    subplot(2,2,3)
    mv=max(max(max(offmat(:)),min(-offmat(:))),max(max(onmat(:)),min(-onmat(:))));
    imshow(offmat',[-mv,mv])
    colormap jet
    title(sprintf('Off, p=%1.3f',pval_off))
    
    off_fit=filter2(ones(3,3)/9,offmat,'same');
    
    subplot(2,2,4)
    hold off
    imshow(off_fit',[-mv,mv])
    colormap jet
    drawnow
    [mv1,mi1]=max(off_fit(:));
    [xi,yi]=ind2sub(size(off_fit),mi1);
    hold all
    
    if pval_off<0.05
        offroi=off_fit'>=0.5*mv1;
        temp=zeros(size(offroi)+2);
        temp(2:end-1,2:end-1)=offroi;
        [roip,proi]=GetROIOutline(temp);
        if iscell(roip)
            for k=1:size(roip,1)
                temp=proi{k,1}(2:end-1,2:end-1);
                roip{k,1}=roip{k,1}-1;
                if temp(yi,xi)
                    plot([roip{k,1}(1:end,1);roip{k,1}(1,1)],[roip{k,1}(1:end,2);roip{k,1}(1,2)],'color','k');
                    off_est=off_fit.*(proi{k,1}(2:end-1,2:end-1)');
                    x_off=sum(X(:).*off_est(:))/sum(off_est(:));
                    y_off=sum(Y(:).*off_est(:))/sum(off_est(:));
                    plot(x_off,y_off,'k*')
                end
            end
        elseif ~isempty(roip)
            roip=roip-1;
            plot([roip(1:end,1);roip(1,1)],[roip(1:end,2);roip(1,2)],'color','k');
            off_est=offmat.*(proi(2:end-1,2:end-1)');
            x_off=sum(X(:).*off_est(:))/sum(off_est(:));
            y_off=sum(Y(:).*off_est(:))/sum(off_est(:));
            plot(x_off,y_off,'k*')
        end
    else
        off_est=NaN*offmat;
        x_off=NaN;
        y_off=NaN;
    end

    subplot(2,2,1)
    imshow(onmat',[-mv,mv])
    colormap jet
    title(sprintf('On, p=%1.3f',pval_on))
    
    on_fit=filter2(ones(3,3)/9,onmat,'same');
    
    subplot(2,2,2)
    hold off
    imshow(on_fit',[-mv,mv])
    colormap jet
    drawnow
    [mv2,mi2]=max(on_fit(:));
    [xi,yi]=ind2sub(size(on_fit),mi2);
    hold all
    
    if pval_on<0.05
        onroi=on_fit'>=0.5*mv2;
        temp=zeros(size(onroi)+2);
        temp(2:end-1,2:end-1)=onroi;
        [roip,proi]=GetROIOutline(temp);
        if iscell(roip)
            for k=1:size(roip,1)
                temp=proi{k,1}(2:end-1,2:end-1);
                roip{k,1}=roip{k,1}-1;
                if temp(yi,xi)
                    plot([roip{k,1}(1:end,1);roip{k,1}(1,1)],[roip{k,1}(1:end,2);roip{k,1}(1,2)],'color','k');
                    on_est=on_fit.*(proi{k,1}(2:end-1,2:end-1)');
                    x_on=sum(X(:).*on_est(:))/sum(on_est(:));
                    y_on=sum(Y(:).*on_est(:))/sum(on_est(:));
                    plot(x_on,y_on,'k*')
                end
            end
        elseif ~isempty(roip)
            roip=roip-1;
            plot([roip(1:end,1);roip(1,1)],[roip(1:end,2);roip(1,2)],'color','k');
            on_est=on_fit.*(proi(2:end-1,2:end-1)');
            x_on=sum(X(:).*on_est(:))/sum(on_est(:));
            y_on=sum(Y(:).*on_est(:))/sum(on_est(:));
            plot(x_on,y_on,'k*')
        end
    else
        on_est=NaN*onmat;
        x_on=NaN;
        y_on=NaN;
    end

    [Y,X]=meshgrid(1:size(on_est,2),1:size(on_est,1));
    onoff_est=on_est+off_est;
    x_onoff=sum(X(:).*onoff_est(:))/sum(onoff_est(:));
    y_onoff=sum(Y(:).*onoff_est(:))/sum(onoff_est(:));
    

    drawnow

end
