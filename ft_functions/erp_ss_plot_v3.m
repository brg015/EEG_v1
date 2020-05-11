% B.R. Geib (Winter 2015)
% Function file
%
% deltAB=erp_ss_plot(f1,f2,chan_time)
%
% Inputs:
%   f1 (char)              -> ga.(f1) contrast (on the left)
%   f2 (char)              -> ga.(f2) contrast (on the right)
%   chan [ch1 ch2 ... chN] -> vector of channels to include (idx number)
%   time [start end]       -> time to average over (in seconds)
% Outputs:
%   deltaAB (vector)       -> difference values for each subject
% Description: Plots the difference between two conditions and two time
% points. Also prefroms a t-test on the difference and displays weather the
% result is significant or not.
function SSData=erp_ss_plot_v3(f,chan,tstruct)
%-------------------------------------------------------------------------%
% Preset variables
%-------------------------------------------------------------------------%
if ~isfield(tstruct,'filter'), tstruct.filter.on=0; end
if ~isfield(tstruct,'delta'), tstruct.delta=0; end
if ~isfield(tstruct,'legend'), tstruct.legend=f; end
if ~isfield(tstruct,'fig'), tstruct.fig=true(1,5); end
if ~isfield(tstruct,'SE'), tstruct.SE.on=0; end
if length(tstruct.fig)<5, tstruct(5)=0; end
if tstruct.filter.on==1, cfg=tstruct.filter; end

lim=0.05; seq=3;
tstruct.fig=[1 1 1 1 1]; 

if tstruct.SE.on==1
    if tstruct.SE.diff==0, SErange=1:(length(f)+1); end
    if tstruct.SE.diff==1, SErange=3; end
end
        
tstruct.subj=logical(tstruct.subj);
time=tstruct.time;
%-------------------------------------------------------------------------%
% Setup contrast
%-------------------------------------------------------------------------%
for ii=1:length(f)
    dat{ii}=tstruct.ga.(f{ii}); 
    if tstruct.filter.on==1, 
        dat{ii}=ft_preprocessing(cfg,dat{ii}); 
        dat{ii}.trial=dat{ii}.trial(tstruct.subj,:,:);
        dat{ii}.individual=dat{ii}.trial;
    else
        dat{ii}.individual=dat{ii}.individual(tstruct.subj,:,:);
    end
    % Statistical test are done on these time points
    dat_time{ii}=find(dat{ii}.time>=time(1) & dat{ii}.time<=time(2));
    % Baselining can be done on these
    if isfield(tstruct,'base')
        dat_base{ii}=find(dat{ii}.time>=tstruct.base(1) & dat{ii}.time<=tstruct.base(2));
    end
end
%=========================================================================%
% Baseline
%=========================================================================%
if isfield(tstruct,'base') 
    for ii=1:length(dat)
        for jj=1:size(dat{ii}.individual,1)
            for kk=1:size(dat{ii}.individual,2)
                % Pull time series
                L1=squeeze(dat{ii}.individual(jj,kk,:));
                L2=L1-mean(L1(dat_base{ii}));
                dat{ii}.individual(jj,kk,:)=L2; clear L1 L2;
            end
        end
    end
end
%=========================================================================%
% Difference in extracted values (BOX PLOT)
%=========================================================================%
figure(1);
% Key Variables Here
% dat{#} contrast data
% dat_base{#} potential baseline
% dat_time{#} time to test (T)
for isub=1:sum(tstruct.subj)
    for ii=1:length(f)
        SSData.val(isub,ii)=mean(mean(squeeze(dat{ii}.individual(isub,chan,dat_time{ii}))));
    end
end
% ttest of 1 vs 2
if (tstruct.fig(1)==1 && length(f)==2)
    try
        subplot(2,2,1); set(gcf,'color','w');
        [~,p,~,~]=ttest(SSData.val(:,1)-SSData.val(:,2),0,0.05);
        boxplot(SSData.val,'labels',tstruct.legend(1:2)); title(['p = ' num2str(p)]); grid;
    catch err
    end
    clear p; % val is returned
end

i1=(dat{1}.time>=tstruct.window(1) & dat{1}.time<tstruct.window(2));
for ii=1:length(dat)
    b{ii}=squeeze(mean(dat{ii}.individual(:,chan,i1),2));
end
if tstruct.delta==1, b{3}=b{1}-b{2}; end
t=dat{1}.time(i1);

% Save the key outputs from this
DATA=b; TIME=t;
%=========================================================================%
% Paper Figures...
%=========================================================================%
if tstruct.fig(2)==1
    figure(2);
if length(f)==2
    % subplot(2,2,2);
    cset={'b','r',[0 0 0]};
    wset=[2 2 5];
    lset={'-','-',':'};
else
    cset={'b' 'r' 'g' 'k' 'c' 'm'};
    wset=[5 5 5 5 5 5];
    lset={'-','-','-','-','-','-'};
end

set(gcf,'color','w');

lsize=14;
hwidth=3;

xl='uV';
for ii=1:length(b)
    if sum(tstruct.subj)>1
        plot(t.*1000,mean(b{ii}),'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
    else
        plot(t.*1000,b{ii},'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
    end
end

% Add SE confidence on
if tstruct.SE.on==1
    for ii=SErange
        xconf=[t.*1000 t(end:-1:1)*1000];
        Y=mean(b{ii}); Yerr=std(b{ii})./sqrt(size(b{ii},1));
        yconf=[Y+Yerr Y(end:-1:1)-Yerr(end:-1:1)];
        p=fill(xconf,yconf,cset{ii});
        p.FaceAlpha=0.2;
        p.EdgeColor='none';
        clear xconf Y yconf p;
    end
end

SSData.b=b;
SSData.time=t;

xlabel('Time (ms)','FontSize',lsize); ylabel(xl,'FontSize',lsize);
legend(tstruct.legend,'location','SouthEast'); 
% set(gca,'Ylim',[-2.5 .5]);
set(gca,'FontSize',lsize);
set(gca,'FontWeight','bold');

%---------------------------------------------------------------------%
SEQ =[];
% * for ERP data also wanna subsample
% b{3} is our difference wave to subsample on
fs=0.05; c=1; % sample every 50ms
for ii=TIME(1):fs:TIME(end)
   A(:,c)=mean(b{3}(:,(TIME>=ii & TIME<ii+fs)),2);
   T(c)=ii;
   c=c+1;
end
[~,p,~,~]=ttest(A);

AB=p; 
AB=AB<lim; 
c=0; found=0;
for ii=1:length(AB)
    if (AB(ii)==0 && found==0), continue;
    elseif (AB(ii)==0 && found==1), found=0;
    elseif (AB(ii)==1 && found==0)
        found=1; c=c+1;
        SEQ(c).L=1;
        SEQ(c).I=ii;
    elseif (AB(ii)==1 && found==1)
        SEQ(c).L=SEQ(c).L+1;
        SEQ(c).I=[SEQ(c).I ii];
    end
end
% SEQ(1).I(end)=[];

IND=false(size(AB));
for ii=1:length(SEQ)
    if SEQ(ii).L>=seq, IND(SEQ(ii).I)=1; end
end
% plot(AB,'b'); hold on; plot(IND,'r'); hold on; plot(K,'k');
for ii=1:length(SEQ)
    if SEQ(ii).L>=seq
        a=get(gca,'Ylim'); alow=T(SEQ(ii).I(1))*1000; ahgh=T(SEQ(ii).I(end))*1000; acol=[.9 .9 .9]; 
        hline=patch([alow ahgh ahgh alow alow],[a(1) a(1) a(2) a(2) a(1)],acol,'EdgeColor','none','FaceAlpha',1);
        % hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
    end
    SEQ(ii).T=[T(SEQ(ii).I(1))*1000, T(SEQ(ii).I(end))*1000];
    % Cluster stat too
    [~,p,~,stat]=ttest(mean(A(:,(SEQ(ii).I(1):SEQ(ii).I(end))),2));
    SEQ(ii).p=p;
    SEQ(ii).stat=stat.tstat;
end

SSData.SEQ=SEQ;
%---------------------------------------------------------------------%
    
% a=get(gca,'Ylim'); alow=210; ahgh=310; acol=[.9 .9 .9];
% hline=patch([alow ahgh ahgh alow alow],[a(1) a(1) a(2) a(2) a(1)],acol,'EdgeColor','none','FaceAlpha',1);
% 
% a=get(gca,'Ylim'); alow=400; ahgh=800; acol=[.9 .9 .9];
% hline=patch([alow ahgh ahgh alow alow],[a(1) a(1) a(2) a(2) a(1)],acol,'EdgeColor','none','FaceAlpha',1); 

hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',hwidth);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',hwidth);



for ii=1:length(b)
    if sum(tstruct.subj)>1
        plot(t.*1000,mean(b{ii}),'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
        % Add SE confidence on
        if tstruct.SE.on==1
            if sum(SErange==ii)>0
                xconf=[t.*1000 t(end:-1:1)*1000];
                Y=mean(b{ii}); Yerr=std(b{ii})./sqrt(size(b{ii},1));
                yconf=[Y+Yerr Y(end:-1:1)-Yerr(end:-1:1)];
                p=fill(xconf,yconf,cset{ii});
                p.FaceAlpha=0.2;
                p.EdgeColor='none';
                clear xconf Y yconf p;

            end
        end
    else
        plot(t.*1000,b{ii},'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
    end
end

box off;
set(gca,'TickDir','out'); set(gca,'XLim',[t(1)*1000 t(end)*1000])
legend(tstruct.legend,'location','SouthEast'); grid;
set(gca,'Layer','Top'); 
end
%=========================================================================%
% T figure
%=========================================================================%
  figure(1);
if (tstruct.fig(3)==1 && length(f)==2)

% Smooth data for t-test
subplot(2,2,3); set(gcf,'color','w');
for ii=1:length(dat)
    c=1;
    for aa=1:size(dat{ii}.individual,1)
        for bb=1:size(dat{ii}.individual,2)
            t=squeeze(dat{ii}.individual(aa,bb,:));
            tm=conv(t,ones(1,26),'same')./26;
            dat{ii}.individual(c,bb,:)=tm;
        end
        c=c+1;
    end
end
i1=(dat{1}.time>=tstruct.window(1) & dat{1}.time<tstruct.window(2));
for ii=1:length(dat)
    b{ii}=squeeze(mean(dat{ii}.individual(:,chan,:),2));
end
t=dat{1}.time(i1);

[~,~,~,c1]=ttest(b{1},b{2});
plot(t*1000,c1.tstat(i1),'k','linewidth',2);
set(gca,'XLim',[t(1)*1000 t(end)*1000])

xlabel('Time (s)','FontSize',14); ylabel('T-value','FontSize',14);
set(gca,'FontSize',14); 


end
%=========================================================================%
% T figure
%=========================================================================%
if (tstruct.fig(4)==1 && length(f)==2)

subplot(2,2,4); set(gcf,'color','w');
cset={'k'};
lset={'-'};
wset=[5 5 1]; lsize=14;
hwidth=3;

% i1=(dat{1}.time>=tstruct.window(1) & dat{1}.time<tstruct.window(2));
% for ii=1:length(dat)
%     b{ii}=squeeze(mean(dat{ii}.individual(:,chan,i1),2));
% end
% if tstruct.delta==1, b{3}=b{1}-b{2}; end

t=dat{1}.time(i1);
xl='uV';

ii=1;
if sum(tstruct.subj)>1
    plot(t.*1000,mean(DATA{3}),'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
else
    plot(t.*1000,DATA{3},'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;
end

xlabel('Time (ms)','FontSize',lsize); ylabel(xl,'FontSize',lsize);
% set(gca,'Ylim',[-2.5 .5]);
set(gca,'FontSize',lsize);
hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',hwidth);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',hwidth);
% hline=line([-240 -240],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7],'linestyle',':'); set(hline,'linewidth',3);

ii=1;
plot(t.*1000,mean(DATA{3}),'color',cset{ii},'linewidth',wset(ii),'linestyle',lset{ii}); hold on;


box off;
set(gca,'TickDir','out'); set(gca,'XLim',[t(1)*1000 t(end)*1000])

end

%=========================================================================%
% T figure
%=========================================================================%
% Important variables by here...
% b -> our contrast of interest
if tstruct.fig(5)==1
 
for ii=1:length(DATA)
    figure(2+ii); set(gcf,'color','w');
    imagesc(DATA{ii}); L=get(gca,'CLim');
    if L(1)<0 && L(2)>0, set(gca,'CLim',[-abs(max(L)) max(L)]);
    elseif L(1)<0 && L(2)<0, set(gca,'Clim',[L(1) 0]);
    elseif L(1)>0 && L(2)>0, set(gca,'Clim',[0 L(2)]);
    end
    h1=get(gca,'XTick');
    h1=(h1./250+min(TIME)); set(gca,'XTickLabel',h1); grid;
    colorbar;
    if ii==3, title('delta'); else title(tstruct.legend{ii}); end
    ylabel('SubjectID'); xlabel('Time');
end

end
