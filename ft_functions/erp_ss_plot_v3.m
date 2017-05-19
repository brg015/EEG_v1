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
function val=erp_ss_plot_v3(f,chan,tstruct)
%-------------------------------------------------------------------------%
% Preset variables
%-------------------------------------------------------------------------%
if ~isfield(tstruct,'filter'), tstruct.filter.on=0; end
if tstruct.filter.on==1, cfg=tstruct.filter; end

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
        for jj=1:size(dat{ii}.trial,1)
            for kk=1:size(dat{ii}.trial,2)
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
% Key Variables Here
% dat{#} contrast data
% dat_base{#} potential baseline
% dat_time{#} time to test (T)
for isub=1:sum(tstruct.subj)
    for ii=1:length(f)
        val(isub,ii)=mean(mean(squeeze(dat{ii}.individual(isub,chan,dat_time{ii}))));
    end
end
% ttest of 1 vs 2
if tstruct.fig(1)==1
    try
        figure(1); set(gcf,'color','w');
        [~,p,~,~]=ttest(val(:,1)-val(:,2),0,0.05);
        boxplot(val,'labels',tstruct.legend); title(['p = ' num2str(p)]); grid;
    catch err
    end
    clear p; % val is returned
end

%=========================================================================%
% Paper Figures...
%=========================================================================%
if tstruct.fig(2)==1

figure(2); set(gcf,'color','w');
cset={'k-','k:','k:','k-'};

i1=(dat{1}.time>=tstruct.window(1) & dat{1}.time<tstruct.window(2));
for ii=1:length(dat)
    b{ii}=squeeze(mean(dat{ii}.individual(:,chan,i1),2));
end
t=dat{1}.time(i1);
xl='uV';

for ii=1:length(dat)
    plot(t.*1000,mean(b{ii}),cset{ii},'linewidth',5); hold on;
end

xlabel('Time (ms)','FontSize',16); ylabel(xl,'FontSize',16);
legend(tstruct.legend,'location','SouthEast'); 
% set(gca,'Ylim',[-2.5 .5]);
set(gca,'FontSize',16);
hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
hline=line([-240 -240],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7],'linestyle',':'); set(hline,'linewidth',3);

for ii=1:length(dat)
    plot(t.*1000,mean(b{ii}),cset{ii},'linewidth',5); hold on;
end
box off;
set(gca,'TickDir','out');

end
%=========================================================================%
% T figure
%=========================================================================%
if tstruct.fig(3)==1

figure(3); set(gcf,'color','w');
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
plot(t,c1.tstat(i1),'k','linewidth',2);
xlabel('Time (s)','FontSize',14); ylabel('T-value','FontSize',14);
set(gca,'FontSize',14); grid;

end


