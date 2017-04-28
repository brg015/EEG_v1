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
function val=erp_ss_plot_v2(f,chan,time,tstruct)
%-------------------------------------------------------------------------%
% Preset variables
%-------------------------------------------------------------------------%
if ~isfield(tstruct,'filter'), tstruct.filter.on=0; end
if tstruct.filter.on==1, cfg=tstruct.filter; end

tstruct.subj=logical(tstruct.subj);
%-------------------------------------------------------------------------%
% Setup contrast
%-------------------------------------------------------------------------%
for ii=1:length(f), 
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
    if isfield(tstruct,'base')
        dat_base{ii}=find(dat{ii}.time>=tstruct.base(1) & dat{ii}.time<=tstruct.base(2));
    end
end

if isfield(tstruct,'base')
    for isub=1:sum(tstruct.subj)
        for ii=1:length(f)
            val2(isub,ii)=mean(mean(squeeze(dat{ii}.individual(isub,chan,dat_base{ii}))));
        end
        dat{1}.individual(isub,:,:,:)=(dat{1}.individual(isub,:,:,:)-val2(isub,1));
        dat{2}.individual(isub,:,:,:)=(dat{2}.individual(isub,:,:,:)-val2(isub,2));
    end
    clear val2
end
%=========================================================================%
% Difference in extracted values
%=========================================================================%
figure(1); set(gcf,'color','w');
for isub=1:sum(tstruct.subj)
    for ii=1:length(f)
        val(isub,ii)=mean(mean(squeeze(dat{ii}.individual(isub,chan,dat_time{ii}))));
    end
end
% ttest of 1 vs 2
try
    [~,p,~,~]=ttest(val(:,1)-val(:,2),0,0.05);
    boxplot(val,'labels',tstruct.legend); title(['p = ' num2str(p)]); grid;
catch err
end
%=========================================================================%
% Paper Figures...
%=========================================================================%

figure(2); set(gcf,'color','w');
cset={'k-','k:','g-','k-'};

% subplot(2,1,1);
i1=(dat{1}.time>=tstruct.window(1) & dat{1}.time<tstruct.window(2));
for ii=1:length(dat)
    b{ii}=squeeze(mean(dat{ii}.individual(:,chan,:),2));
end
t=dat{1}.time(i1);
xl='uV';

for ii=1:length(dat)
    plot(t.*1000,mean(b{ii}(:,i1)),cset{ii},'linewidth',5); hold on;
end
xlabel('Time (ms)','FontSize',14); ylabel(xl,'FontSize',14);
legend(tstruct.legend,'location','SouthEast'); 
set(gca,'Ylim',[-2.5 .5]);
set(gca,'FontSize',14);
hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',2.5);
hline=line([1 1],[-2.5 .5]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',2.5);
for ii=1:length(dat)
    plot(t.*1000,mean(b{ii}(:,i1)),cset{ii},'linewidth',5); hold on;
end
box off;
set(gca,'TickDir','out');

for ii=1:length(dat)
    c=1;
    for aa=1:size(dat{ii}.individual,1)
        for bb=1:size(dat{ii}.individual,2)
            t=squeeze(dat{ii}.individual(aa,bb,:));
            tm=conv(t,ones(1,26),'same')./26;
%             tm=conv(t,ones(1,13),'same')./13;

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

% subplot(2,1,2);
% [~,~,~,c1]=ttest(b{1},b{2});
% plot(t,c1.tstat(i1),'k','linewidth',2);
% xlabel('Time (s)','FontSize',14); ylabel('T-value','FontSize',14);
% set(gca,'FontSize',14); grid;



