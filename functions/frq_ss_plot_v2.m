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
function val=frq_ss_plot_v2(cfg)
%-------------------------------------------------------------------------%
% Preset variables
%-------------------------------------------------------------------------%
% Unpackage cfg a bit
chan=cfg.elec;
time=cfg.time;
freq=cfg.foi;
tstruct=cfg.tstruct;
ga=cfg.ga;
% cfg.layout
%-------------------------------------------------------------------------%
% Setup contrast
%-------------------------------------------------------------------------%
for ii=1:length(cfg.contrast)
    dat{ii}=ga.(cfg.contrast{ii}); 
    dat{ii}.powspctrm=dat{ii}.powspctrm(cfg.L,:,:,:);
    % Statistical test are done on these time points
    dat_test{ii}=find(dat{ii}.time>=cfg.test(1) & dat{ii}.time<=cfg.test(2));
    dat_time{ii}=find(dat{ii}.time>=cfg.time(1) & dat{ii}.time<=cfg.time(2));
    dat_freq{ii}=find(dat{ii}.freq>=cfg.foi(1) & dat{ii}.freq<=cfg.foi(2)); 
end
%=========================================================================%
% Difference in extracted values
%=========================================================================%
figure(1); set(gcf,'color','w');
for isub=1:sum(cfg.L)
    for ii=1:length(cfg.contrast)
        val(isub,ii)=mean(mean(mean(squeeze(dat{ii}.powspctrm(isub,cfg.elec,dat_freq{ii},dat_test{ii})))));
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
cset={'b','r','g','k'};

subplot(2,2,1);
for ii=1:length(dat)
    % cast to logical indices
    X=size(dat{ii}.powspctrm);
    x1=false(1,X(2)); x1(cfg.elec)=1;
    x2=false(1,X(3)); x2(dat_freq{ii})=1;
    x3=false(1,X(4)); x3(dat_time{ii})=1;
    y1=squeeze(mean(dat{ii}.powspctrm,1));
    y2=squeeze(mean(y1(x1,:,:),1));
    y3=squeeze(mean(y2(x2,x3),1));
    b{ii}=y3; clear X x1 x2 x3 y1 y2 y3;
end
t=dat{1}.time(dat_time{1});
xl='Z';

for ii=1:length(dat)
    plot(t.*1000,b{ii},cset{ii},'linewidth',2); hold on;
end
grid;
xlabel('Time (ms)','FontSize',14); ylabel(xl,'FontSize',14);
legend(tstruct.legend,'location','SouthEast'); 
set(gca,'FontSize',14);

subplot(2,2,3);
X1=squeeze(mean(dat{1}.powspctrm(:,cfg.elec,:,:),2));
X2=squeeze(mean(dat{2}.powspctrm(:,cfg.elec,:,:),2));
[~,p,~,c1]=ttest(X1-X2);
c1.tstat(p>0.05)=NaN;
imagesc(dat{1}.time,dat{1}.freq,squeeze(c1.tstat)); 
surf(dat{1}.time,dat{1}.freq,squeeze(c1.tstat),'LineStyle','none');
view(0,90); set(gca,'yscale','log'); 
set(gca,'YTickLabel',[4 10 20 30 50])
set(gca,'YTick',[4 10 20 30 50])
set(gca,'YLim',[0 30]);
set(gca,'Xlim',[dat{1}.time(1) dat{1}.time(end)]);
xlabel('Time (s)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
title('T-Test','FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
% colorbar;
keyboard;
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

CLIM=[-.1 .1];
subplot(2,2,2);
surf(dat{1}.time,dat{1}.freq,squeeze(mean(X1)),'LineStyle','none');
xlabel('Time (s)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
title(tstruct.legend{1},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
view(0,90); set(gca,'yscale','log'); 
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
set(gca,'YLim',[0 30]);
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
% set(gca,'YLim',[4 12]);
set(gca,'XLim',[dat{1}.time(1) dat{1}.time(end)]);
set(gca,'CLim',CLIM); 
% colorbar;

subplot(2,2,4);
surf(dat{2}.time,dat{2}.freq,squeeze(mean(X2)),'LineStyle','none');
xlabel('Time (s)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
title(tstruct.legend{2},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
view(0,90); set(gca,'yscale','log'); 
set(gca,'YTickLabel',[4 10 20 30 50])
set(gca,'YTick',[4 10 20 30 50])
set(gca,'YLim',[0 30]);
set(gca,'YLim',[0 30]);
set(gca,'XLim',[dat{2}.time(1) dat{2}.time(end)]);
set(gca,'CLim',CLIM); 
% colorbar;


