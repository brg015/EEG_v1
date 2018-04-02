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

if isfield(cfg,'clim'), CLIM=cfg.clim; else CLIM=[]; end
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
%-------------------------------------------------------------------------%
% Paper figures
%-------------------------------------------------------------------------%
figure(2); set(gcf,'color','w');
cset={'k-','k:','k-.'}; % Same as erp_ss_plot_v3
wset=[4 4 1];

subplot(2,1,1);

if cfg.diff_wave==1
    dat{3}=dat{2}; dat{3}.powspctrm=dat{2}.powspctrm-dat{1}.powspctrm; 
    dat_freq{3}=dat_freq{2};
    dat_time{3}=dat_time{2};
    tstruct.legend{3}='Difference';
end

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
xl=cfg.type;
for ii=1:length(dat)
    plot(t.*1000,b{ii},cset{ii},'linewidth',wset(ii)); hold on;
end
set(gca,'FontSize',16); set(gca,'XLim',[t(1)*1000 t(end)*1000]);
xlabel('Time (ms)','FontSize',16); ylabel(xl,'FontSize',16);
legend(tstruct.legend,'location','SouthEast'); box off;
set(gca,'TickDir','out'); set(gca,'FontWeight','bold');
% set(gca,'Ylim',[.72 .84]);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
if cfg.diff_wave==1, hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3); end
% hline=line([-240 -240],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7],'linestyle','--'); set(hline,'linewidth',3);
for ii=1:length(dat)
    plot(t.*1000,b{ii},cset{ii},'linewidth',wset(ii)); hold on;
end

if length(dat)>1
    subplot(2,1,2);
    X=size(dat{ii}.powspctrm);
    x3=false(1,X(4)); x3(dat_time{ii})=1;
    X1=squeeze(mean(dat{1}.powspctrm(:,cfg.elec,:,x3),2));
    X2=squeeze(mean(dat{2}.powspctrm(:,cfg.elec,:,x3),2));
    [~,p,~,c1]=ttest(X1-X2);
    c1.tstat(p>0.05)=NaN; % c1.tstat(c1.tstat<0)=NaN;
    % imagesc(dat{1}.time(x3),dat{1}.freq,squeeze(c1.tstat)); 
    surf(dat{1}.time(x3)*1000,dat{1}.freq,squeeze(c1.tstat),'LineStyle','none');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30 50])
    set(gca,'YTick',[4 10 20 30 50])
    set(gca,'YLim',[0 30]);
    set(gca,'XLim',[t(1)*1000 t(end)*1000]);
    xlabel('Time (ms)','FontSize',16); ylabel('Freq (Hz)','FontSize',16);
    title('T-Test','FontSize',16); set(gca,'FontSize',16); set(gca,'YDir','normal');
    set(gca,'CLim',[-4 4]); set(gca,'FontWeight','bold');
    colormap('jet'); colorbar;
end
%-------------------------------------------------------------------------%
% Difference
%-------------------------------------------------------------------------%
figure(3); set(gcf,'color','w');
subplot(2,1,1);
for ii=1:length(dat)
    plot(t.*1000,b{1}-b{2},'k-.','linewidth',1); hold on;
end
hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
% hline=line([-240 -240],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7],'linestyle','--'); set(hline,'linewidth',3);
for ii=1:length(dat)
    plot(t.*1000,b{1}-b{2},'k','linewidth',5); hold on;
end
xlabel('Time (ms)','FontSize',16); ylabel([xl 'difference'],'FontSize',16);
set(gca,'FontSize',16); set(gca,'XLim',[t(1)*1000 t(end)*1000]);
box off; set(gca,'TickDir','out'); set(gca,'FontWeight','bold');

subplot(2,1,2);
for ii=1:length(dat)
    % cast to logical indices
    X=size(dat{ii}.powspctrm);
    x1=false(1,X(2)); x1(cfg.elec)=1;
    x2=false(1,X(3)); x2(dat_freq{ii})=1;
    x3=false(1,X(4)); x3(dat_time{ii})=1;
    y1=dat{ii}.powspctrm(:,x1,x2,x3);
    y2{ii}=squeeze(mean(mean(y1,2),3));
    clear X x1 x2 x3 y1 
end

[~,~,~,stat]=ttest(y2{1}-y2{2});
for ii=1:length(dat)
    plot(t.*1000,stat.tstat,'k','linewidth',5); hold on;
end
hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
hline=line([0 0],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
% hline=line([-240 -240],[get(gca,'Ylim')]); set(hline,'color',[.7 .7 .7],'linestyle','--'); set(hline,'linewidth',3);
for ii=1:length(dat)
    plot(t.*1000,stat.tstat,'k','linewidth',5); hold on;
end
xlabel('Time (ms)','FontSize',16); ylabel('T-Test','FontSize',16);
set(gca,'FontSize',16); set(gca,'XLim',[t(1)*1000 t(end)*1000]);
set(gca,'FontWeight','bold'); colormap('jet');
box off;
set(gca,'TickDir','out'); set(gca,'FontWeight','bold');
%-------------------------------------------------------------------------%
% Raw Spectrograms
%-------------------------------------------------------------------------%
% cast to logical indices
X=size(dat{1}.powspctrm);
x1=false(1,X(2)); x1(cfg.elec)=1;
x3=false(1,X(4)); x3(dat_time{ii})=1;
X1=squeeze(mean(dat{1}.powspctrm(:,x1,:,x3),2));
if length(dat)>1, X2=squeeze(mean(dat{2}.powspctrm(:,x1,:,x3),2)); end

figure(4); set(gcf,'color','w');
T=dat{1}.time(x3)*1000;

subplot(2,1,1);
surf(T,ga.(cfg.contrast{1}).freq,squeeze(mean(X1)),'LineStyle','none');
xlabel('Time (ms)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
title(tstruct.legend{1},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
view(0,90); set(gca,'yscale','log'); 
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
set(gca,'YLim',[0 30]);
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
set(gca,'XLim',[T(1) T(end)]);
if ~isempty(CLIM), set(gca,'CLim',CLIM); end
CLim=get(gca,'CLim'); colorbar; % Yoke axes
set(gca,'FontWeight','bold'); colormap('jet');

if length(dat)>1
    subplot(2,1,2);
    surf(T,ga.(cfg.contrast{1}).freq,squeeze(mean(X2)),'LineStyle','none');
    xlabel('Time (s)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
    title(tstruct.legend{2},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'YLim',[0 30]);
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'XLim',[T(1) T(end)]);
    set(gca,'CLim',CLim); colorbar;
    set(gca,'FontWeight','bold'); colormap('jet');

    figure(4); set(gcf,'color','w');
    T=dat{1}.time(x3)*1000;

    figure(5); set(gcf,'color','white')
    surf(T,ga.(cfg.contrast{1}).freq,squeeze(mean(X1)-mean(X2)),'LineStyle','none');
    xlabel('Time (ms)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
    title([tstruct.legend{1} ' minus ' tstruct.legend{2}],'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'YLim',[0 30]);
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'XLim',[T(1) T(end)]);
    if ~isempty(CLIM), set(gca,'CLim',[-.03 .03]); end
    set(gca,'CLim',[-.03 .03]);
    % CLim=get(gca,'CLim'); colorbar; % Yoke axes
    set(gca,'FontWeight','bold'); colormap('jet');

end

end


