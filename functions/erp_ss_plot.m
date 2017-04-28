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
function deltaAB=erp_ss_plot(f,chan,time,freq,sn,tstruct,L1)

global RUN;
w=3;
yr=[-5 5]; yr2=[-3 3];
Nex=[];

ColorSet=jet(length(RUN.dir.subjects));
Legend=RUN.dir.subjects(setdiff(1:length(RUN.dir.subjects),Nex));

if ~isempty(freq)
    freqOn=1;
    for ii=1:length(f), 
        B{ii}=RUN.ga_freq_Fz.(f{ii});
        timeselB{ii}=find(B{ii}.time>=time(1) & B{ii}.time<=time(2));
        freqselB{ii}=find(B{ii}.freq>=freq(1) & B{ii}.freq<=freq(2));
    end
else
    freqOn=0;
    for ii=1:length(f), 
        A{ii}=RUN.ga.(f{ii}); 
        timeselA{ii}=find(A{ii}.time>=time(1) & A{ii}.time<=time(2));
    end
end


% for ii=1:length(f)
%     c=1;
%     for aa=1:8
% %         if (aa~=4 && aa~=8)
%             for bb=1:62
%                 t=squeeze(A{ii}.individual(aa,bb,:));
%                 tm=conv(t,ones(1,26),'same')./26;
%                 A{ii}.individual(c,bb,:)=tm;
%             end
%             c=c+1;
% %         end
%     end
% end

if ~exist('sn','var'), sn='temp'; end
save1=fullfile(RUN.dir.QAL,[sn '_1.png']);
save2=fullfile(RUN.dir.QAL,[sn '_2.png']);
save3=fullfile(RUN.dir.QAL,[sn '_3.png']);
save4=fullfile(RUN.dir.QAL,[sn '_4.png']);
save5=fullfile(RUN.dir.QAL,[sn '_5.png']);


figure(1);
%=========================================================================%
% Subplot 1
%=========================================================================%
subplot(1,3,1); c=1;
for isub=1:length(L1)
    if L1(isub)==0, continue; end
    for ii=1:length(f)
        switch freqOn
            case 1
                valuesA{ii}(c)=mean(mean(mean(squeeze(B{ii}.powspctrm(c,chan,freqselB{ii},timeselB{ii})))));
            case 0
                valuesA{ii}(c)=mean(mean(squeeze(A{ii}.individual(c,chan,timeselA{ii}))));
        end
        if isempty(freq), N1{ii}(c)=sum(A{ii}.cfg.previous{c}.trials); end
    end
    c=c+1;
end
% ttest of 1 vs 2
deltaAB=[valuesA{1};valuesA{2}]; [h,p,ci,stats]=ttest(valuesA{1}-valuesA{2},0,0.05); display(p)

M=[]; for ii=1:length(f), M=[M,valuesA{ii}']; end

for ii=1:length(M)
   plot(M(ii,:),'o-','linewidth',w,'Color',ColorSet(ii,:)); hold on;
end
xlim([0.75 length(f)+0.25]); set(gca,'xticklabel',[]); legend(Legend); 
title(['p = ' num2str(p)]); hold on;

if isempty(freq),
%=========================================================================%
% Subplot 2
%=========================================================================%
% subplot(1,3,2); c=1;
% for ii=1:sum(RUN.dir.plot)
%    if RUN.dir.plot(isub)==0, continue; end
%    for jj=1:length(f)
%        plot(jj,N1{jj}(c),'o','Color',ColorSet(c,:),'markersize',10,'markerfacecolor',ColorSet(c,:)); hold on;
%    end
%    c=c+1;
% end
% xlim([0.75 length(f)+0.25]); set(gca,'XTickLabel',[])
%=========================================================================%
% Subplot 3
%=========================================================================%
subplot(1,3,3);
N=[];
for ii=1:length(f), N=[N;N1{ii}]; end; N=N';
boxplot(N,'plotstyle','traditional');
set(gcf,'position',[0 0 1280 1024]); export_fig(save1);
%=========================================================================%
end
%=========================================================================%
% Paper Figures...
%=========================================================================%
figure(8); set(gcf,'color','w');
switch freqOn
    case 0
        i1=(A{1}.time>=tstruct.window(1) & A{1}.time<tstruct.window(2));
        b1=squeeze(mean(A{1}.individual(find(L1),chan,:),2));
        b2=squeeze(mean(A{2}.individual(find(L1),chan,:),2));
        t=A{1}.time(i1);
        xl='uV';
    case 1
        i1=(B{1}.time>=tstruct.window(1) & B{1}.time<tstruct.window(2));
        b1=squeeze(mean(mean(B{1}.powspctrm(find(L1),chan,freqselB{1},:),2),3));
        b2=squeeze(mean(mean(B{2}.powspctrm(find(L1),chan,freqselB{2},:),2),3));
        t=B{1}.time(i1);
        xl='z-score(log)';
end
% b3=b1-b2;

plot(t.*1000,mean(b1(:,i1)),'b','linewidth',2); hold on;
plot(t.*1000,mean(b2(:,i1)),'r','linewidth',2); hold on;
% plot(t.*1000,mean(b3(:,i1)),'k','linewidth',2);
grid;
xlabel('Time (ms)','FontSize',14); ylabel(xl,'FontSize',14);
legend(tstruct.legend,'location','SouthEast'); 
set(gca,'FontSize',14);

figure(9); set(gcf,'color','w');
[~,~,~,c1]=ttest(b1,b2);
plot(t,c1.tstat(i1),'k','linewidth',2);
xlabel('Time (s)','FontSize',14); ylabel('T-value','FontSize',14);
set(gca,'FontSize',14);
grid;






%=========================================================================%
% Figure 2) Individual Subjects
%=========================================================================%
figure(2);
cset={'b' 'r' 'g' 'k'}; c=1;
for isub=1:length(RUN.dir.plot)
    if RUN.dir.plot(isub)==0, continue; end
    subplot(ceil(length(RUN.dir.subjects)/3),3,isub);

    for ii=1:length(f)
        switch freqOn
            case 1
                y1=squeeze(mean(mean(B{ii}.powspctrm(c,chan,freqselB{ii},:),2),3));
                x1=B{ii}.time;
            case 0
                y1=squeeze(mean(A{ii}.individual(c,chan,:),2));
                x1=A{ii}.time;
        end
        plot(x1,y1,'linewidth',2,'color',cset{ii}); hold on; 
    end
    
    title(RUN.dir.subjects{isub}); 
    switch freqOn
        case 1, grid;
        case 0, ylim(yr); grid; grid minor;
    end
    c=c+1;
end
set(gcf,'position',[0 0 1280 1024]); export_fig(save2);
%=========================================================================%
% Figure 3) Difference wave & individual waves
%=========================================================================%
% figure(3);
% switch freqOn
%     case 1, tvect=1:length(B{1}.time);
%     case 0, tvect=1:25:length(A{1}.time);
% end
% cc=1;
% for isub=1:length(RUN.dir.subjects),
%     if sum(isub==Nex)>0, continue; end
%     for ii=1:length(f)
%         c=1;
%         for jj=tvect(1:end-1)
%             twin=tvect(c):tvect(c+1)-1;
%             if isempty(freq), tp(c)=mean(A{1}.time(twin)); end
%             switch freqOn
%                 case 1
%                     valuesA{ii}(cc,c)=mean(mean(mean(squeeze(B{ii}.powspctrm(isub,chan,freqselB{ii},twin)))));
%                 case 0
%                     valuesA{ii}(cc,c)=mean(mean(squeeze(A{ii}.individual(isub,chan,twin))));
%             end
%             c=c+1;
%         end
%     end
%     cc=cc+1;
% end
% switch freqOn
%     case 1, tp=diff(B{1}.time)./2+B{1}.time(1:end-1);
%     case 0, tp=round(tp*100)/100;
% end
% for ii=1:length(tp)
%     if mod(ii,2)==0
%         tpl{ii}=num2str(tp(ii));
%     else
%         tpl{ii}=[];
%     end
% end
% subplot(3,1,1);
% notBoxPlot(valuesA{1}); grid; title(f{1});
% set(gca,'xticklabel',tpl);
% switch freqOn, case 0, ylim(yr); end
% subplot(3,1,2);
% notBoxPlot(valuesA{2}); grid; title(f{2});
% set(gca,'xticklabel',tpl);
% switch freqOn, case 0, ylim(yr); end
% 
% subplot(3,1,3);
% notBoxPlot(valuesA{1}-valuesA{2}); grid; title([f{1} ' minus ' f{2}]);
% set(gca,'xticklabel',tpl);
% switch freqOn, case 0, ylim(yr); end
% 
% set(gcf,'position',[0 0 1280 1024]); export_fig(save3);
%=========================================================================%
% Figure 5)
%=========================================================================%
% figure(5);
% m1=mean(valuesA{1}); std1=std(valuesA{1});
% m2=mean(valuesA{2}); std2=std(valuesA{2});
% plot(tp,m1,'-','color','b','linewidth',w); hold on;
% plot(tp,m1+std1,':','color','c','linewidth',w); hold on;
% plot(tp,m1-std1,':','color','c','linewidth',w); hold on;
% plot(tp,m2,'-','color','r','linewidth',w); hold on;
% plot(tp,m2+std2,':','color','k','linewidth',w); hold on;
% plot(tp,m2-std2,':','color','k','linewidth',w); hold on;
% xlabel('Time (s)'); grid; grid minor;
% set(gcf,'position',[0 0 1280 1024]); export_fig(save5);
% figure(6);
% m3=mean(valuesA{3}); std3=std(valuesA{3});
% m4=mean(valuesA{4}); std4=std(valuesA{4});
% plot(tp,m1,'-','color','b','linewidth',w); hold on;
% plot(tp,m2,'-','color','r','linewidth',w); hold on;
% plot(tp,m3,'-','color','g','linewidth',w); hold on;
% plot(tp,m4,'-','color','k','linewidth',w); hold on;
% xlabel('Time (s)'); grid; grid minor;
%=========================================================================%
% Figure 4)
%=========================================================================%
% figure(4);
% subplot(2,1,1);
% [h,~,ci,stats]=ttest(valuesA{1}-valuesA{2});
% % [h,~,ci,stats]=ttest(valuesA{1});
% 
% h(isnan(h))=0;
% plot(tp,nanmean(ci),'g-','linewidth',w); hold on;
% plot(tp,nanmean(ci),'.','markersize',w*8); hold on;
% plot(tp,ci(1,:),'r:','linewidth',w); hold on;
% plot(tp,ci(2,:),'r:','linewidth',w); 
% switch freqOn, case 0, ylim(yr2); end
% ylabel('CI'); xlabel('Time(s)'); grid; grid minor;
% subplot(2,1,2);
% plot(tp,stats.tstat,'linewidth',w); hold on;
% plot(tp(logical(h)),stats.tstat(logical(h)),'r.','markersize',w*8); 
% ylabel('T-value'); xlabel('Time(s)'); grid; grid minor;
% set(gcf,'position',[0 0 1280 1024]); export_fig(save4);
% close all;



