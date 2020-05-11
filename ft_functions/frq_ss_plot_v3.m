% B.R. Geib (Summer 2019)
% Function file
%
% O = freq_ss_plot_v3(cfg)
%
% cfg.L           -> logical array: subjects to include  
% cfg.ga          -> ga structure for ft
% cfg.elec	      -> double array: list of electrodes to include
% cfg.time        -> time in ms to plot [min max]
% cfg.foi         -> frequency of interest band [low high]
% cfg.test        -> temporal test interval [min max]
% cfg.contrast    -> things to contrast (from ga)
% cfg.double_axis -> true or false
% cfg.tstruct.legend  -> legend label
%
% Description: Plots the difference between two conditions and two time
% points. Also prefroms a t-test on the difference and displays weather the
% result is significant or not. Makes paper quality figures :)
function O=frq_ss_plot_v3(cfg)
%-------------------------------------------------------------------------%
% Preset variables
%-------------------------------------------------------------------------%
SEQ = [];
lim=0.05; % critical p-value
seq=3;    % length of sequence

if length(cfg.contrast)>1, DIFF=true; else DIFF=false; end


%-------------------------------------------------------------------------%
% Do required maths for all the later plots
%-------------------------------------------------------------------------%
% 1) Extract data from ga for wanted subjects
for ii=1:length(cfg.contrast)
    dat{ii}=cfg.ga.(cfg.contrast{ii}); 
    dat{ii}.powspctrm=dat{ii}.powspctrm(cfg.L,:,:,:);
end
% 2) Use the first contrast to find indices for test, time, and freq
% -> these are identical, so first is arbit
dat_time=find(dat{1}.time>=cfg.time(1) & dat{1}.time<=cfg.time(2));
dat_freq=find(dat{1}.freq>=cfg.foi(1) & dat{1}.freq<=cfg.foi(2));

X=size(dat{1}.powspctrm); % Again 1 is arbit
elec_ind=false(1,X(2)); elec_ind(cfg.elec)=1; % Elec logical
freq_ind=false(1,X(3)); freq_ind(dat_freq)=1; % Freq logical
time_ind=false(1,X(4)); time_ind(dat_time)=1; % Time logical

O.time=dat{1}.time(dat_time); % Save out the time vector
clear X;

% 3) Do a simple t-test based upon time, and save out the subjectwise
% values as well
if DIFF
    dat_test=find(dat{1}.time>=cfg.test(1) & dat{1}.time<=cfg.test(2));
    for isub=1:sum(cfg.L)
        for ii=1:length(cfg.contrast)
            O.val(isub,ii)=mean(mean(mean(squeeze(dat{ii}.powspctrm(isub,cfg.elec,dat_freq,dat_test)))));
        end
    end
    clear dat_test;
    % Create a difference variable as well
    dat{3}=dat{2}; dat{3}.powspctrm=dat{1}.powspctrm-dat{2}.powspctrm; 
    cfg.tstruct.legend{3}='Difference';   
end
% Output variable: 
% -> O.val(subject,contrast)
% -> dat (updated)
%
% 4) create some time series based upon toi and foi
for ii=1:length(dat)
    y1=squeeze(mean(dat{ii}.powspctrm(:,elec_ind,:,:),2)); % electrodes
    y2=squeeze(mean(y1(:,freq_ind,:),2));                  % freq
    y3=squeeze(mean(y2(:,time_ind),1));                    % time-series
    O.b{ii}=y2(:,time_ind);                                % subject time-series
    b{ii}=y3; 
    clear y1 y2 y3;
end
% Output variable: 
% -> O.b{contrast}(subject,time)
% -> b{contrast}(time) -> average time-series
%
% 5) find sequences of signifcant effects

if DIFF
    % Paired t-test across subjects for contrast 1 and 2
    [~,p,~,~]=ttest(O.b{1}-O.b{2});
    %---------------------------------------------------------------------%
    % the clever little algorithm
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
        
    IND=false(size(AB));
    for ii=1:length(SEQ)
        if SEQ(ii).L>=seq, IND(SEQ(ii).I)=1; end
    end
    % IND (logical array) -> marks significant time-periods
    % This is used later in plotting
    for ii=1:length(SEQ)
        % Temporal interval
        SEQ(ii).T=[O.time(SEQ(ii).I(1))*1000, O.time(SEQ(ii).I(end))*1000];
        % Cluster stat too
        [~,p,~,stat]=ttest(mean(O.b{3}(:,(SEQ(ii).I(1):SEQ(ii).I(end))),2));
        SEQ(ii).p=p;
        SEQ(ii).stat=stat.tstat;
    end
    O.SEQ=SEQ;
end   
% Output variable: 
% -> O.SEQ
% -> IND    -> time periods to mark in plots
%=========================================================================%
% Figure 1: Boxplot
%=========================================================================%
if DIFF
    figure(1); set(gcf,'color','w');
    [~,p,~,~]=ttest(O.val(:,1)-O.val(:,2),0,0.05);
    boxplot(O.val,'labels',[cfg.tstruct.legend(1) cfg.tstruct.legend(2)]); title(['p = ' num2str(p)]); grid;
    set(gca,'fontsize',14);
end
%=========================================================================%
% Figure 2: Raw Spectograms
%=========================================================================%
% X1 and X2 spectograms
% T time variable
X1=squeeze(mean(dat{1}.powspctrm(:,elec_ind,:,time_ind),2));
if DIFF, X2=squeeze(mean(dat{2}.powspctrm(:,elec_ind,:,time_ind),2)); end

figure(2); set(gcf,'color','w');
T=O.time*1000; % cast to ms

if DIFF, subplot(2,1,1); end
surf(T,cfg.ga.(cfg.contrast{1}).freq,squeeze(mean(X1)),'LineStyle','none');
xlabel('Time (ms)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
title(cfg.tstruct.legend{1},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
view(0,90); set(gca,'yscale','log'); 
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
set(gca,'YLim',[0 30]);
set(gca,'YTickLabel',[4 10 20 30])
set(gca,'YTick',[4 10 20 30])
set(gca,'XLim',[T(1) T(end)]);
CLim=get(gca,'CLim'); colorbar; % Yoke axes
set(gca,'FontWeight','bold'); colormap('jet');

if DIFF
    subplot(2,1,2);
    surf(T,cfg.ga.(cfg.contrast{1}).freq,squeeze(mean(X2)),'LineStyle','none');
    xlabel('Time (s)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
    title(cfg.tstruct.legend{2},'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'YLim',[0 30]);
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'XLim',[T(1) T(end)]);
    set(gca,'CLim',CLim); colorbar;
    set(gca,'FontWeight','bold'); colormap('jet');
end
%=========================================================================%
% Figure 3: Difference Spectrograms
%=========================================================================%
if DIFF
    figure(3); set(gcf,'color','w');
    subplot(1,2,1);   
    [~,p,~,c1]=ttest(X1-X2);
    c1.tstat(p>0.05)=NaN; % c1.tstat(c1.tstat<0)=NaN;
    % imagesc(dat{1}.time(x3),dat{1}.freq,squeeze(c1.tstat)); 
    surf(T,dat{1}.freq,squeeze(c1.tstat),'LineStyle','none');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30 50])
    set(gca,'YTick',[4 10 20 30 50])
    set(gca,'YLim',[0 30]);
    set(gca,'XLim',[T(1) T(end)]);
    xlabel('Time (ms)','FontSize',16); ylabel('Freq (Hz)','FontSize',16);
    title('T-Test','FontSize',16); set(gca,'FontSize',16); set(gca,'YDir','normal');
    set(gca,'CLim',[-4 4]); set(gca,'FontWeight','bold');
    colormap('jet'); colorbar;
   
    subplot(1,2,2);
    surf(T,cfg.ga.(cfg.contrast{1}).freq,squeeze(mean(X1)-mean(X2)),'LineStyle','none');
    xlabel('Time (ms)','FontSize',14); ylabel('Freq (Hz)','FontSize',14);
    title([cfg.tstruct.legend{3}],'FontSize',14); set(gca,'FontSize',14); set(gca,'YDir','normal');
    view(0,90); set(gca,'yscale','log'); 
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'YLim',[0 30]);
    set(gca,'YTickLabel',[4 10 20 30])
    set(gca,'YTick',[4 10 20 30])
    set(gca,'XLim',[T(1) T(end)]);
    CLim=get(gca,'CLim'); colorbar; 
    a=max(abs(CLim)); % Ensure symetry
    set(gca,'CLim',[-a a]); clear a;
    set(gca,'FontWeight','bold'); colormap('jet');
end
%=========================================================================%
% Figure 4: Difference traces
%=========================================================================%
cset={'b-','r-','k:'}; % Same as erp_ss_plot_v3
wset=[2 2 5];

% Setup the y-axis limits
if DIFF
    if cfg.double_axis
        y1=min([b{1} b{2}]); 
        y2=max([b{1} b{2}]);
        y3=min(b{3});
        y4=max(b{3});
        % Force symetry into the diff wave
        y5=max([abs(y3),abs(y4)]);
        y3=-y5; y4=y5; clear y5;
    else
        y1=min([b{1} b{2} b{3}]); 
        y2=max([b{1} b{2} b{3}]);
        y3=y1;
        y4=y2;
    end
    y3=y3-(y4-y3)*0.1;
    y4=y4+(y4-y3)*0.1;
else
    y1=min(b{1});
    y2=max(b{1});
end
y1=y1-(y2-y1)*0.1;
y2=y2+(y2-y1)*0.1;
clear p;

h=figure(4); set(gcf,'color','w');
set(h,'defaultAxesColorOrder',[0 0 0; 0 0 0]); 

for ii=1:length(dat)
       
    if ii==1
        % needed to setup x axis
        p{ii}=plot(T,b{ii},cset{ii},'linewidth',wset(ii)); hold on;

        if cfg.double_axis, yyaxis left; end
        for jj=1:length(SEQ)
            if SEQ(jj).L>=seq
                % difference yaxis limits
                a=[y3 y4]; alow=T(SEQ(jj).I(1)); ahgh=T(SEQ(jj).I(end)); acol=[.9 .9 .9]; 
                hline=patch([alow ahgh ahgh alow alow],[a(1) a(1) a(2) a(2) a(1)],acol,'EdgeColor','none','FaceAlpha',1);
                hold on;
            end
        end
        % Line to demarcate 0 ms
        hline=line([0 0],[y3 y4]); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3);
        % Line to demarcate intercept
        hline=refline(0,0); set(hline,'color',[.7 .7 .7]); set(hline,'linewidth',3); 
    end

    if ii<=2
        if cfg.double_axis, yyaxis right; end
        p{ii}=plot(T,b{ii},cset{ii},'linewidth',wset(ii)); hold on;
        ylabel('Memory: Log_1_0(Power)','FontSize',16);
        set(gca,'YLim',[y1 y2]);
        if (cfg.SE.on==1 && cfg.SE.diff~=1)
            xconf=[T T(end:-1:1)];
            Y=mean(O.b{ii}); Yerr=std(O.b{ii})./sqrt(size(O.b{ii},1));
            yconf=[Y+Yerr Y(end:-1:1)-Yerr(end:-1:1)];
            pp=fill(xconf,yconf,cset{ii});
            pp.FaceAlpha=0.2;
            pp.EdgeColor='none';
            clear xconf Y yconf pp;
        end
    end

    if ii==3
        if cfg.double_axis, yyaxis left; end           

        p{ii}=plot(T,b{ii},cset{ii},'linewidth',wset(ii)); hold on;
        
        if (cfg.SE.on==1)
            xconf=[T T(end:-1:1)];
            Y=mean(O.b{ii}); Yerr=std(O.b{ii})./sqrt(size(O.b{ii},1));
            yconf=[Y+Yerr Y(end:-1:1)-Yerr(end:-1:1)];
            pp=fill(xconf,yconf,cset{ii});
            pp.FaceAlpha=0.2;
            pp.EdgeColor='none';
            clear xconf Y yconf pp;
        end
        
        ylabel('Difference: Log_1_0(Power)','FontSize',16);
        set(gca,'YLim',[y3 y4]);
        set(gca,'TickDir','out'); set(gca,'FontWeight','bold');
    end

end

set(gca,'FontSize',16); set(gca,'XLim',[T(1) T(end)]); xlabel('Time (ms)'); 
if DIFF
    legend([p{1},p{2},p{3}],cfg.tstruct.legend,'location','SouthEast'); box off; 
else
    legend([p{1}],cfg.tstruct.legend,'location','SouthEast'); box off; 
end
grid; set(gca,'Layer','Top'); 




