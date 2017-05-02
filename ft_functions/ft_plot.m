% B.R. Geib (Winter 2015)
% Function file
%
% ft_plot(ctrast,ID,scale,save_as,interactive)
%
% Inputs:
%   ctrast{}        -> cell array of contrast to plot, up to a maximum of
%                      four. Plotted in order: 1=blue, 2=red, 3=green,
%                      4=black. Keep this at 2 for topo-maps. Will plot
%                      ctrast{1} - ctrast{2}
%   ID              -> Subject to plot (0 = all subjects)
%   scale           -> scale of ERP maps
%   save_as         -> save name, by default all images are saved to
%                      RUN.dir.QA
%   interactive     -> interactive contrast maps, if turned off, topos are
%                      made and saved. by defualt, topos are set to be
%                      between [-2.5uV and 2.5uV]
% Outputs:
%   None
% Description: Plots the difference between two conditions and two time
% points. Also prefroms a t-test on the difference and displays weather the
% result is significant or not.


function ft_plot(ctrast,ID,scale,save_as,interactive)

global RUN;
cfg_topo = []; cfg_erp = [];

Ga=RUN.ga;
filt=0;
bs=RUN.pre.baseline;
fs='';
template=RUN.template2;
freq=0;

f1=[]; f2=[]; extra=[]; extra2=[];
for ii=1:length(ctrast)
    if ii==1, 
        f1=ctrast{1}; 
        if ~isfield(Ga.(f1),'avg'), Ga.(f1).avg=squeeze(mean(Ga.(f1).individual,1)); end
    end
    if ii==2, 
        f2=ctrast{2}; 
        if ~isfield(Ga.(f2),'avg'), Ga.(f2).avg=squeeze(mean(Ga.(f2).individual,1)); end
    end
    if ii==3, 
        extra=ctrast{3}; 
        if ~isfield(Ga.(f2),'avg'), Ga.(f2).avg=squeeze(mean(Ga.(f2).individual,1)); end
    end
    if ii==4, 
        extra2=ctrast{4}; 
        if ~isfield(Ga.(f2),'avg'), Ga.(f2).avg=squeeze(mean(Ga.(f2).individual,1)); end
    end
end

if isempty(f2), f2=f1; end
close all; 
if ~exist('scale','var'), scale=[-4 4]; end

cfg_topo.parameter = 'individual';
cfg_topo.parameter = 'avg';

if isfield(Ga.(f1),'beta'), 
    cfg_erp.parameter='beta'; 
    cfg_topo.parameter='beta';
end
if ~isfield(Ga.(f1),'dimord'),Ga.(f1).dimord='chan_time'; end
if ~isfield(Ga.(f2),'dimord'),Ga.(f2).dimord='chan_time'; end



%-------------------------------------------------------------------------%
% Set Parameters
%-------------------------------------------------------------------------%
if interactive==1
    group_topo=0;
else
    group_topo=1;
end

subj_topo=0;
overwrite=1;

ts=-.3; tstep=0.025; te=0.6;

% Plot topos
if ~strcmp(f1,f2)
    cfg_topo.zlim=[-1 1];
    if freq==1, cfg_topo.zlim=[-1 1]; end
else
    cfg_topo.zlim=scale;
end
cfg_topo.commentpos='title';
cfg_topo.comment='no';

cfg_topo.layout = template{2}; 
% cfg_topo.interplimits='electrodes';

cfg_erp.layout = template{1}; 
cfg_erp.showlabels  = 'yes'; 
if ~isempty(bs), cfg_erp.baseline = bs; end
cfg_erp.ylim=scale;
cfg_erp.linewidth=2;

cfg_erp2.layout = template{1}; 
cfg_erp2.showlabels  = 'yes'; 
if ~isempty(bs), cfg_erp2.baseline = bs; end
cfg_erp2.ylim=scale;
cfg_erp2.linewidth=2;
% cfg_erp.xlim=[-0.3 0.3];


%-------------------------------------------------------------------------%
if ID~=0
    ID_str=RUN.dir.inc_subjects{ID};
    Ga.(f1).individual(setdiff(1:length(RUN.dir.inc_subjects),ID),:,:)=[];
    if ~strcmp(f1,f2)
        Ga.(f2).individual(setdiff(1:length(RUN.dir.inc_subjects),ID),:,:)=[];
    end
    save_dir=fullfile(RUN.dir.QAL,'Subject',ID_str,'ERP');
else
    ID_str='Group';
%     save_dir=fullfile(RUN.dir.QAL,ID_str,'ERP');
    save_dir=RUN.dir.sav;
end


sdisp(save_as,2);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

save1=fullfile(save_dir,[save_as '_' fs '_' ID_str '_both.png']);
save2=fullfile(save_dir,[save_as '_' fs '_' ID_str '_diff.png']); 
%-------------------------------------------------------------------------%
% Filter Data if asked
%-------------------------------------------------------------------------%
filt=0;
if filt==1
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 10;
    plot1 = ft_preprocessing(cfg,Ga.(f1)); 
    plot1.individual=plot1.trial;
    plot2 = ft_preprocessing(cfg,Ga.(f2));
    plot2.individual=plot2.trial;
else
    plot1 = Ga.(f1);
    plot2 = Ga.(f2);
    if ~isempty(extra)
        plot4 = Ga.(extra);
    end
end
%-------------------------------------------------------------------------%
% Setup diff plot
%-------------------------------------------------------------------------%
plot3=plot1; % Copy params & subtract if needed
plot3.avg=plot1.avg-plot2.avg;
plot3.individual=plot1.individual-plot2.individual;

% [~,~,~,a]=ttest(plot1.avg-plot2.avg);
% plot3.avg=a.tstat;
% [~,~,~,a]=ttest(plot1.individual-plot2.individual);
% plot3.individual=a.tstat;

if isempty(extra2)
    if isempty(extra)
        if (~exist(save1,'file') || overwrite==1)
            figure(1); 
            set(gcf,'position',[0 0 1280 1024]);
            ft_multiplotER(cfg_erp,plot1,plot2);  
            if interactive==0, close(1); end
        end
    else
        if (~exist(save1,'file') || overwrite==1)
            figure(1); 
            set(gcf,'position',[0 0 1280 1024]);
            ft_multiplotER(cfg_erp,plot1,plot2,plot4);  
            if interactive==0, close(1); end
        end
        return; % If extra var included, diff plot never reached
    end
else
    plot5=Ga.(extra2);
    if (~exist(save1,'file') || overwrite==1)
        figure(1); 
        set(gcf,'position',[0 0 1280 1024]);
        ft_multiplotER(cfg_erp,plot1,plot2,plot4,plot5);  
        % export_fig(save1); 
        if interactive==0, close(1); end
    end
    return;
end

if strcmp(f1,f2)~=1
    cfg_erp2.ylim=[-4 4];
    figure(2); ft_multiplotER(cfg_erp2,plot3);
    set(gcf,'position',[0 0 1280 1024]);
    if interactive==0, 
        close(2); 
    else
        return;
    end
end

if interactive==0
    if (ID==0 && subj_topo==1)
        if (~exist(fullfile(save_dir,[save_as '_' fs '_group_topo' num2str(1) '.png']),'file') || overwrite==1)
            % Special now
            s1=plot3;
            s1.individual=s1.individual(1:4,:,:);
            eeg_plot_group_erp(s1,f1,f2,fs,ts,tstep,te,RUN.dir.subjects(1:4),1,save_dir,cfg_topo);

            s2=plot3;
            s2.individual=s2.individual(5:8,:,:);       
            eeg_plot_group_erp(s2,f1,f2,fs,ts,tstep,te,RUN.dir.subjects(5:8),2,save_dir,cfg_topo);   
            
            s2=plot3;
            s2.individual=s2.individual(9:12,:,:);       
            eeg_plot_group_erp(s2,f1,f2,fs,ts,tstep,te,RUN.dir.subjects(9:12),3,save_dir,cfg_topo); 
        end
    end
end

if freq==1
    K=~isnan(plot3.avg(1,:));
    plot3.avg=plot3.avg(:,K);
    plot3.individual=plot3.individual(:,:,K);
    plot3.time=plot3.time(K);
end

% Group over time
% BADc=[5:7 12:14 19:21];
BADc=[6:9 15:18 24:27 33:36];
if group_topo==1
    if (~exist(fullfile(save_dir,[save_as '_' fs '_' ID_str '_topo' num2str(1) '.png']),'file') || overwrite==1)
        c=1; f=1;
        for ii=ts:tstep:te
            cfg_topo.xlim=[ii ii+tstep]; 
            if c==1, 
                save3=fullfile(save_dir,[save_as '_' fs '_' ID_str '_topo' num2str(f) '.png']);
                figure(f); 
                set(gcf,'position',[0 0 1280 1024]); 
                set(gcf,'color','w');
                subaxis(4,7,c, 'Spacing', 0, 'Padding', 0, 'Margin', 0,'SpacingVert',0.05,'MarginTop',0.2);
                ft_topoplotER(cfg_topo, plot3); 
            end
            subaxis(5,9,c, 'Spacing', 0, 'Padding', 0, 'Margin', 0,'SpacingVert',0.05,'MarginTop',0.2);
            ft_topoplotER(cfg_topo, plot3); 
            T=int32((ii+tstep/2)*1000);
            title([num2str(T) 'ms'],'FontSize',20,'FontWeight','bold');
            axis tight
            axis off
            c=c+1;
            if (c>=36 || abs(ii-te)<0.01), c=1; 
                colorbar;
                export_fig(save3);
                close(f); f=f+1;
            else
                while sum(c==BADc)>0
                    if (c>=36 || abs(ii-te)<0.01), c=1; 
                       export_fig(save3); 
                        close(f); f=f+1; break;
                    end
                    c=c+1;
                end
            end
        end
    end
end

% catch err 
%     return;
% end

