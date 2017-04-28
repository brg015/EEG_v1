function eeg_plot_SEFER()
%=========================================================================%
% BVD - BRG edits Summer 2016
%=========================================================================%
% Extensive edits needed from the initial version, the original was a huge
% pain to work with easily
global RUN; 
if strcmp(RUN.dir.sess,'enc'), 
    pdir='1'; RUN.dir.save='D:\Data\SEFER\EEG\Winter2016\Enc\';
              RUN.pre=RUN.enc.pre;
else
    pdir='2'; RUN.dir.save='D:\Data\SEFER\EEG\Winter2016\Ret\';
              RUN.pre=RUN.ret.pre;
end

for iSubj=1:length(RUN.dir.subjects)
    if RUN.dir.analyze(iSubj)==0, continue; end
    
    if exist(fullfile(RUN.dir.save,'timelock',[RUN.dir.subjects{iSubj} '_timelock.mat']),'file')
       continue;
    end
    
    post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}];
    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,filesep);
    load(fullfile(save_pro_dir,post_save),'data'); % var=data;

    [booleans,RT] = createConditions_v2(data.trialStruct,iSubj); 
    save(fullfile(RUN.dir.save,'booleans',[RUN.dir.subjects{iSubj} '.mat']),'booleans','RT');
    
    fn = fieldnames(booleans);
    for iConditions=1:length(fn)
        % ERP
        cfg=[];             
        cfg.trials = booleans.(fn{iConditions});
        cfg.keeptrials='no';
        if sum(cfg.trials)>0
            timelock.(fn{iConditions}){1}=ft_timelockanalysis(cfg, data); 
         else
            timelock.(fn{iConditions}){1}=[];
        end
    end
    save(fullfile(RUN.dir.save,'timelock',[RUN.dir.subjects{iSubj} '_timelock.mat']),'timelock');
end

RUN.dir.plot=ones(1,length(RUN.dir.subjects));
RUN.dir.pro=fullfile(RUN.dir.save,'timelock');

sfile=fullfile(RUN.dir.save,'group',['N_' num2str(sum(RUN.dir.plot)) '_ERP.mat']);
if ~exist(sfile,'file')
    ga=eeg_ga_erp();
    save(sfile,'ga');
else
    load(sfile);
end

sdisp('Plotting',1);
load(RUN.template.plt);
template{2}=layoutmw64; clear layoutmw64;
template{2}.scale=[1 1];
load(RUN.template.plt2);
template{1}=layoutmw64_martyPlot; clear layoutmw64_martyPlot;
template{1}.scale=[2.3 1.5];

RUN.ga=ga;
RUN.template2=template;

keyboard;
%=========================================================================%
% Subject List
%=========================================================================%
L0=ones(1,17); % Encoding all
L0=ones(1,16);
%=========================================================================%
% Plot all channels
%=========================================================================%
cfg.scale=[-10 15]; cfg.filt=0; cfg.subj=L0;
ft_plot_v2({'Conc_Hit','Conc_Miss'},cfg);

for ii=1:17
    figure(ii);
    plot(squeeze(ga.Conc_Hit.individual(ii,29,:)));
end
%=========================================================================%
% Plot topos
%=========================================================================%
cfg=[];
cfg.latencyInc=0.2;
cfg.latency = -.2:cfg.latencyInc:0.8;
cfg.latency = [cfg.latency' cfg.latency'+cfg.latencyInc];
cfg.foi = []; cfg.perspective = {'top'}; cfg.clim = [-3 3]; cfg.ncontours=10;

X=load(RUN.template.plt3); cfg.layout = X.layoutmw64.elec;
X=load('F:\Data2\SEA\elec_pnt.mat'); cfg.layout.pnt=X.X;
cfg.lpfilter = 'yes'; cfg.lpfreq = 30; % Gentle filter

% d=ga.For_ContraMinusIpsi;   plot1 = ft_preprocessing(cfg,d); 
% e=ga.Rem_ContraMinusIpsi;   plot2 = ft_preprocessing(cfg,e); 
d=ga.Conc_Hit;   plot1 = ft_preprocessing(cfg,d); 
e=ga.Conc_Miss;   plot2 = ft_preprocessing(cfg,e);

d.individual=plot1.trial-plot2.trial;
d.individual=d.individual(logical(L0),:,:);
jp_topoplotFT(cfg,d)
%=========================================================================%
% Plot channels
%=========================================================================%
L1=zeros(1,17); L1(14)=1; L1(1)=1;

RUN.GAL=RUN.ga.CatchHit.label;
RUN.TML=RUN.template2{1}.label;
RUN.dir.QAL=RUN.dir.save;

tstruct.window=[-.5 2];
tstruct.legend={'Remembered' 'Forgotten'};
tstruct.filter.on=0;
tstruct.filter.lpfreq=30;
tstruct.filter.lpfilter='no';
tstruct.subj=L1;

tstruct.ind=0;
tstruct.write='D:\Data\SEFER\EEG\Winter2016\Enc\group\';
tstruct.folder='N2pc_TC_101216';

O1=erp_ss_plot_v2({'Conc_Hit','Conc_Miss'},...
    get_el({'38'}),[0.700 0.900],tstruct);



