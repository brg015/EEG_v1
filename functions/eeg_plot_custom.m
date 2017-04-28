function eeg_plot_custom()
%=========================================================================%
% BVD - BRG edits Summer 2016
%=========================================================================%
% Extensive edits needed from the initial version, the original was a huge
% pain to work with easily
global RUN; load layoutmw64.mat;

% fn = fieldnames(booleans);
% % report_bool_trial(data,booleans,'R2miss')
% % For all those condtions we found
% for iConditions = 1 : length(fn)
%     % ERP
%     cfg=[];             
%     cfg.trials = booleans.(fn{iConditions});
%     cfg.keeptrials='no';
%     if sum(cfg.trials)>0
%         if DatOn==1, 
%             timelock.(fn{iConditions}){1}=ft_timelockanalysis(cfg, data); 
%         end
%      else
%         if DatOn==1, timelock.(fn{iConditions}){1}=[]; end
%     end
% end

%=========================================================================%
%% Freq Data
%=========================================================================%
F_method='v3';
% cfg_freq.method='mtmconvol';
% cfg_freq.output='fourier';
% cfg_freq.taper='hanning';
% cfg_freq.channel='all';
% cfg_freq.keeptrials='yes';
% cfg_freq.keeptapers='yes';
% cfg_freq.foi=3:0.5:30;
% cfg_freq.t_ftimwin=7./cfg_freq.foi; % 7 cycles per window
% cfg_freq.toi=-1:0.05:2.5;
% 
% % All frequency information has been moved to 'Connectivity'
% for iSubj=1:length(RUN.dir.subjects), 
% %=========================================================================%
% % 1) Load in subject data
% %=========================================================================%
%     if strcmp(RUN.dir.sess,'ret')
%         switch RUN.dir.retset
%             case 2, SV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
%             case 3, SV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
%                     continue;
%         end
%     else
%         SV=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
%     end
%     if (exist(SV,'file') && RUN.dir.overwrite==0), continue; end
%     continue;
%     sdisp(['Frequing Subject: ' RUN.dir.subjects{iSubj}],2);
%     % These subjects have 2 retrieval files which are later concatenated
%     save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,filesep);
% 
%     post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]; % Based on trials*
%     try
%         load(fullfile(save_pro_dir,post_save),'data'); % var=data;
%     catch err
%         display('No Data');
%     end
%     
%     csd_data=ft_scalpcurrentdensity([],data);
%     freqlock=ft_freqanalysis(cfg_freq, csd_data); 
%     freqlock.crsspctrm=[];
%     save(SV,'freqlock'); clear freqlock;
%     continue;
% %=========================================================================%
% % 2) Manipulate frequency information
% %=========================================================================%    
%     % Too expensive to save the crsspectral information, so let's make an
%     % => Fast baseline to get z-scored
%     cfg_f.baselinetype='ztransformSession';
%     cfg_f.baseline=[-1 2.5]; % Not used to z-transform
%     cfg_f.mu=[];
%     cfg_f.sigma=[];
%     Fz=ft_freqbaseline(cfg_f,freqlock);
% %=========================================================================%    
%     % => Get trial types
%     % Lets focus on semantics for the moment R1 mem
%     booleans = createConditions(data.trialStruct,iSubj,[]); 
%     R1hit=booleans.stm_p_1_R1hit;
%     R1miss=booleans.stm_p_1_R1miss;
%     Fz.hit=R1hit;
%     Fz.miss=R1miss;
%     Fdata{iSubj}=Fz;
%     clear freqlock Fz;
% %     % Coherence craziness
% %     cfg_coh.method='coh';
% %     cfg_coh.trials=R1hit;
% %     Xhit = ft_connectivityanalysis(cfg_coh, freqlock);
%     Shit = Fz.powspctrm(R1hit,:,:,:);
% %     fancy_eeg_adj(Xhit,Shit,[],[400 800],[],[RUN.dir.subjects{iSubj} '_R1hit']);
% %     
% %     cfg_coh.trials=R1miss;
% %     Xmiss = ft_connectivityanalysis(cfg_coh, freqlock);
% %     Smiss = Fz.powspctrm(R1miss,:,:,:);
% %     fancy_eeg_adj(Xmiss,Smiss,[],[400 800],[],[RUN.dir.subjects{iSubj} '_R1miss']);
% % 
% %     freqlock.crsspctrm=[]; % For now just clear it 
% %     save(SV,'freqlock');
% %     clear freqlock;
% end
% 
% % save(fullfile(RUN.adv.save,'Fdata.mat'),'Fdata');
% 
% 

%=========================================================================%
%% Grand Averages
%=========================================================================%
if ~exist(fullfile(RUN.dir.pro,'group'),'dir'),
    mkdir(fullfile(RUN.dir.pro,'group'));
end
% Save group averages based upon N - impercise but it'll do for now
[~,u_save_str,~]=fileparts(RUN.save_str{5});
N=length(RUN.dir.subjects);
if strcmp(RUN.dir.sess,'ret')
    switch RUN.dir.retset
        case 2
            gfreq_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_freq_' F_method '.mat']);
            gerps_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_erps.mat']);
        case 3
            gfreq_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_perc_freq_' F_method '.mat']);
            gerps_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_perc_erps.mat']);
    end
else
    gfreq_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_freq_' F_method '.mat']);
    gerps_save=fullfile(RUN.dir.pro,'group',[u_save_str '_N' num2str(N) '_erps.mat']);
end

keyboard;
if RUN.plt.GA==1
    [ga,ga_freq,ga_bfreq,ga_phase]=eeg_ga(F_method,pdir,pre);
    if RUN.plt.FreqOn==1, 
        save(gfreq_save,'ga_freq','ga_bfreq','ga_phase','-v6'); 
    end
    if RUN.plt.DatOn==1,  save(gerps_save,'ga','-v6'); end   
else
    if RUN.plt.FreqOn==1, load(gfreq_save); end
    if RUN.plt.DatOn==1, load(gerps_save); end
end     
%=========================================================================%
%% QA micro
%=========================================================================%\
% eeg_beh_summary()
%=========================================================================%
%% Define Contrasts
%=========================================================================%
sdisp('Plotting',1);
load(RUN.template.plt);
template{2}=layoutmw64; clear layoutmw64;
template{2}.scale=[1 1];
load(RUN.template.plt2);
template{1}=layoutmw64_martyPlot; clear layoutmw64_martyPlot;
template{1}.scale=[2.3 1.5];

RUN.ga=ga;
RUN.template2=template;
RUN.plt.name=u_save_str;

fs=u_save_str;
switch RUN.dir.sess
    case 'enc'
        pre.baseline=[-.400,-.200];
        phase=[1];
    case 'ret'
        pre.baseline=[-.200,0];
        phase=[2 3];
end

RUN.pre=pre;

%=========================================================================%
% Save_Dir=RUN.dir.QAL;
% Get a quick summary from GA
% data{1}.header='subject';
% for ii=1:length(RUN.dir.subjects),
%     data{1}.col{ii}=RUN.dir.subjects{ii};
% end
% 
% fn = fieldnames(ga);
% for ii=1:length(fn)
%     data{ii+1}.header=fn{ii};
%     for jj=1:length(RUN.dir.subjects)
%         try
%             N=sum(ga.(fn{ii}).cfg.previous{jj}.trials);
%         catch err
%             N=0;
%         end
%         data{ii+1}.col(jj)=N; clear N;
%     end
% end
% write_struct(data,fullfile(Save_Dir,'Ga_Trial_Counts.csv')); clear data;
%=========================================================================%
% User plotting area
%=========================================================================%
% for ii=1:18,
%     display(['Subject: ' RUN.dir.subjects{ii}]);
%     A=length(ga.stm_p_2_stm_all.cfg.previous{ii}.trials);
%     B=sum(ga.stm_p_2_stm_all.cfg.previous{ii}.trials);
%     display(['  Detected: ' num2str(A)]);
%     display(['  Kept: ' num2str(B)]);
% end
keyboard;
%=========================================================================%
% R translation area...
%=========================================================================%
RUN.GAL=RUN.ga.stm_p_1_R1hit.label;
RUN.TML=RUN.template2{1}.label;

for iSubj=1:length(RUN.dir.subjects)
    % Raw data
    post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]; % Based on trials*
    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},'1',filesep);
    X2=load(fullfile(save_pro_dir,post_save),'data'); % var=data;
    
    inc_trials=~X2.data.trialStruct.reject.thresh_indx';
           
    A=excel_reader(fullfile('D:\Data\SEFER\behav\eeg\pre',...
        ['subject' RUN.dir.subjects{iSubj} '.csv']));

    X2.data.trialStruct.descrip{26}='RT1';
    X2.data.trialStruct.descrip{27}='RT2';
    X2.data.trialStruct.descrip{28}='invRT1';
    X2.data.trialStruct.descrip{29}='invRT2';
    
    ID=cell2num(A{1}.col); 
    RT1=cell2num(A{23}.col); 
    RT2=cell2num(A{29}.col); 
    Hc=cell2num(A{18}.col);
    R1=cell2num(A{19}.col);
    Ms=cell2num(A{20}.col);
    
    for ii=1:length(X2.data.trialStruct.rProfile(:,19))
        I=find(X2.data.trialStruct.rProfile(ii,19)==ID);
        X2.data.trialStruct.rProfile(ii,26)=RT1(I);
        X2.data.trialStruct.rProfile(ii,27)=RT2(I);
        X2.data.trialStruct.rProfile(ii,28)=1/RT1(I);
        X2.data.trialStruct.rProfile(ii,29)=1/RT2(I);
        X2.data.trialStruct.rProfile(ii,3)=Hc(I);
        X2.data.trialStruct.rProfile(ii,4)=R1(I);
        X2.data.trialStruct.rProfile(ii,5)=Ms(I);
    end
    
%     csd_data=ft_scalpcurrentdensity([],X2.data);
    freqlock=ft_freqanalysis(cfg_freq, X2.data); 
    
    freqlock.dimord='subj_chan_freq_time';
    freqlock.powspctrm=abs(freqlock.fourierspctrm).^2;
    freqlock.phase=angle(freqlock.fourierspctrm);
                
    cfg_f.baselinetype='ztransformSession';
    cfg_f.inc_trials=inc_trials;
    cfg_f.baseline=[-1 2.5]; % Not used to z-transform
    cfg_f.mu=[];
    cfg_f.sigma=[];
    Fz=ft_freqbaseline(cfg_f,freqlock);
    
    cfg_f=[];
    cfg_f.inc_trials=inc_trials;
    cfg_f.baselinetype='db';
    cfg_f.baseline=[-.8 -.4]; % DO NOT CHANGE!
    cfg_f.mu=[];
    cfg_f.sigma=[];
    FdB=ft_freqbaseline(cfg_f,freqlock); 
    
    % Make a data frame
    ft_create_dataframe(X2.data,Fz,FdB,RUN.dir.subjects{iSubj});   
end
keyboard;
RUN.ga_freq=ga_freq;   % Zscoring method
RUN.ga_bfreq=ga_bfreq; % dB baseline
RUN.ga_phase=ga_phase; % Phase information 
% RUN.ga -> ERP information

cfg.data=RUN.ga_freq;
cfg.dim={'ID' 'Channel' 'Freq' 'Time'};
cfg.dim_label={{},{RUN.ga_freq.stm_p_1_R1hit.label},{RUN.ga_freq.stm_p_1_R1hit.freq},{RUN.ga_freq.stm_p_1_R1hit.time}};
cfg.stm_set=fieldnames(RUN.ga_freq);


RUN.template.plt2='D:\Data\SEFER\EEG\template\layoutmw64_martyPlot.mat';
X=load(RUN.template.plt2);
% L=X.layoutmw64_martyPlot.label(X.layoutmw64_martyPlot.pos(:,1)<0);
% R=X.layoutmw64_martyPlot.label(X.layoutmw64_martyPlot.pos(:,1)>0);
% B=X.layoutmw64_martyPlot.label(X.layoutmw64_martyPlot.pos(:,2)<0);
% F=X.layoutmw64_martyPlot.label(X.layoutmw64_martyPlot.pos(:,2)>0);
keyboard

tstruct.window=[-.5 1];
tstruct.legend={'No Response' 'Response' 'Difference'};
O1=erp_ss_plot({'stm_p_1_stm','stm_p_1_Resp'},...
    get_el({'21'}),[.25 0.4],[],'temp',tstruct);

tstruct.window=[-.5 1];
tstruct.legend={'Congruent' 'Incongruent' 'Difference'};
O1=erp_ss_plot({'stm_p_1_Con','stm_p_1_Inc'},...
    get_el({'36'}),[.4 0.5],[],'temp',tstruct);

tstruct.window=[-.5 1];
tstruct.legend={'Remembered' 'Forgotten'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'21' '15' 'LVEOG'}),[.25 .40],[],'temp',tstruct);

RUN.GAL=RUN.ga.stm_p_1_R1hit.label;
RUN.TML=RUN.template2{1}.label;

ft_plot({'stm_p_1_stm','stm_p_1_Resp'},0,[-10 10],'CatchNotT',1);
% Language
ft_plot({'stm_p_1_LL_WL','stm_p_1_LL_WH','stm_p_1_LH_WL','stm_p_1_LH_WH'},0,[-5 5],'Letter',1);
ft_plot({'stm_p_1_LL_WL','stm_p_1_LL_WH'},0,[-5 5],'LetterL_W_LvH_T',0);
ft_plot({'stm_p_1_LH_WL','stm_p_1_LH_WH'},0,[-5 5],'LetterH_W_LvH_T',0);
ft_plot({'stm_p_1_LL_WL','stm_p_1_LH_WL'},0,[-5 5],'WordL_L_LvH_T',0);
ft_plot({'stm_p_1_LL_WH','stm_p_1_LH_WH'},0,[-5 5],'WordH_L_LvH_T',0);
ft_plot({'stm_p_1_Con','stm_p_1_Inc'},0,[-5 5],'ConInc_NotT',1);


ft_plot({'stm_p_1_stm_all','stm_p_1_stm_all'},0,[-18 18],'All',0);

ft_plot({'stm_p_1_stm','stm_p_1_Resp'},0,[-18 18],'CatchNotT',0);
ft_plot({'stm_p_1_stm','stm_p_1_Resp'},0,[-18 18],'CatchNotT',0);

ft_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},0,[-18 18],'R1g',0);
ft_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},0,[-18 18],'R2g',0);

tstruct.window=[-1 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'21'}),[.7 .8],[],'temp',tstruct);
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};

tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'42' '35' '41'}),[0.6 0.8],[],'temp',tstruct);

ft_plot_freq('stm_p_1_R2hit','stm_p_1_R2miss',0,'stm_g',1)


%=========================================================================%
% Load Regression Templates
%=========================================================================%
RUN.ga_freq.Resp=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.Lemma=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.Letter=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.RT=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.Mem=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.Comp=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.WxL=RUN.ga_freq.stm_p_1_stm_all;
RUN.ga_freq.WxM=RUN.ga_freq.stm_p_1_stm_all;

for ii=1:length(RUN.dir.subjects)
    t1=load(fullfile('D:\Data\SEFER\EEG\RegressionResults\',[RUN.dir.subjects{ii} '.mat']));
    RUN.ga_freq.Resp.powspctrm(ii,:,:,:)=t1.RegSave(1,:,:,:);
    RUN.ga_freq.Lemma.powspctrm(ii,:,:,:)=t1.RegSave(2,:,:,:);
    RUN.ga_freq.Letter.powspctrm(ii,:,:,:)=t1.RegSave(3,:,:,:);
    RUN.ga_freq.RT.powspctrm(ii,:,:,:)=t1.RegSave(4,:,:,:);
    RUN.ga_freq.Mem.powspctrm(ii,:,:,:)=t1.RegSave(5,:,:,:);
    RUN.ga_freq.Comp.powspctrm(ii,:,:,:)=t1.RegSave(6,:,:,:);
    RUN.ga_freq.WxL.powspctrm(ii,:,:,:)=t1.RegSave(7,:,:,:);
    RUN.ga_freq.WxM.powspctrm(ii,:,:,:)=t1.RegSave(8,:,:,:);
end
%=========================================================================%
% Gamma!
%=========================================================================%
ft_plot_freq('stm_p_1_R1hit','stm_p_1_R1miss',0,'stm_g',1)
ft_plot_freq('stm_p_1_stm_all','stm_p_1_stm_all',0,'stm_g',1)


f1=10:30; f2='ga_freq';
RUN.ga.lg_R1hit=RUN.ga.stm_p_1_R1hit; RUN.ga.lg_R1hit.time=RUN.(f2).stm_p_1_R1hit.time;
RUN.ga.lg_R1hit.individual=squeeze(mean(RUN.(f2).stm_p_1_R1hit.powspctrm(:,:,f1,:),3));

RUN.ga.lg_R2hit=RUN.ga.stm_p_1_R1hit; RUN.ga.lg_R2hit.time=RUN.(f2).stm_p_1_R1hit.time;
RUN.ga.lg_R2hit.individual=squeeze(mean(RUN.(f2).stm_p_1_R2hit.powspctrm(:,:,f1,:),3));

RUN.ga.lg_R1miss=RUN.ga.stm_p_1_R1hit; RUN.ga.lg_R1miss.time=RUN.(f2).stm_p_1_R1hit.time;
RUN.ga.lg_R1miss.individual=squeeze(mean(RUN.(f2).stm_p_1_R1miss.powspctrm(:,:,f1,:),3));

RUN.ga.lg_R2miss=RUN.ga.stm_p_1_R1hit; RUN.ga.lg_R2miss.time=RUN.(f2).stm_p_1_R1hit.time;
RUN.ga.lg_R2miss.individual=squeeze(mean(RUN.(f2).stm_p_1_R2miss.powspctrm(:,:,f1,:),3));

ft_plot({'lg_R1hit','lg_R1miss'},0,[-.1 .1],'R1g',1);
ft_plot({'lg_R2hit','lg_R2miss'},0,[-.1 .1],'R2g',0);



tstruct.window=[-1 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'15' '21' 'LVEOG'}),[0.1 .2],[40 60],'temp',tstruct);
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};

tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'42' '35' '41'}),[0.4 0.5],[40 60],'temp',tstruct);

tstruct.window=[-1 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'15' '21' 'LVEOG'}),[.1 0.2],[40 60],'temp',tstruct);

tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'lg_R2hit','lg_R2miss'},...
    get_el({'42' '35' '41'}),[-.2 -.1],[],'temp',tstruct);

%=========================================================================%
% Low Frequency Results
%=========================================================================%
ft_plot({'lg_R1hit','lg_R1miss'},0,[-.1 .1],'R1g',0);
ft_plot({'lg_R2hit','lg_R2miss'},0,[-.1 .1],'R2g',0);

ft_plot_freq('stm_p_1_R2crhit','stm_p_1_R2famiss',0,'R2_lf_cr',1)


ft_plot_freq('stm_p_1_stm_all','stm_p_1_stm_all',0,'R1_lf',1)

ft_plot_freq('stm_p_1_R1hit','stm_p_1_R1miss',0,'R1_lf',1)
ft_plot_freq('stm_p_1_R2hit','stm_p_1_R2miss',0,'R2_lf',0)

ft_plot_freq('stm_p_1_stm_all','stm_p_1_stm_all',0,'stm',0)

tstruct.window=[-1 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'36' '43' '44'}),[.5 1],[3 6],'temp',tstruct);


tstruct.window=[-1 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'36' '43' '44'}),[-.7 -.5],[],'temp',tstruct);

tstruct.window=[-1 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'36' '43' '44'}),[.8 .9],[12 20],'temp',tstruct);



R{1}.name='PeakTheta';
R{1}.ROI={'35' '36' '37' '44' '46' '43' '45'};
R{2}.nama='PeakAlpha';
R{2}.ROI={'35' '36' '42' '44' '46' '54' '41' '43' '45' '53'};
R{3}.nama='PeakLAlpha';
R{3}.ROI={'41' '43' '45' '53'};
R{4}.nama='PeakRAlpha';
R{4}.ROI={'42' '44' '46' '54'};

% left frontal 14 7 58 9
% frontal (1 5 6) 1 41 52
% parietal (45 36 46) 27 36 37
%=========================================================================%
% Important Results
%=========================================================================%
ft_plot({'stm_p_1_R2hit','stm_p_1_R2cr','stm_p_1_R2miss','stm_p_1_R2fa'},0,[-15 15],'R2cr',1);

%=========================================================================%
% Important Results ERPs
%=========================================================================%
% LPFC
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'15' '21' 'LVEOG'}),[0.250 0.400],[],'temp',tstruct);
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'15' '21' 'LVEOG'}),[1.5 2],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'15' '21' 'LVEOG'}),[0.250 0.400],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2cr','stm_p_1_R2fa'},...
    get_el({'15' '21' 'LVEOG'}),[0.250 0.400],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2crhit','stm_p_1_R2famiss'},...
    get_el({'15' '21' 'LVEOG'}),[0.250 0.400],[],'temp',tstruct);


% Central Effect (DM)
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[0.700 0.900],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[-0.025 0.025],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[0.700 0.900],[],'temp',tstruct);
tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2cr','stm_p_1_R2fa'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[0.700 0.900],[],'temp',tstruct);
tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2crhit','stm_p_1_R2famiss'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[0.700 0.900],[],'temp',tstruct);

tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'49' '47' '4' '38' '37' '50' '48'}),[-0.025 0.025],[],'temp',tstruct);


% Persistent frontal negativity
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'1' '6' '16' '5' '15'}),[1.4 1.8],[],'temp',tstruct);

% P100 effect
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'42'}),[0.075 0.125],[],'temp',tstruct);
tstruct.window=[-0.5 2];

tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'42'}),[0.075 0.125],[],'temp',tstruct);
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2cr','stm_p_1_R2fa'},...
    get_el({'42'}),[0.075 0.125],[],'temp',tstruct);
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2crhit','stm_p_1_R2famiss'},...
    get_el({'42'}),[0.075 0.125],[],'temp',tstruct);

% 29 L. Temporal T7
tstruct.window=[-0.5 2];
tstruct.legend={'Conceptual Hit' 'Conceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R1hit','stm_p_1_R1miss'},...
    get_el({'29'}),[1 1.5],[],'temp',tstruct);
tstruct.window=[-0.5 2];
tstruct.legend={'Perceptual Hit' 'Perceptual Miss' 'Difference'};
O1=erp_ss_plot({'stm_p_1_R2hit','stm_p_1_R2miss'},...
    get_el({'29'}),[1 1.5],[],'temp',tstruct);





