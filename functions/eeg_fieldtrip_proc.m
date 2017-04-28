function eeg_fieldtrip_proc(iSubj)
global RUN; pdir='';
%=========================================================================%
% Setup file names and other presets
%=========================================================================%
keyboard;

[~,save_prefix,~]=fileparts(['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]);
if strcmp(RUN.dir.sess,'ret')
    switch RUN.dir.retset
        case 2, data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix '_conc_timelock.mat']);
        case 3, data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix '_perc_timelock.mat']);
    end
    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},'2');
else
    data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix 'timelock.mat']);
    save_pro_dir=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},'1');
end

% Setup save files and skip if these already exist
DatOn=RUN.plt.DatOn;

if (exist(data_file,'file') && RUN.dir.overwrite==0 && RUN.plt.DatOn==1),      
    DatOn=0;
end

if (DatOn==0)
    return;
end

%=========================================================================%
post_save=['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]; % Based on trials*
[~,tstr,~]=fileparts(RUN.save_str{4});
pre_save_reref=['subj_' RUN.dir.subjects{iSubj} '_' tstr '_reref.set'];
save_dir=fullfile(RUN.dir.pre,RUN.dir.subjects{iSubj},pdir,filesep);

if ~exist(fullfile(save_dir,pre_save_reref),'file'),
    display('Subject not preprocessed - Skipping'); return;
end
root_dir='D:\Data\';
addpath(genpath(fullfile(root_dir,'Geib\Scripts\EEG\eeglab13_2_2b\')));
EEG=pop_loadset('filename',pre_save_reref,'filepath',save_dir);
%=========================================================================%
% Load EEG measures into data for fieldtrip
%=========================================================================%
reject_data.thresh_ch=EEG.reject.rejthreshE;
reject_data.thresh_indx=EEG.reject.rejthresh;
reject_data.trend_ch=EEG.reject.rejconstE;
reject_data.trend_indx=EEG.reject.rejconst;
    
%------------------%
% Minifiles to look at removed events
%------------------%
%list of events
trialStruct = extractTrialStruct(EEG,iSubj);

%create a column of 1's, 0's (reject, don't reject) for each epoch
trialStruct.rejthresh = EEG.reject.rejthresh';
trialStruct.reject=reject_data;
%save in fieldtrip
data = eeglab2fieldtrip(EEG,'preprocessing');
data.trialStruct = trialStruct;
data.EEGepoch = EEG.epoch;
data.EEGurevent= EEG.urevent;
% dt_plot(data,[],trialStruct,3,3);

if ~exist(save_pro_dir,'dir'), mkdir(save_pro_dir); end
save(fullfile(save_pro_dir,post_save),'data');

rmpath(genpath(fullfile(root_dir,'Geib\Scripts\EEG\eeglab13_2_2b\')));
%=========================================================================%
% Process that data
%=========================================================================%    
% Creat those trial conditions
if RUN.dir.hc==0
    % Right now hc=0 is only needed to divide HC and LC trials, more
    % complicated measures are unnecessary here
    booleans_n = createConditions(data.trialStruct,iSubj,[]); % LC 
    fn_n=fieldnames(booleans_n);
     for jj=1:length(fn_n),
        if ~isempty(strfind(fn_n{jj},'_R1')) || ~isempty(strfind(fn_n{jj},'Freq'))
            booleans.([fn_n{jj}])=booleans_n.(fn_n{jj});
        end
    end
    
    booleans_LC = createConditions(data.trialStruct,iSubj,1); % LC 
    fn_LC = fieldnames(booleans_LC);
    for jj=1:length(fn_LC),
        if ~isempty(strfind(fn_LC{jj},'_R1')) || ~isempty(strfind(fn_LC{jj},'Freq'))
            booleans.(['LC_' fn_LC{jj}])=booleans_LC.(fn_LC{jj});
        end
    end
    booleans_HC = createConditions(data.trialStruct,iSubj,2); % HC
    fn_HC = fieldnames(booleans_HC);
    for jj=1:length(fn_HC),
        if ~isempty(strfind(fn_HC{jj},'_R1')) || ~isempty(strfind(fn_HC{jj},'Freq'))
            booleans.(['HC_' fn_HC{jj}])=booleans_HC.(fn_HC{jj});
        end
    end

    fn = fieldnames(booleans);
    for ii=1:4
        switch ii
            case 1, M=booleans.(fn{5});   str_M='HF_';    wind=[1:2,17:18,33:34];
            case 2, M=booleans.(fn{6});  str_M='LF_';    wind=[1:2,17:18,33:34];
            case 3, M=booleans.(fn{7});   str_M='L_HF_';  wind=[1:2,17:18,33:34];
            case 4, M=booleans.(fn{8});   str_M='L_LF_';  wind=[1:2,17:18,33:34];
        end
        for jj=wind % Only R1 components
            booleans.([str_M fn{jj}])=and(booleans.(fn{jj}),M);
        end
    end
    
else
    booleans = createConditions(data.trialStruct,iSubj,[]); 
    
    booleans.stm_p_1_R2crhit=or(booleans.stm_p_1_R2hit,booleans.stm_p_1_R2cr);
    booleans.stm_p_1_R2famiss=or(booleans.stm_p_1_R2miss,booleans.stm_p_1_R2fa);

    fn = fieldnames(booleans);
    %============%
    % Create conditions for HF and LF
    %============%
%     for ii=1:6,
%         switch ii
%             case 1, M=booleans.(fn{12});   str_M='HF_';    wind=[4:11];
%             case 2, M=booleans.(fn{13});   str_M='LF_';    wind=[4:11];
%             case 3, M=booleans.(fn{14});   str_M='L_HF_';  wind=[4:13];
%             case 4, M=booleans.(fn{15});   str_M='L_LF_';  wind=[4:13];
%         end
%         for jj=wind % Only R1 components
%             booleans.([str_M fn{jj}])=and(booleans.(fn{jj}),M);
%         end
%     end
%     clear fn M;
end

fn = fieldnames(booleans);
% report_bool_trial(data,booleans,'R2miss')
% For all those condtions we found
for iConditions = 1 : length(fn)
    % ERP
    cfg=[];             
    cfg.trials = booleans.(fn{iConditions});
    cfg.keeptrials='no';
    if sum(cfg.trials)>0
        if DatOn==1, 
            timelock.(fn{iConditions}){1}=ft_timelockanalysis(cfg, data); 
        end
     else
        if DatOn==1, timelock.(fn{iConditions}){1}=[]; end
    end
end

%=========================================================================%
% Save the output
%=========================================================================%
if DatOn==1, save(data_file,'timelock','-v7.3'); end


    
