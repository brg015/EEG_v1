function eeg_icrun()
%=========================================================================%
% BVD - BRG edits Summer 2016
%=========================================================================%
global RUN; eeglab -redraw;

pre.sf=RUN.set.sf;
pre.filter_ica=RUN.set.filter_ica;
pre.filter_dta=RUN.set.filter_dta;
pre.thresh_ica=RUN.set.thresh_ica;
pre.thresh_dta=RUN.set.thresh_dta;
pdir='';
for iSubj =1:length(RUN.dir.subjects) % loop through subjects
    if RUN.dir.plot(iSubj)==0, continue; end
    sdisp(RUN.dir.subjects{iSubj},2);
    if (exist(fullfile(RUN.dir.pre,[RUN.dir.subjects{iSubj} '_ica.set']),'file') && RUN.dir.overwrite==0)
        display('...skipping');
        continue;
    end
%         
    % Create subject directory

    save_1=fullfile(RUN.dir.QAL,'ICA');
    if ~exist(save_1,'dir'), mkdir(save_1); end
    
    %---------------------------------------------------------------------%
    % EEG Preprocess 
    %---------------------------------------------------------------------% 
    pre_save=fullfile(RUN.dir.pre,[RUN.dir.subjects{iSubj} '_filter.set']);
    if (~exist(pre_save,'file'))
        EEG=eeg_basic(0,pre,iSubj,pdir);
        EEG=pop_saveset(EEG,'filename',[RUN.dir.subjects{iSubj} '_filter.set'],'filepath',RUN.dir.pre);
    else
        EEG=pop_loadset('filename',[RUN.dir.subjects{iSubj} '_filter.set'],'filepath',RUN.dir.pre);
    end

    %-------------------------%
    % Modify epoched events by study
    %-------------------------%
    switch RUN.dir.study
        case 'SEA'
            EEG = pop_epoch(EEG,{'5'},RUN.pre.epoch, 'newname', ...
                ' resampled epochs', 'epochinfo', 'yes');
    end
    
    %remove baseline [in miliseconds]
    EEG = pop_rmbase( EEG, RUN.pre.baseline);

    %mark trials with extreme values for rejection ([electrodes:],min,max,timestart,timeend,)
    %leave eye blinks intact
    EEG = pop_eegthresh(EEG,1,setdiff(1:63,31),pre.thresh_ica(1),pre.thresh_ica(2),RUN.pre.epoch(1),RUN.pre.epoch(2),0,0);

    ICrun.artfThreshIC = EEG.reject.rejthresh;      % is the list of epochs to remove
    EEG = pop_rejepoch(EEG,ICrun.artfThreshIC,0);   % Removes the epochs from EEG
    %-------------------------%
    %ICA Calculations
    %-------------------------%
    EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
        RUN.template.elp),'load',{RUN.template.ica 'filetype' 'autodetect'});
    
    if ~isempty(RUN.subj{iSubj}.pre.interp)
        EEG = pop_select(EEG,'nochannel',RUN.subj{iSubj}.pre.interp);
    end
    EEG = pop_runica( EEG, 'runica' ,'chanind',[]);

    %save ICrun for eah subject 
   
    pop_topoplot(EEG,0,1:EEG.nbchan);
    set(gcf,'position',[0 0 1280 1024]);
    export_fig(fullfile(save_1,[RUN.dir.subjects{iSubj} '.png']));
    IC_file=[RUN.dir.subjects{iSubj} '_ica.set'];
    EEG = pop_saveset(EEG,'filename',IC_file,'filepath',RUN.dir.pre);

end


