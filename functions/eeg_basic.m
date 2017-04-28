function EEG=eeg_basic(ica_time,pre,iSubj,pdir)

global RUN;

chopSui=fullfile(RUN.dir.raw,RUN.dir.subjects{iSubj},pdir,[RUN.dir.subjects{iSubj} '_chopped.set']);
if ~exist(chopSui,'file')
    display('You must examine data and ensure goodness');
    keyboard;
else
    EEG=pop_loadset('filename',[RUN.dir.subjects{iSubj} '_chopped.set'],'filepath',fullfile(RUN.dir.raw,RUN.dir.subjects{iSubj},pdir));
end
% Strip S from codes
EEG = bv2event(EEG);

% ica_time allows for data and ICA to be filtered at slightly different
% frequencies
filtorder=33000; % pulled from pop_eegfiltnew suggested order
switch ica_time
    case 0, 
        EEG = pop_eegfiltnew(EEG,[],pre.filter_ica(2)); 
        EEG=pop_firws(EEG,'ftype','highpass','minphase',1,'wtype','hamming',...
            'fcutoff',pre.filter_ica(1),'forder',filtorder);
    case 1, 
        EEG = pop_eegfiltnew(EEG,[],pre.filter_dta(2)); 
        EEG=pop_firws(EEG,'ftype','highpass','minphase',1,'wtype','hamming',...
            'fcutoff',pre.filter_dta(1),'forder',filtorder);
end
    
EEG = pop_resample(EEG,pre.sf);
