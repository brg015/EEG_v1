%% initiate glob settings
set_study
addpath(eeglabDir)
eeglab -redraw;
dataPath = 'data\';
eeglabDir = 'matlab_apps\eeglab13_2_2b\';
addpath('matlab_apps/functions')

rand('state',2000);
%% 
%DO ICrun BEFORE PREPROCESSING

for iSubj = 1:length(subjectID)
    %load header file
    EEG = pop_loadbv(fullfile(dataPath,subjectID{iSubj}), [subjectID{iSubj} '.vhdr']);
    %remove "S" before each port code
    EEG = bv2event(EEG);
    EEG = pop_chanedit(EEG, 'lookup', fullfile(fileparts(which('eeglab.m')), ...
            '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'), ...
        'load',{'mw64.ced' 'filetype' 'autodetect'});


%filter data, resample at new rate
    EEG = pop_eegfiltnew(EEG, 0, 70, 10); %min - 0.5?
    EEG = pop_resample(EEG,250);

%remove bad channels
    tmpChn = EEG.chanlocs
     if(length(interp{iSubj})~=0)
        EEG = pop_select(EEG,'nochannel' ,interp{iSubj});
     end
     
     
     %loads dataset after ICA run is complete
    load([dataPath subjectID{iSubj} '\' 'ICrun' subjectID{iSubj} '.mat'])
    
    EEG.icawinv = ICrun.icawinv
    EEG.icasphere = ICrun.icasphere
    EEG.icaweights = ICrun.icaweights
    EEG.icachansind = ICrun.icachansind
    
    %subtract eye blink component
    EEG = pop_subcomp(EEG,ICA{iSubj})
    
    
        
%interpolate channels (set in "set_study.m")
    
    EEG = pop_interp(EEG,tmpChn,'spherical');


%create epochs based on these events (targets and cues)
%[-.5 1.2] OR [-1 2]? <- which is better
    EEG = pop_epoch( EEG, {  '101'  '102'  '103'  '104' ...
        '105'  '106'  '107'  '108'  '109'  '110'  '111' ...
        '112'  '113'  '114'  '115'  '116'  '117'  '118'  '119' ...
        '120'  '121'  '122'  '123'  '124'  '125'  '126'  '127'  ...
        '128'  '129'  '130'  '131'  '132'  '41'  '42'  } ...
        ,[-.5 1.5], 'newname', ' resampled epochs', 'epochinfo', 'yes');
    

%remove baseline [in miliseconds]
    EEG = pop_rmbase( EEG, [-200 0]);
    
    
%list of events
    trialStruct = extractTrialStruct(EEG);
    
%reject trials with extreme values ([electrodes:],min,max,timestart,timeend,)   
    EEG = pop_eegthresh(EEG,1,[1:63] ,-100,100,-0.5,1.2,0,0);

    %create a column of 1's, 0's (reject, don't reject) for each epoch
    trialStruct.rejthresh = EEG.reject.rejthresh';
    trialStruct.artfThreshIC = ICrun.artfThreshIC';

%save in fieldtrip
    data = eeglab2fieldtrip(EEG,'preprocessing');
    data.trialStruct = trialStruct;
    data.EEGepoch = EEG.epoch;
    data.EEGurevent= EEG.urevent;

    save([dataPath subjectID{iSubj}  '/' subjectID{iSubj} '.mat'],'data')
end

