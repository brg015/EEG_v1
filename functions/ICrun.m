%% initiate glob settings
set_study
addpath(eeglabDir)
eeglab -redraw;
dataPath = 'data/';
addpath('matlab_apps/functions')



%what does this do?
rand('state',2000);
%% 
for iSubj = 1
    EEG = pop_loadbv(fullfile(dataPath,subjectID{iSubj}), [subjectID{iSubj} '.vhdr']);
    EEG = bv2event(EEG);
    EEG = pop_chanedit(EEG, 'lookup', fullfile(fileparts(which('eeglab.m')), ...
            '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'), ...
        'load',{'mw64.ced' 'filetype' 'autodetect'});


%filter data, resample at new rate
    EEG = pop_eegfiltnew(EEG, 0.5, 70); %min - 0.5?
    EEG = pop_resample(EEG,250);

%interpolate any bad channels (set in "set_study.m")

    if(length(interp{iSubj})~=0)
         EEG = pop_interp(EEG,interp{iSubj},'spherical');
    end
%create epochs based on these port codes (targets and cues)
%[-.5 1.2] OR [-1 2]? <- which is better
    EEG = pop_epoch( EEG, {  '101'  '102'  '103'  '104' ...
        '105'  '106'  '107'  '108'  '109'  '110'  '111' ...
        '112'  '113'  '114'  '115'  '116'  '117'  '118'  '119' ...
        '120'  '121'  '122'  '123'  '124'  '125'  '126'  '127'  ...
        '128'  '129'  '130'  '131'  '132'  '41'  '42'  } ...
        ,[-.5 1.5], 'newname', ' resampled epochs', 'epochinfo', 'yes');
    

%remove baseline [in miliseconds]
    EEG = pop_rmbase( EEG, [-200 0]);

%mark trials with extreme values for rejection ([electrodes:],min,max,timestart,timeend,)
%leave eye blinks intact
    EEG = pop_eegthresh(EEG,1,setdiff(1:63,31),-150,1500,-1,2,0,0);
    ICrun.artfThreshIC = EEG.reject.rejthresh;
    EEG = pop_rejepoch(EEG,ICrun.artfThreshIC,0);
    
   
    %Remove bad channels before ICA 
      if(length(interp{iSubj})~=0)
        EEG = pop_select(EEG,'nochannel' ,interp{iSubj});
      end
    
      
    EEG = pop_runica( EEG, 'runica' )
    
    
    ICrun.icaact = EEG.icaact;
    ICrun.icawinv = EEG.icawinv;
    ICrun.icasphere = EEG.icasphere;
    ICrun.icaweights = EEG.icaweights;
    ICrun.icachansind = EEG.icachansind;
    
    %save ICrun for eah subject
    save([dataPath subjectID{iSubj}  '/' ['ICrun' subjectID{iSubj}] '.mat'],'ICrun')
    EEG = pop_saveset( EEG, 'filename',['ICrun' subjectID{iSubj} '.set'],'filepath',[dataPath subjectID{iSubj} '\']);
end

