function eeg_coherence_calc

global RUN;

% Berry Code for coherence


dataPath = '~/storage/projects/visual_search';
projectFile
%% define channel combinations
channelcmbRight ={}
for i = 1:length(layoutmw64.elec.label)
    if length(channelcmbRight) == 0
   channelcmbRight = {'56' layoutmw64.elec.label{i}}
    else
       channelcmbRight = vertcat(channelcmbRight,{'56' layoutmw64.elec.label{i}})
    end
    
end

channelcmbLeft ={}
for i = 1:length(layoutmw64.elec.label)
    if length(channelcmbLeft) == 0
   channelcmbLeft = {'55' layoutmw64.elec.label{i}}
    else
       channelcmbLeft = vertcat(channelcmbLeft,{'55' layoutmw64.elec.label{i}})
    end
end

for j=5:17;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(['participant: ' int2str(j)])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pp = participants{j}
    load(fullfile(dataPath,pp, 'fieldtripRawFeb14.mat'));
    load(fullfile(dataPath,pp, 'trialStruct.mat'));
    booleans = createBooleans(trialStruct);
   
    
   for iTrials = 1:length(data.trial)
 
        data.trial{iTrials}(ismember(data.label,channels.occipitalLeft{1}),:) = mean(data.trial{iTrials}(ismember(data.label,channels.occipitalLeft),:));
        data.trial{iTrials}(ismember(data.label,channels.occipitalRight{1}),:) = mean(data.trial{iTrials}(ismember(data.label,channels.occipitalRight),:));
    end    
    

    fn = fieldnames(booleans);
    for iConditions = 1 : length(fn)
        cfg              = [];
        cfg.channelcmb = channelcmbRight;
        cfg.channels = unique(cfg.channelcmb);    
        cfg.trials=  booleans.(fn{iConditions});
        cfg.keeptrials = 'yes';
        cfg.output       = 'powandcsd';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 4:2:30;  
        cfg.t_ftimwin    = 5./cfg.foi;  
        cfg.toi          = -0.5:0.1:1;               
        connectivityDataRight.(fn{iConditions}){j} = ft_freqanalysis(cfg,data);      
        
        cfg = [];
        cfg.method = 'coh';
        connectivityDataRight.(fn{iConditions}){j} = ft_connectivityanalysis(cfg,connectivityDataRight.(fn{iConditions}){j});
        
    end
    
    for iConditions = 1 : length(fn)
        cfg              = [];
        cfg.channelcmb = channelcmbLeft;
        cfg.channels = unique(cfg.channelcmb);    
        cfg.trials=  booleans.(fn{iConditions});
        cfg.keeptrials = 'yes';
        cfg.output       = 'powandcsd';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 4:2:30;  
        cfg.t_ftimwin    = 5./cfg.foi;  
        cfg.toi          = -0.5:0.1:1;               
        connectivityDataLeft.(fn{iConditions}){j} = ft_freqanalysis(cfg,data);      

        cfg = [];
        cfg.method = 'coh';
        connectivityDataLeft.(fn{iConditions}){j} = ft_connectivityanalysis(cfg,connectivityDataLeft.(fn{iConditions}){j});
    end
    
    for iConditions = 1 : length(fn)
        connectivityDataLeft.(fn{iConditions}){j}.powspctrm =connectivityDataLeft.(fn{iConditions}){j}.cohspctrm;
        connectivityDataLeft.(fn{iConditions}){j}.dimord = 'chan_freq_time';
        connectivityDataLeft.(fn{iConditions}){j}.label = connectivityDataLeft.(fn{iConditions}){j}.labelcmb(:,2);
        connectivityDataLeft.(fn{iConditions}){j} =  rmfield(connectivityDataLeft.(fn{iConditions}){j},'labelcmb');
        
        connectivityDataRight.(fn{iConditions}){j}.powspctrm =connectivityDataRight.(fn{iConditions}){j}.cohspctrm;
        connectivityDataRight.(fn{iConditions}){j}.dimord = 'chan_freq_time';
        connectivityDataRight.(fn{iConditions}){j}.label = connectivityDataRight.(fn{iConditions}){j}.labelcmb(:,2);
        connectivityDataRight.(fn{iConditions}){j} =  rmfield(connectivityDataRight.(fn{iConditions}){j},'labelcmb');
        
    end
end

save('connectivityDataLeft20.mat','connectivityDataLeft')
save('connectivityDataRight20.mat','connectivityDataRight')

%%
load('connectivityDataLeft20.mat')
load('connectivityDataRight20.mat')


fn = fieldnames(connectivityDataLeft);
for iConditions = 1 : length(fn)
    cfg = [];
    cfg.keepindividual = 'yes';
    gaConLeft.(fn{iConditions}) = ft_freqgrandaverage(cfg, connectivityDataLeft.(fn{iConditions}){:});
    gaConLeft.(fn{iConditions}).layout = layoutmw64;
    
    %%replace the seed channels with ones
    idx  = find(strcmp(layoutmw64.label,'55'));
   temp =  gaConLeft.(fn{iConditions}).powspctrm;
    sizeTmp = size(temp);
    sizeTmp(2) = 1;
    temp = gaConLeft.(fn{iConditions}).powspctrm ;
    gaConLeft.(fn{iConditions}).powspctrm =     cat(2, temp(:,1:(idx-1),:,:),ones(sizeTmp),temp(:,(idx):size(temp,2),:,:));
    gaConLeft.(fn{iConditions}).label = layoutmw64.label(1:64);
    
    
    
    
    gaConRight.(fn{iConditions}) = ft_freqgrandaverage(cfg, connectivityDataRight.(fn{iConditions}){:});
    gaConRight.(fn{iConditions}).layout = layoutmw64;
    idx  = find(strcmp(layoutmw64.label,'56'));
   temp =  gaConRight.(fn{iConditions}).powspctrm;
    sizeTmp = size(temp);
    sizeTmp(2) = 1;
    temp = gaConRight.(fn{iConditions}).powspctrm ;
    gaConRight.(fn{iConditions}).powspctrm =     cat(2, temp(:,1:(idx-1),:,:),ones(sizeTmp),temp(:,(idx):size(temp,2),:,:));
    gaConRight.(fn{iConditions}).label = layoutmw64.label(1:64);
    
    
end