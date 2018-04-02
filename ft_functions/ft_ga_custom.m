function ga=ft_ga_custom(cfg)
ga={};
sdisp('Grand Averaging',1);  
%-------------------------------------------------------------------------%
% ERP
%-------------------------------------------------------------------------%
if strcmp(cfg.type,'ERP')
    
for iSubj = 1:length(cfg.file) 
    data_file=cfg.file{iSubj}; load(data_file); fn = fieldnames(timelock); 
    % Baseline [-100 0]
    cfg2.baseline=[-.1 0];
%     for ii=1:length(fn), if ~isempty(timelock.(fn{ii}){1})
%         timelock.(fn{ii}){1}=ft_timelockbaseline(cfg2,timelock.(fn{ii}){1});
%     end; end
%     plot(timelock.Rem{1}.time,timelock.Rem{1}.avg(10,:))

    if iSubj==1
        gtimelock=timelock;
    else
        for ii=1:length(fn), gtimelock.(fn{ii}){iSubj}=timelock.(fn{ii}){1}; end    
    end
    clear timelock;
end

cfg_ga = []; cfg_ga.keepindividual = 'yes'; cfg_ga.parameter='avg';

for iConditions = 1 : length(fn)
    % Grand average
    GAset=~cellfun('isempty',gtimelock.(fn{iConditions}));
    if sum(GAset)>0
        sdisp(['Averaging: ' fn{iConditions}],1);
        ga.(fn{iConditions}) = ft_timelockgrandaverage(cfg_ga, gtimelock.(fn{iConditions}){GAset}); 
        clear GAset;
    else
        continue;
    end
    % If a subject has no data, then flush it with NaNs
    I=cellfun('isempty',gtimelock.(fn{iConditions}));
    if sum(I)>0
        dat_temp=ga.(fn{iConditions}).individual;
        % This is gonna break if 1st subject is empty 8|
        c2=1; for ii=1:length(I)
            if I(ii)==1
                display('...Clear data');
                dat_repl(ii,:,:)=NaN(size(dat_temp,2),size(dat_temp,3));
            else
                dat_repl(ii,:,:)=dat_temp(c2,:,:); c2=c2+1;
            end
        end
        ga.(fn{iConditions}).individual=dat_repl;
        clear dat_repl dat_temp;
    end
    if size(ga.(fn{iConditions}).individual,1)<length(cfg.file)
        keyboard;
    end
    clear I;
end

end
%-------------------------------------------------------------------------%
% FREQ
%-------------------------------------------------------------------------%
if strcmp(cfg.type,'FRQ')
    
for iSubj = 1:length(cfg.file) 
    load(cfg.file{iSubj});
    if iSubj==1 
        gfreqlock=sfdata; fn = fieldnames(gfreqlock); 
    else
        for jj=1:length(fn), gfreqlock.(fn{jj}){iSubj}=sfdata.(fn{jj}){1}; end
    end  
end

cfg_freq=[]; cfg_freq.keepindividual = 'yes';
for iConditions = 1 : length(fn)
    % Grand average
    GAset=~cellfun('isempty',gfreqlock.(fn{iConditions}));
    if sum(GAset)>0
        ga.(fn{iConditions}) = ft_freqgrandaverage(cfg_freq, gfreqlock.(fn{iConditions}){GAset}); 
    else
        % Changed on 8/2/2017
        continue; % ga.(fn{iConditions}) = []; 
    end

    % If a subject has no data, then flush it with NaNs
    I=cellfun('isempty',gfreqlock.(fn{iConditions}));
    if sum(I)>0
        dat_temp=ga.(fn{iConditions}).powspctrm;
        % This is gonna break if 1st subject is empty 8|
        c2=1; for ii=1:length(I)
            if I(ii)==1
                display('...Clear data');
                dat_repl(ii,:,:,:)=NaN(size(dat_temp,2),size(dat_temp,3),size(dat_temp,4));
            else
                dat_repl(ii,:,:,:)=dat_temp(c2,:,:,:); c2=c2+1;
            end
        end
        ga.(fn{iConditions}).powspctrm=dat_repl;
        clear dat_repl dat_temp;   
    end
end   

end
%-------------------------------------------------------------------------%
% Coherence
%-------------------------------------------------------------------------%
if strcmp(cfg.type,'COH')

for iSubj=1:length(cfg.file)
    load(cfg.file{iSubj});
    fn=fieldnames(ITC.itpc);
    for jj=1:length(fn)
        ga.(fn{jj}).powspctrm(iSubj,:,:,:)=ITC.itpc.(fn{jj});
        ga.(fn{jj}).label=ITC.label;
        ga.(fn{jj}).freq=ITC.freq;
        ga.(fn{jj}).time=ITC.time;
        ga.(fn{jj}).dimord='subj_chan_freq_time';
    end
end

end
    