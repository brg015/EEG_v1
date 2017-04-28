function ga=eeg_ga_erp()
global RUN;
ga={};
sdisp('Grand Averaging',1);  
f=0; c=0;
for iSubj = 1:length(RUN.dir.subjects) 
    data_file=fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_timelock.mat']);
    %-----------------------------------------------------------------%
    if RUN.dir.plot(iSubj)==1
        load(data_file); fn = fieldnames(timelock); c=c+1;
        if f==0, f=1; end
    else
        continue;
    end
    
    % Group timelock
    if f==1
        gtimelock=timelock; f=2;
    else
        for ii=1:length(fn) 
            gtimelock.(fn{ii}){c}=timelock.(fn{ii}){1}; 
        end    
    end
    clear timelock 
end % Subject Loop

cfg = []; 
cfg.keepindividual = 'yes';
cfg.parameter='avg';
cfg.baseline = RUN.pre.baseline;
%=====================================================================%

for iConditions = 1 : length(fn)
    % Grand average
    GAset=~cellfun('isempty',gtimelock.(fn{iConditions}));
    if sum(GAset)>0
        sdisp(['Averaging: ' fn{iConditions}],1);
        ga.(fn{iConditions}) = ft_timelockgrandaverage(cfg, gtimelock.(fn{iConditions}){GAset}); 
        clear GAset;
    else
        continue;
    end
    % If a subject has no data, then flush it with NaNs
    I=cellfun('isempty',gtimelock.(fn{iConditions}));
    if sum(I)>0
        dat_temp=ga.(fn{iConditions}).individual;
        % This is gonna break if 1st subject is empty 8|
        for ii=1:length(I)
            if I(ii)==1,
                dat_repl(ii,:,:)=dat_temp(ii,:,:);
            else
                dat_repl(ii,:,:)=NaN(size(dat_temp,2),size(dat_temp,3));
            end
        end
        ga.(fn{iConditions}).individual=dat_repl;
        clear dat_repl dat_temp;
    end
    clear I;
end
    
