function [ga_freq_dB,ga_freq_Fz]=eeg_ga_freq(suffix)
global RUN;

for ii=1:2
    f=0; c=0;
    for iSubj = 1:length(RUN.dir.subjects)  
        
        if RUN.dir.plot(iSubj)==0
            continue;
        else
            if f==0, f=1; end; c=c+1;
        end
        
        switch ii
            case 1, load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_Fz' suffix '.mat']));
            case 2, load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_dB' suffix '.mat']));
        end
    
        
        
        % Group timelock
        if f==1 
            gfreqlock=sfdata; 
            fn = fieldnames(gfreqlock); 
            f=2;
        else
            for jj=1:length(fn),
                gfreqlock.(fn{jj}){c}=sfdata.(fn{jj}){1}; 
            end
        end  
    end  

    cfg_freq=[];
    cfg_freq.keepindividual = 'yes';
    for iConditions = 1 : length(fn)
        % Grand average
        GAset=~cellfun('isempty',gfreqlock.(fn{iConditions}));
        if sum(GAset)>0
            ga_freq.(fn{iConditions}) = ft_freqgrandaverage(cfg_freq, gfreqlock.(fn{iConditions}){GAset}); 
        else
            ga_freq.(fn{iConditions}) = []; 
        end
    end   
        
    switch ii
        case 1, ga_freq_Fz=ga_freq;
        case 2, ga_freq_dB=ga_freq;
    end
    clear ga_freq;   
end
    
    