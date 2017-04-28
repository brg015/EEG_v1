function [ga,ga_freq,ga_bfreq,ga_phase]=eeg_ga(F_method,pdir,pre)
global RUN;
ga={};
ga_phase={};
ga_freq={};
ga_bfreq={};

csd_flag=0;

FreqOn=RUN.plt.FreqOn; % Force freq GA
DatOn=RUN.plt.DatOn;  % Force GA

sdisp('Grand Averaging',1); c=1; 
if (DatOn || FreqOn)
    for iSubj = 1:length(RUN.dir.subjects) 
        tic;
        %-----------------------------------------------------------------%
        % Setup load strings (t=0s
        %-----------------------------------------------------------------%
        [~,save_prefix,~]=fileparts(['subj_' RUN.dir.subjects{iSubj} '_' RUN.save_str{5}]);
        if strcmp(RUN.dir.sess,'ret')
            switch RUN.dir.retset
                case 2, data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix '_conc_timelock.mat']);
                        data_file_freq=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);

                case 3, data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix '_perc_timelock.mat']);
                        data_file_freq=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);

            end
        else
            data_file=fullfile(RUN.dir.pro,RUN.dir.subjects{iSubj},pdir,[save_prefix 'timelock.mat']);
            data_file_freq=fullfile(RUN.adv.save,[RUN.dir.subjects{iSubj} '_freqinfo_' F_method RUN.dir.sess '.mat']);
        end

        %-----------------------------------------------------------------%
        % Load data file ~5s
        if DatOn==1, load(data_file); fn = fieldnames(timelock); end
        X=whos('timelock');
        sdisp(['Read: ' RUN.dir.subjects{iSubj} ' - Size: ' num2str(round(X.bytes./(1024^2))) ' MB'],1); 
        
        if FreqOn==1,
            t1=load([data_file(1:end-12) '.mat']);
        inc_trials=~t1.data.trialStruct.reject.thresh_indx';
            % Load freq file ~1.3s, preproc ~7s
            load(data_file_freq); X=whos('freqlock');
            sdisp(['Freq Size: ' num2str(round(X.bytes./(1024^2))) ' MB'],2);

            if csd_flag==1
                freqlock.dimord='subj_chan_freq_time';
                freqlock.powspctrm=abs(freqlock.fourierspctrm).^2;
                freqlock.phase=angle(freqlock.fourierspctrm);
            end

            cfg_f=[];
            cfg_f.baselinetype='ztransformSession';
            cfg_f.baseline=[-1 2.5]; % DO NOT CHANGE!
            cfg_f.mu=[];
            cfg_f.sigma=[];
            cfg_f.inc_trials=inc_trials;
            % Ztransform of time data
            Fz=ft_freqbaseline(cfg_f,freqlock);           
 
            cfg_f=[];
            cfg_f.baselinetype='db';
            cfg_f.baseline=[-.8 -.4]; % DO NOT CHANGE!
            cfg_f.mu=[];
            cfg_f.sigma=[];
            cfg_f.inc_trials=inc_trials;
            FdB=ft_freqbaseline(cfg_f,freqlock); 
            
            if csd_flag==1
               Fang=FdB;
               Fang.powspctrm=freqlock.phase;
            end
%             figure(1);
%             for ii=1:2
%                 switch ii
%                     case 1, D=Fz; 
%                     case 2, D=FdB; 
%                 end
%                 subplot(2,2,ii);
%                 freq_method_check(D,freqlock); 
%             end
            
            % Epoch the data based on fn, ~280s
            for ii=1:length(fn)
                % Fz basline
                if ~isempty(timelock.(fn{ii}){1})
                    TINC=timelock.(fn{ii}){1}.cfg.trials;
                    freqlock.(fn{ii}){1}=Fz;
                    freqlock.(fn{ii}){1}.powspctrm(~TINC,:,:,:)=[];
                    % Can't save trials here, so simplify
                    freqlock.(fn{ii}){1}=ft_freqdescriptives([],freqlock.(fn{ii}){1});
                    if csd_flag==1
                       afreqlock.(fn{ii}){1}= freqlock.(fn{ii}){1};
                       ITPCraw=Fang.powspctrm(TINC,:,:,:);
                       ITPC=abs(mean(exp(i*ITPCraw)));
                       afreqlock.(fn{ii}){1}.powspctrm=squeeze(ITPC);
                       clear ITPCraw ITPC;
                    end
                else
                    freqlock.(fn{ii}){1}=[];
                end
                % Fzb basline
                if ~isempty(timelock.(fn{ii}){1})
                    TINC=timelock.(fn{ii}){1}.cfg.trials;
                    bfreqlock.(fn{ii}){1}=FdB;
                    bfreqlock.(fn{ii}){1}.powspctrm(~TINC,:,:,:)=[];
                    % Can't save trials here, so simplify
                    bfreqlock.(fn{ii}){1}=ft_freqdescriptives([],bfreqlock.(fn{ii}){1});
                else
                    bfreqlock.(fn{ii}){1}=[];
                end
            end 
            clear Fz FdB
        end
        sdisp(['Processing: ' save_prefix],1);
            
        % Group timelock
        if c==1 
            if FreqOn==1, 
                gfreqlock=freqlock;  
                gbfreqlock=bfreqlock; 
                gafreqlock=afreqlock;
                fn = fieldnames(gbfreqlock); 
            end
            if DatOn==1, gtimelock=timelock;  fn = fieldnames(gtimelock); end   
        else
            for ii=1:length(fn), 
                if DatOn==1, 
                    switch csd_flag
                        case 0
                            gtimelock.(fn{ii}){c}=timelock.(fn{ii}){1}; 
                        case 1
                            cfg=[];
                            cfg.elec=t1.data.elec;
                            tmp1=timelock.(fn{ii}){1};
                            tmp2=ft_scalpcurrentdensity([],tmp1);
                            tmp2.dof=tmp1.dof;
                            gtimelock.(fn{ii}){c}=tmp2;
                    end
                end
                
                if FreqOn==1, 
                    gfreqlock.(fn{ii}){c}=freqlock.(fn{ii}){1}; 
                    gbfreqlock.(fn{ii}){c}=bfreqlock.(fn{ii}){1};   
                    gafreqlock.(fn{ii}){c}=afreqlock.(fn{ii}){1};
                end
            end    
        end
        clear timelock freqlock; c=c+1; toc
    end % Subject Loop

    cfg = []; 
    cfg.keepindividual = 'yes';
    cfg.parameter='avg';
    cfg.baseline = pre.baseline;
    
    cfg_freq=[];
    cfg_freq.keepindividual = 'yes';
    %cfg_freq.parameter = 'powspctrm{1}.powspctrm';

    cfg_gco = [];
    cfg_gco.keepindividual = 'yes';
    %=====================================================================%
    
    for iConditions = 1 : length(fn)
        % Grand average
        GAset=~cellfun('isempty',gtimelock.(fn{iConditions}));
         if sum(GAset)>0
            sdisp(['Averaging: ' fn{iConditions}],1);
            if DatOn==1, 
                ga.(fn{iConditions}) = ft_timelockgrandaverage(cfg, gtimelock.(fn{iConditions}){GAset}); 
            end
            if FreqOn==1, 
                % Take subject means
                ga_freq.(fn{iConditions}) = ft_freqgrandaverage(cfg_freq, gfreqlock.(fn{iConditions}){GAset}); 
                ga_bfreq.(fn{iConditions}) = ft_freqgrandaverage(cfg_freq, gbfreqlock.(fn{iConditions}){GAset});
            end
            clear GAset;
        else
            if DatOn==1, ga.(fn{iConditions}) = []; end
            if FreqOn==1, 
                ga_freq.(fn{iConditions}) = []; 
                ga_bfreq.(fn{iConditions}) = [];
            end
        end
    end
    
    % Add phase angle stuff
    ga_phase=ga_freq;
    if csd_flag==1
        for ii=1:length(fn)
            for jj=1:length(gafreqlock.(fn{ii}))
                ga_phase.(fn{ii}).powspctrm(jj,:,:,:)=gafreqlock.(fn{ii}){jj}.powspctrm;
            end
        end
    end
end