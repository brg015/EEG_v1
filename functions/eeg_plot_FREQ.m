function eeg_plot_FREQ()
%=========================================================================%
% BRG 2016
%=========================================================================%
global RUN; 
load(RUN.template.plt3);
%=========================================================================%
%% Freq Data
%=========================================================================%
% cfg_freq.method='mtmconvol';
% cfg_freq.output='pow';
% cfg_freq.taper='hanning';
% cfg_freq.channel='all';
% cfg_freq.keeptrials='yes';
% cfg_freq.keeptapers='no';
% cfg_freq.foi=3:1:30;
% cfg_freq.t_ftimwin=3./cfg_freq.foi; % 5 cycles per window
% cfg_freq.toi=-1:0.05:1.8;
% suffix='_ENS_HAN_v2';

% keyboard;
% 
% cfg_freq            = [];
% cfg_freq.keeptrials = 'yes';
% cfg_freq.output     = 'pow';
% cfg_freq.method     = 'mtmconvol';
% cfg_freq.taper = 'dpss';
% cfg_freq.foi =  logspace(log10(4),log10(80),50);
% cfg_freq.toi        = -0.2:0.05:1.5;
% cfg_freq.tapsmofrq  = 5*log10(cfg_freq.foi);
% cfg_freq.t_ftimwin =   [3./cfg_freq.foi(cfg_freq.foi<=3) ...
%     4./cfg_freq.foi(cfg_freq.foi> 3 & cfg_freq.foi<=7) ...
%     5./cfg_freq.foi(cfg_freq.foi> 7 & cfg_freq.foi<=14) ...
%     7./cfg_freq.foi(cfg_freq.foi> 14 &  cfg_freq.foi<=20) ...
%     10./cfg_freq.foi(cfg_freq.foi> 20 &  cfg_freq.foi<=50) ...
%     15./cfg_freq.foi(cfg_freq.foi> 50)];
% cfg_freq.keeptapers = 'no';
% cfg_freq.pad        = 'maxperlen';
% suffix='_ENS_DPSS';

% cfg_freq            = [];
% cfg_freq.keeptrials = 'yes';
% cfg_freq.output     = 'pow';
% cfg_freq.method     = 'mtmconvol';
% cfg_freq.taper = 'dpss';
% cfg_freq.foi =  logspace(log10(4),log10(40),50);
% cfg_freq.toi        = -1:0.05:1.7;
% cfg_freq.tapsmofrq  = 5*log10(cfg_freq.foi);
% cfg_freq.t_ftimwin =   [3./cfg_freq.foi(cfg_freq.foi<=3) ...
%     4./cfg_freq.foi(cfg_freq.foi> 3 & cfg_freq.foi<=7) ...
%     5./cfg_freq.foi(cfg_freq.foi> 7 & cfg_freq.foi<=14) ...
%     7./cfg_freq.foi(cfg_freq.foi> 14 &  cfg_freq.foi<=20) ...
%     10./cfg_freq.foi(cfg_freq.foi> 20 &  cfg_freq.foi<=50)];
% cfg_freq.keeptapers = 'no';
% cfg_freq.pad        = 'maxperlen';
% suffix='_ENS_DPSS_v2';

cfg_freq            = [];
cfg_freq.keeptrials = 'yes';
cfg_freq.output     = 'pow';
cfg_freq.method     = 'mtmconvol';
cfg_freq.taper = 'dpss';
cfg_freq.foi =  logspace(log10(4),log10(30),34);
cfg_freq.toi        = -1:0.05:1.7;
cfg_freq.tapsmofrq  = 5*log10(cfg_freq.foi);
cfg_freq.t_ftimwin =   [3./cfg_freq.foi(cfg_freq.foi<=3) ...
    4./cfg_freq.foi(cfg_freq.foi> 3 & cfg_freq.foi<=7) ...
    5./cfg_freq.foi(cfg_freq.foi> 7 & cfg_freq.foi<=14) ...
    7./cfg_freq.foi(cfg_freq.foi> 14 &  cfg_freq.foi<=20) ...
    10./cfg_freq.foi(cfg_freq.foi> 20 &  cfg_freq.foi<=50)];
cfg_freq.keeptapers = 'no';
cfg_freq.pad        = 'maxperlen';
suffix='_ENS_DPSS_v3';

for iSubj=2:length(RUN.dir.subjects)
    if RUN.dir.plot(iSubj)==0, continue; end
    if ~exist(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_dB' suffix '.mat']),'file')
        % Load in data
        load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
        load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_booleans.mat']));
        fn=fieldnames(booleans);
        
        booleans.Targets=(booleans.Rem | booleans.For);

        data=ft_remove_ems(data,(data.reject | data.lat_rej30));
%         csd_data=ft_scalpcurrentdensity([],data);
        freqlock=ft_freqanalysis(cfg_freq, data);
        
        cfg_f.baselinetype='ztransformSession';
        cfg_f.inc_trials=~logical(data.reject | data.lat_rej30);
        cfg_f.baseline=[-.5 .2]; % Not used to z-transform
        cfg_f.ztrials=(booleans.Rem | booleans.For);
        Fz=ft_freqbaseline(cfg_f,freqlock);
        
        cfg_f.baselinetype='logpower';
        cfg_f.inc_trials=~logical(data.reject | data.lat_rej30);
        cfg_f.baseline=[-.5 .2]; % Not used to z-transform
        dB=ft_freqbaseline(cfg_f,freqlock);
        frej=logical(data.reject | data.lat_rej30);
        for ii=1:2
            switch ii
                case 1, fdata=Fz; save_file=fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_Fz' suffix '.mat']);
                        save_file2=fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_Fz_all' suffix '.mat']);
                case 2, fdata=dB; save_file=fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_dB' suffix '.mat']);
                        save_file2=fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_dB_all' suffix '.mat']);
            end
            save(save_file2,'fdata','frej');
            
            for jj=1:length(fn)
                sfdata.(fn{jj}){1}=fdata;
                if sum(~booleans.(fn{jj}))~=length(booleans.(fn{jj}))
                    sfdata.(fn{jj}){1}.powspctrm(~booleans.(fn{jj}),:,:,:)=[];
                    % Can't save trials here, so simplify
                    sfdata.(fn{jj}){1}=ft_freqdescriptives([],sfdata.(fn{jj}){1});
                else
                    sfdata.(fn{jj}){1}=[]; 
                end
            end
            save(save_file,'sfdata'); clear sfdata fdata;
        end
    end
end

N=sum(RUN.dir.plot);
if ~exist(fullfile(RUN.dir.sav,['GA_N' num2str(N) '_freq_ENS.mat']),'file')
    [ga_freq_dB,ga_freq_Fz]=eeg_ga_freq(suffix);
    fn=fieldnames(ga_freq_dB);
    L={'Rem','RemHi','RemLo','For','ForHi','ForLo','Catch','Targets'};
    for ii=1:8
        ga_freq_dB.([L{ii} '_CI'])  = contraMinusIpsi(layoutmw64,ga_freq_dB.(fn{ii+8}),ga_freq_dB.(fn{ii}),'power'); 
        ga_freq_Fz.([L{ii} '_CI'])  = contraMinusIpsi(layoutmw64,ga_freq_Fz.(fn{ii+8}),ga_freq_Fz.(fn{ii}),'power'); 
         ga_freq_dB.([L{ii} '_L_CI'])  = contraMinusIpsi(layoutmw64,ga_freq_dB.(fn{ii+8}),[],'power'); 
        ga_freq_Fz.([L{ii} '_L_CI'])  = contraMinusIpsi(layoutmw64,ga_freq_Fz.(fn{ii+8}),[],'power'); 
         ga_freq_dB.([L{ii} '_R_CI'])  = contraMinusIpsi(layoutmw64,[],ga_freq_dB.(fn{ii}),'power'); 
        ga_freq_Fz.([L{ii} '_R_CI'])  = contraMinusIpsi(layoutmw64,[],ga_freq_Fz.(fn{ii}),'power'); 
    end
    save(fullfile(RUN.dir.sav,['GA_N' num2str(N) '_freq' suffix '.mat']),'ga_freq_dB','ga_freq_Fz');
end

