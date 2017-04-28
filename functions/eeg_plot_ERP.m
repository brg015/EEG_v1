function eeg_plot_ERP()
%=========================================================================%
% BRG 2016
%=========================================================================%
global RUN; 
load(RUN.template.plt3);
%=========================================================================%
% Timlockage...
%=========================================================================%
% for iSubj=1:length(RUN.dir.subjects)
%     iSubj
%     load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
%     cfg.sstep=50; cfg.smax=40; cfg.tmax=100; cfg.win=[-.2 .5];
%     [lrej,~,~,~,~] = ft_reject(data,cfg);
%     data.lat_rej40=lrej;
%     
%      cfg.sstep=50; cfg.smax=30; cfg.tmax=100; cfg.win=[-.2 .5];
%     [lrej,~,~,~,~] = ft_reject(data,cfg);
%     data.lat_rej30=lrej;
%     
%      cfg.sstep=50; cfg.smax=20; cfg.tmax=100; cfg.win=[-.2 .5];
%     [lrej,~,~,~,~] = ft_reject(data,cfg);
%     data.lat_rej20=lrej;
%     
%     save(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']),'data');
% end

% for iSubj=1:length(RUN.dir.subjects)
%     load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
%     LatRej40(iSubj)=sum(data.lat_rej40);
%     LatRej30(iSubj)=sum(data.lat_rej30);
%     LatRej20(iSubj)=sum(data.lat_rej20);
%     SpkRej(iSubj)=sum(data.reject);
%     All30Rej(iSubj)=sum([data.lat_rej30 | data.reject]);
%     All20Rej(iSubj)=sum([data.lat_rej20 | data.reject]);
%     All40Rej(iSubj)=sum([data.lat_rej40 | data.reject]);
%     N(iSubj)=length(data.trial);
% end
% Rmat=[LatRej40;LatRej30;LatRej20;SpkRej;All40Rej;All30Rej;All20Rej;N]';

% for iSubj=1:length(RUN.dir.subjects)
% %     if RUN.dir.plot(iSubj)==0, continue; end
% %     if ~exist(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_timelock.mat']),'file'),
%         % Load in data
%         load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
%       
%         % -------------- Study dependent ---------------------------------%
%         % If subject has events removed, this must be detected and there
%         % behave data needs to be updated as well
%         %
%         % Picks up missed events, but vector of 'real' events is still
%         % shorter - cause epochs are corrupted by bounds
%         % 
%         % Need to figure out what additional trials are being tossed ...
%         % those with
%         tkill=[]; % trials to remove from behave for alignment
%         if length(data.trial)<640
%             display('Missing epochs detected');
%             if iSubj==6 % 3451 by context
%                 tkill=setdiff(1:640,data.accepted_events);
%             elseif iSubj==17 
%                 % 3486 by context (636 accepted out of 638)
%                 % Removed 337 and 350 manually
%                 I=1:640; I([337,350])=[]; % Manually exclude
%                 I2=setdiff(1:max(data.accepted_events),data.accepted_events);
%                 I(I2)=[]; % Excluded via too little overlap
%                 tkill=setdiff(1:640,I); clear I I2;
%             elseif iSubj==30
%                 % 3517 by context (634 accepted out of 638)
%                 I=1:640; I([85,220])=[]; % Manually exclude
%                 I2=setdiff(1:max(data.accepted_events),data.accepted_events);
%                 I(I2)=[]; % Excluded via too little overlap
%                 tkill=setdiff(1:640,I); clear I I2;
%             elseif iSubj==31
%                 % 3519 by context (619 accepted out of 620)
%                 I=1:640; I([181:200])=[]; % Manually exclude
%                 I2=setdiff(1:max(data.accepted_events),data.accepted_events);
%                 I(I2)=[]; % Excluded via too little overlap
%                 tkill=setdiff(1:640,I); clear I I2;
%             end
%         end
%         reject=[data.lat_rej30 | data.reject];
%         [booleans,N,RT,target]=SEA_bool(iSubj,reject,tkill);
%         save(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_booleans.mat']),'booleans','RT','target');
% 
%         % Timelock via booleans - keeping this independent for now
%         fn = fieldnames(booleans);
%         for iConditions = 1 : length(fn)
%             % ERP
%             cfg=[];             
%             cfg.trials = booleans.(fn{iConditions});
%             cfg.keeptrials='no';
%             if sum(cfg.trials)>0
%                 timelock.(fn{iConditions}){1}=ft_timelockanalysis(cfg, data); 
%             else
%                 timelock.(fn{iConditions}){1}=[]; 
%             end
%         end
%         save(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_timelock.mat']),'timelock');
% %     end
% end
%=========================================================================%
%% Grand Averages
%=========================================================================%
% Micro to decide which subjects to include based upon including trial
% counts
% for ii=1:length(RUN.dir.subjects)
%     load(fullfile(RUN.dir.pro,[RUN.dir.subjects{ii} '_booleans.mat']));
%     fn=fieldnames(booleans);
%     for jj=1:length(fn)
%         N(ii,jj)=sum(booleans.(fn{jj}));
%     end
% end
% % 
keyboard;
ga=eeg_ga_erp();
fn=fieldnames(ga);
L={'Rem','RemHi','RemLo','For','ForHi','ForLo','Catch','Targets'};
for ii=1:length(L)
    ga.([L{ii} '_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+8}),ga.(fn{ii}),'erp'); 
end

for ii=1:length(L)
    ga.([L{ii} '_R_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,[],ga.(fn{ii}),'erp'); 
end

for ii=1:length(L)
    ga.([L{ii} '_L_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+8}),[],'erp'); 
end

% for ii=1:length(L)
%     [ga.([L{ii} '_Contra']),ga.([L{ii} '_Ipsi'])]  = contraAndIpsi(layoutmw64,ga.(fn{ii+8}),ga.(fn{ii}),'erp'); 
% end
RUN.dir.sav='F:\Data2\SEA\Enc_EEG\group_ERP\';
% N=sum(RUN.dir.plot);
save(fullfile(RUN.dir.sav,['GA_N' num2str(25) '_erps_CNS_lat30.mat']),'ga');
N=32;
load(fullfile(RUN.dir.sav,['GA_N' num2str(N) '_erps_lat30.mat']),'ga');
%=========================================================================%
%% Define Contrasts
%=========================================================================%
sdisp('Plotting',1);
load(RUN.template.plt);
template{2}=layoutmw64; clear layoutmw64;
template{2}.scale=[1 1];
load(RUN.template.plt2);
template{1}=layoutmw64_martyPlot; clear layoutmw64_martyPlot;
template{1}.scale=[2.3 1.5];

RUN.ga=ga;
RUN.template2=template;
%=========================================================================%
% All plotting functionality starts from this point on
%
% Key functions are...
% fn has all contrast names
% ft_plot({'contrast A','contrast B',subjectID,[min max],'title',interactive)
% -> set subjectID = 0 to plot all subjects
% -> set interactive = 1 for plots to stay displayed
% -> plots BLUE,RED,GREEN, then BLACK wrt contrasts
% -> if interactive = 0 saves plots to RUN.dir.sav
%=========================================================================%
% Subject List
%=========================================================================%
keyboard;
RUN.dir.plot=ones(1,32);
% 32 already removed for Tourettes
RUN.dir.plot(23)=0; % No memory
RUN.dir.plot(25)=0; % No memory
RUN.dir.plot(6)=0;  % low memory, bad data, low catch
RUN.dir.plot(9)=0;  % low memory, bad data
RUN.dir.plot(3)=0;  % artifects
RUN.dir.plot(10)=0; % artifects
RUN.dir.plot(32)=0;

L0=RUN.dir.plot; % Initial optimistic set

L1=RUN.dir.plot;
L1([15,16,18,19,24,27])=0; % Soft rejections
    
% Weird Subjects
L3=zeros(1,31);
L3(10)=1;

%=========================================================================%
% Plot all channels
%=========================================================================%
cfg.scale=[-10 10]; cfg.filt=0; cfg.subj=L3;
ft_plot_v2({'Targets','Catch'},cfg);

cfg.scale=[-5 5]; cfg.filt=0; cfg.subj=Lin;
ft_plot_v2({'Rem_ContraMinusIpsi','For_ContraMinusIpsi'},cfg);
%=========================================================================%
% Plot topos
%=========================================================================%
Lin=logical(ones(1,25));
% Basic Targets
cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-4 4]; 
cfg.perspective = {'top'};  cfg.contrast={'Targets'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-4 4]; 
cfg.perspective = {'left'};  cfg.contrast={'Targets_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-3 3]; 
cfg.perspective = {'back'};  cfg.contrast={'Left_Targets'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'};  cfg.contrast={'Targets_L_ContraMinusIpsi'};
jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-3 3]; 
cfg.perspective = {'back'};  cfg.contrast={'Right_Targets'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'};  cfg.contrast={'Targets_R_ContraMinusIpsi'};
jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-3 3]; 
cfg.perspective = {'back'};  cfg.contrast={'RemHi_R_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'};  cfg.contrast={'RemHi_L_ContraMinusIpsi'};
jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'};  cfg.contrast={'RemHi_ContraMinusIpsi'};
jp_topoplotFT_v2(cfg);


% Memory
cfg=[]; cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-2 2]; 
cfg.perspective = {'top'};  cfg.contrast={'Rem','For'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg=[]; cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-2 2]; 
cfg.perspective = {'top'};  cfg.contrast={'RemHi','For'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-.5 .5]; 
cfg.perspective = {'back'};  cfg.contrast={'Rem_ContraMinusIpsi','For_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'left'}; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-.5 .5]; 
cfg.perspective = {'back'};  cfg.contrast={'RemHi_ContraMinusIpsi','For_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'left'}; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-.5 .5]; 
cfg.perspective = {'back'};  cfg.contrast={'Rem_R_ContraMinusIpsi','For_R_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-.5 .5]; 
cfg.perspective = {'back'};  cfg.contrast={'Rem_L_ContraMinusIpsi','For_L_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);


%=========================================================================%
% Plot channels
%=========================================================================%

O1=erp_ss_plot_v2({'RemHi_ContraMinusIpsi' 'For_ContraMinusIpsi'},...
    get_el(ROI2),[0.300 0.400],tstruct);

O1=erp_ss_plot_v2({'Rem','For'},...
    get_el({'13' '14' '3' '38' '49' '50' 'Cz'}),[0.600 .800],tstruct);






