function eeg_plot_ERP_SEA_Ret()
%=========================================================================%
% BRG 2016
%=========================================================================%
global RUN; 
load(RUN.template.plt3);
keyboard;
%=========================================================================%
% Timlockage...
%=========================================================================%
% for iSubj=1:length(RUN.dir.subjects)
%     load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
%     LatRej40(iSubj)=sum(data.lat_rej40);
%     LatRej20(iSubj)=sum(data.lat_rej20);
%     SpkRej(iSubj)=sum(data.reject);
%     All40Rej(iSubj)=sum([data.lat_rej40 | data.reject]);
%     All20Rej(iSubj)=sum([data.lat_rej20 | data.reject]);
%     N(iSubj)=length(data.trial);
% end
% Rmat=[LatRej40;LatRej20;SpkRej;All40Rej;All20Rej]';
% for iSubj=1:length(RUN.dir.subjects)
%     if RUN.dir.plot(iSubj)==0, continue; end
% %     if ~exist(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_timelock.mat']),'file'),
%         % Load in data
%         load(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '.mat']));
%         % -------------- Study dependent ---------------------------------%
%         % If subject has events removed, this must be detected and there
%         % behave data needs to be updated as well
%         %
%         % Picks up missed events, but vector of 'real' events is still
%         % shorter - cause epochs are corrupted by bounds
%         % 
%         % Need to figure out what additional trials are being tossed ...
%         % those with
%         %
%         % This was taken from the encoding setup, but I don't see why the
%         % general setup wouldn't work for retrieval data as well
%         tkill=[]; % trials to remove from behave for alignment
%         if length(data.trial)<768
%             display('Missing epochs detected');
%             % Missing events are detectable via 'boundary' events within
%             % blocks
%             % 3517 by context (634 accepted out of 638)
%             if iSubj==12 % % 3457
%                 I=1:768; I(10)=[]; % Manually exclude
%                 I2=setdiff(1:max(data.accepted_events),data.accepted_events);
%                 I(I2)=[]; % Excluded via too little overlap
%                 tkill=setdiff(1:768,I); clear I I2;
%             elseif iSubj==3 % 3448
%                 load(fullfile(RUN.dir.pro,[RUN.dir.subjects{4} '.mat']));
%                 % Cruddy data likely
%                 % Block 8 1:10, last 3
%                 % Block 9 9:10
%                 % Block 10 1,5,10,11
%                 % Block 11 1
%                 % Block 19 6,23,24,27,28,31,32
% %                 tkill=1:65; % Data is garbage - but preserves indi
% %                 mEX=[[1:10,22,23,24]+8*24,[9:10]+9*24,[1,5,10,11]+10*24,1+24*11,[6,23,24,27,28,31,32]+19*24];
% %                 I=1:768; I(mEX)=[]; % Manually exclude
% %                 I2=setdiff(1:max(data.accepted_events),data.accepted_events);
% %                 I(I2)=[]; % Excluded via too little overlap
% %                 tkill=setdiff(1:768,I); clear I I2 mEX;
%             end
%         end
% 
%         [booleans,N,RT,targets]=SEA_bool_ret(iSubj,(data.lat_rej30 | data.reject),tkill);
%         save(fullfile(RUN.dir.pro,[RUN.dir.subjects{iSubj} '_booleans.mat']),'booleans','RT','targets');
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

%3509 is empty 'cause they have no trials...
ga=eeg_ga_erp();
fn=fieldnames(ga);
L={'Rem','RemHi','RemLo','For','ForHi','ForLo'};
for ii=1:length(L)
    ga.([L{ii} '_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+12}),ga.(fn{ii}),'erp'); 
    ga.([L{ii} '_L_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,ga.(fn{ii+12}),[],'erp'); 
    ga.([L{ii} '_R_ContraMinusIpsi'])  = contraMinusIpsi(layoutmw64,[],ga.(fn{ii}),'erp'); 
end
keyboard;
N=sum(RUN.dir.plot);
RUN.dir.plot(3)=1;
save(fullfile(RUN.dir.sav,['GA_N' num2str(N) '_lat30_CNS_erps.mat']),'ga');
% load(fullfile(RUN.dir.sav,['GA_N' num2str(N) '_erps.mat']),'ga');

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
RUN.dir.sav='F:\Data2\SEA\Ret_EEG\group_ERP\';
%=========================================================================%
% Subject List
%=========================================================================%
L0=ones(1,17); L0(4)=0; L0(7)=0;
L2=zeros(1,17); L2([1,3,7,8,11,14,15])=1;

L=false(1,31);
RUN.dir.plot=true(1,31);
% 32 already removed for Tourettes
RUN.dir.plot(23)=0; % No memory
RUN.dir.plot(25)=0; % No memory
RUN.dir.plot(6)=0;  % low memory, bad data, low catch
RUN.dir.plot(9)=0;  % low memory, bad data
RUN.dir.plot(3)=0;  % artifects
RUN.dir.plot(10)=0; % artifects
L=RUN.dir.plot;
%=========================================================================%
% Plot all channels
%=========================================================================%
Lin=L; 

cfg=[]; cfg.Inc=0.025; cfg.t=[-.1 .2]; cfg.clim = [-.3 .3]; 
cfg.perspective = {'top'};  cfg.contrast={'Rem_ContraMinusIpsi' 'For_ContraMinusIpsi'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);
cfg.perspective = {'left'}; jp_topoplotFT_v2(cfg);



Lin=L; 
cfg.scale=[-.5 .5]; cfg.filt=0; cfg.subj=Lin;
ft_plot_v2({'Rem_ContraMinusIpsi' 'For_ContraMinusIpsi'},cfg);


cfg=[]; cfg.Inc=0.05; cfg.t=[-.1 .5]; cfg.clim = [-4 4]; 
cfg.perspective = {'top'};  cfg.contrast={'Rem'};
cfg.L=Lin; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);

cfg.scale=[-5 5]; cfg.filt=0; cfg.subj=Lin;
ft_plot_v2({'Rem_ContraMinusIpsi','For_ContraMinusIpsi'},cfg);

cfg.scale=[-10 10]; cfg.filt=1; cfg.subj=L0;
ft_plot_v2({'Rem','For','CR','FA'},cfg);

cfg.scale=[-2 2]; cfg.filt=0; cfg.subj=L2;
ft_plot_v2({'Rem_ContraMinusIpsi','For_ContraMinusIpsi'},cfg);

RUN.GAL=RUN.ga.Rem.label;
RUN.TML=RUN.template2{1}.label;
RUN.dir.QAL=RUN.dir.sav;
RUN.dir.sav='F:\Data2\SEA\Ret_EEG\group_ERP\';

tstruct.window=[-.2 1];
tstruct.legend={'Remembered' 'Forgotten'};
tstruct.filter.on=0;
tstruct.filter.lpfreq=10;
tstruct.filter.lpfilter='no';
tstruct.subj=Lin;

ROI1={'55','57','53','25','63','47','45'};
ROI2={'55' '53'};
ROI3={'53'};
ROI4={'23' '25' '27'};
ROI5={'26' '24' '28'};

tstruct.ind=0;
tstruct.write='F:\Data2\SEA\Ret_EEG\QA\v0\';
tstruct.folder='Memory_013117';
O1=erp_ss_plot_v2({'Rem_ContraMinusIpsi' 'For_ContraMinusIpsi'},...
    get_el(ROI4),[0.05 .1],tstruct);



%=========================================================================%
% Plot topos
%=========================================================================%
cfg=[];
cfg.latencyInc=0.1;
cfg.latency = -.05:cfg.latencyInc:0.25;
cfg.latency = [cfg.latency' cfg.latency'+cfg.latencyInc];
cfg.foi = []; cfg.perspective = {'back'}; cfg.clim = [-0.5 .5]; cfg.ncontours=10;
cfg.elecsize=2;

X=load(RUN.template.plt3); cfg.layout = X.layoutmw64.elec;
X=load('F:\Data2\SEA\elec_pnt.mat'); cfg.layout.pnt=X.X;
cfg.lpfilter = 'no'; cfg.lpfreq = 30; % Gentle filter

d=ga.For_ContraMinusIpsi;   %plot1 = ft_preprocessing(cfg,d); 
e=ga.Rem_ContraMinusIpsi;   % plot2 = ft_preprocessing(cfg,e); 
% d=ga.CR;   plot1 = ft_preprocessing(cfg,d); 
% e=ga.Right_RemHi;   plot2 = ft_preprocessing(cfg,e);

d.individual=e.individual-d.individual;
d.individual=d.individual(logical(L2),:,:);
% Kick damn mastoid
d.label(end)=[];
d.individual(:,62,:)=[];
jp_topoplotFT_v2(cfg,d)
%=========================================================================%
% Plot channels
%=========================================================================%
RUN.GAL=RUN.ga.Rem.label;
RUN.TML=RUN.template2{1}.label;
RUN.dir.QAL=RUN.dir.sav;
RUN.dir.sav='F:\Data2\SEA\Ret_EEG\group_ERP\';

tstruct.window=[-.8 0.8];
tstruct.legend={'Rem' 'For'};
tstruct.filter.on=0;
tstruct.filter.lpfreq=10;
tstruct.filter.lpfilter='no';
tstruct.subj=L;
tstruct.base=[0 .1];

tstruct.ind=0;
tstruct.write='F:\Data2\SEA\Ret_EEG\QA\v0\';
tstruct.folder='ContraLatReact_021617_lp10';

ROI1={'55','57','53','25','63','47','45'};
ROI2={'55' '53'};
ROI3={'53'};
ROI4={'36' '45' '37' '46'}; % Occip
ROI5={'3' '13'}; %midfront
ROI6={'23' '25' '27'}; %left temporal
ROI7={'24' '26' '28'}; %right temporal

O1=erp_ss_plot_v2({'Rem_ContraMinusIpsi' 'For_ContraMinusIpsi'},get_el(ROI6),[.1 .2],tstruct);

O1=erp_ss_plot_v2({'Rem_L_ContraMinusIpsi','Rem_R_ContraMinusIpsi',...
    'For_L_ContraMinusIpsi','For_R_ContraMinusIpsi'},get_el(ROI7),[-.4 0],tstruct);






