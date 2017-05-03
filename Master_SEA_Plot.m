%=========================================================================%
%% Plotting Setup
%=========================================================================%
clc; clear all; close all;
global RUN;
addpath('C:\Users\brg13\Desktop\SEA\');
%=========================================================================%
%% Encoding & Retrieval (CNS)
%=========================================================================%
wrk_dir=fullfile('C:\Users\brg13\Desktop\SEA\');
% 1a) Encoding ERPs
E_ERP=load(fullfile(wrk_dir,'CNS_2017','DataSets','GA_N25_erps_CNS_lat30.mat'));
LEerp=true(1,25);
% 1b) Encoding FRQs
% suffix='_ENS_DPSS_v2';
% suffix='_ENS_HAN_v2';
% E_FRQ=load(fullfile(wrk_dir,'CNS_2017','DataSets','GA_N25_freq_ENS_DPSS_v2.mat'));
E_FRQ=load(fullfile(wrk_dir,'CNS_2017','DataSets','GA_N25_freq_ENS.mat'));
LEfrq=true(1,25);
E_FRQ.ga_freq_Fz=[];
G=E_FRQ.ga_freq_dB.Catch;
clear E_FRQ;

S=load(fullfile('F:\SEA\Enc_EEG\pro\v0\3446_dB_all_ENS_HAN_v2.mat'));
E_COH=load(fullfile(wrk_dir,'CNS_2017','DataSets','GA_N25_ITC_ENS_COH_v4.mat'));
fn=fieldnames(E_COH.GA);
for ii=1:length(fn)
    E_COH2.GA.(fn{ii}).powspctrm=E_COH.GA.(fn{ii}); E_COH.GA.(fn{ii})=[];
    E_COH2.GA.(fn{ii}).label=S.fdata.label;
    E_COH2.GA.(fn{ii}).freq=S.fdata.freq;
    E_COH2.GA.(fn{ii}).time=S.fdata.time;
    E_COH2.GA.(fn{ii}).dimord=G.dimord;
    E_COH2.GA.(fn{ii}).cfg=G.cfg; 
end
clear S G;
%2a) Retrieval ERPs
% R_ERP=load(fullfile(wrk_dir,'CNS_2017','DataSets','GAret_N25_lat30_CNS_erps.mat'));
% LRerp=true(1,25);
% 2b) Retrieval FRQs
% R_FRQ=load(fullfile(wrk_dir,'CNS_2017','DataSets','GAret_N25_freq_ENS_DPSS.mat'));
% LRfrq=true(1,25);


L=true(1,25);
%=========================================================================%
% Plotter
%=========================================================================%
flush_run(); % load needed RUN variables
%=========================================================================%
ROI1={'55','57','53','25','63','47','45'};
ROI2={'55' '53'};
ROI3={'53'};
ROI4={'36' '45' '37' '46'}; % Occip
ROI5={'3'}; %midfront
ROI6={'38','49','47','57'}; % Lpar
ROI7={'Cz' '3' '14' '50' '38' '49' '13'}; % Cz cluster
ROI8={'23' '25' '27'};
ROI9={'37' '38' '48' '46' '36' '43' '45'}; % Pz cluster

%=========================================================================%
% ERP
%=========================================================================%
cfg.scale=[-.5 .5]; cfg.filt=0; cfg.subj=L; cfg.ga=E_ERP.ga;
ft_plot_v2({'Rem_ContraMinusIpsi' 'For_ContraMinusIpsi'},cfg);

cfg=[]; cfg.Inc=0.2; cfg.t=[0 1]; cfg.clim = [-2 2]; cfg.ga=E_ERP.ga;
cfg.perspective = {'top'};  cfg.contrast={'Catch' 'Targets'};
cfg.L=L; jp_topoplotFT_v2(cfg);
cfg=[]; cfg.Inc=0.2; cfg.t=[0 1]; cfg.clim = [-2 2]; cfg.ga=E_ERP.ga;
cfg.perspective = {'top'};  cfg.contrast={'Rem' 'For'};
cfg.L=L; jp_topoplotFT_v2(cfg);

cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);

tstruct.window=[-.2 1]; tstruct.legend={'Remembered (c-i)' 'Forgotten (c-i)'}; 
tstruct.subj=L; tstruct.ga=E_ERP.ga;
tstruct.base=[-.2 0];
tstruct.filter.on=1;
tstruct.filter.lpfilter='yes';
tstruct.filter.lpfreq=20;
O1=erp_ss_plot_v2({'Catch' 'Targets'},...
    get_el(ROI7),[0 .1],tstruct);
%=========================================================================%
% FRQ
%=========================================================================%
cfg=[]; cfg.Inc=0.1; cfg.t=[-.4 .5]; cfg.clim = [0 .075]; cfg.foi=[3 7];
cfg.type='pow';  cfg.contrast={'Rem'}; 
cfg.L=L; 
% cfg.ga=E_FRQ.ga_freq_Fz;
cfg.ga=E_COH2.GA;
cfg.perspective = {'left'}; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);
cfg.perspective = {'top'}; jp_topoplotFT_v2(cfg);

cfg=[]; cfg.Inc=0.1; cfg.t=[-.5 1]; cfg.foi=[3 7]; cfg.clim = [];
cfg.type='pow';  cfg.contrast={'Targets'}; 
cfg.L=L; cfg.ga=E_FRQ.ga_freq_dB;
cfg.perspective = {'left'}; jp_topoplotFT_v2(cfg);
cfg.perspective = {'back'}; jp_topoplotFT_v2(cfg);
cfg.perspective = {'top'}; jp_topoplotFT_v2(cfg);

cfg=[];
cfg.elec=get_el({'3' '2'}); cfg.time=[-1 1]; cfg.foi=[3 7]; cfg.test=[0 .5];
cfg.contrast={'Rem' 'For'};
cfg.L=L; 
% cfg.ga=E_FRQ.ga_freq_Fz;
cfg.ga=E_COH2.GA;
tstruct.legend={'Remembered' 'Forgotten'};
cfg.tstruct=tstruct; frq_ss_plot_v2(cfg);

cfg=[];
cfg.elec=get_el({'53' '55'}); cfg.time=[-1 1]; cfg.foi=[8 14]; cfg.test=[0 .5];
cfg.contrast={'Rem_CI' 'For_CI'};
cfg.L=L; cfg.ga=E_FRQ.ga_freq_Fz;
tstruct.legend={'Remembered' 'Forgotten'};
cfg.tstruct=tstruct; frq_ss_plot_v2(cfg);



ft_plot_freq('Rem','For',R_FRQ.ga_freq_Fz);














