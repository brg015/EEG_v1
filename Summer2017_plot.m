%=========================================================================%
% Scratch Space
%=========================================================================%
cod_dir='E:\SEFER\MatlabCode\';
% Specify where your function files are
code_dir{1}=fullfile(cod_dir,'EEG_v1-master\ft_functions');
code_dir{2}=fullfile(cod_dir,'EEG_v1-master\functions');
code_dir{3}=fullfile(cod_dir,'funtion_files-master');
for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end

%=========================================================================%
load(fullfile('E:\SEFER\EEG\Summer2017\','GA_25_ERP.mat'),'ga'); 
% if N=25 then flip 13 for flipping keys
R=ga.Rem.individual(13,:,:);
F=ga.For.individual(13,:,:);
ga.Rem.individual(13,:,:)=F;
ga.For.individual(13,:,:)=R;
clear R F;
load(fullfile('E:\SEFER\EEG\Summer2017\','GA_15_Z.mat'),'ga'); 
load(fullfile('E:\SEFER\EEG\Summer2017\','GA_15_LP.mat'),'ga'); 
load(fullfile('E:\SEFER\EEG\Summer2017\','GA_15_COH.mat'),'ga'); 

L=true(1,25); LABEL=ga.Rem.label;
%=========================================================================%
% Plotter
%=========================================================================%
% template.ced='mw64_withMastoids_fixed.ced';
% template.ica='mw64.ced';
% template.elp='/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp';
% template.plt='E:\SEFER\MatlabCode\EEG_v1-master\templates\layoutmw64_ActiChamp.mat';
% template.plt='E:\SEFER\MatlabCode\EEG_v1-master\templates\layoutmw64.mat';
% template.plt2='E:\SEFER\MatlabCode\EEG_v1-master\templates\layoutmw64_martyPlot.mat';
% template.plt3='E:\SEFER\MatlabCode\EEG_v1-master\templates\layoutmw64.mat';

template.erp='E:\SEFER\MatlabCode\EEG_v1-master\templates\layoutMartyFixed.mat';
template.topo='E:\SEFER\MatlabCode\EEG_v1-master\templates\mw64_withMastoids_fixed.ced';

%=========================================================================%
% ERP
%=========================================================================%
ROI1={'21','15'};
ROI2={'Cz' '49' '50' '38'};
ROI3={'53','51'};
ROI4={'36','55'};

cfg=[]; cfg.scale=[-10 10]; cfg.filt=0; 
cfg.subj=L; cfg.ga=ga; cfg.template=template.erp;
ft_plot_v3({'Rem' 'For'},cfg);

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.2; cfg.t=[-.4 1]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'Rem' 'For'}; cfg.L=L; 
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg);
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg.perspective = {'left'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'right'}; jp_topoplotFT_v3(cfg);

L=true(1,25);
L([2,12,20,21,22])=0; % <30 trials N=5
L([11,3,9,12,22])=0; % Bad Mem <0
L([1,13])=0;         % Pilot/Wrong responses
sum(L)

tstruct.legend={'Remembered' 'Forgotten'}; 
tstruct.window=[-.5 .5]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.25 .45]; tstruct.fig=[1 1 1];
tstruct.filter.on=1; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
tstruct.base=[-.24 0];
O1=erp_ss_plot_v3({'CatchHit' 'RemFor'},get_el(ROI1,LABEL),tstruct);

tstruct=[];
tstruct.legend={'Remembered' 'Forgotten'}; 
tstruct.window=[-.5 .5]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.25 .35]; tstruct.fig=[1 1 1];
tstruct.filter.on=1; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
O1=erp_ss_plot_v3({'Rem' 'For'},get_el(ROI1,LABEL),tstruct);

tstruct=[];
tstruct.legend={'Remembered' 'Forgotten'}; 
tstruct.window=[-.5 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.5 1];
tstruct.filter.on=1; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
tstruct.fig=[1 1 1];
O1=erp_ss_plot_v3({'Rem' 'For'},get_el(ROI2,LABEL),tstruct);

tstruct.legend={'Remembered' 'Forgotten'}; 
tstruct.window=[-.5 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.5 1];
tstruct.filter.on=0; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
tstruct.fig=[1 1 1];
O1=erp_ss_plot_v3({'Rem' 'For'},get_el(ROI4,LABEL),tstruct);
%=========================================================================%
% FRQ
%=========================================================================%
L=true(1,15);
cfg=[]; cfg.Inc=0.2; cfg.t=[-.4 1]; cfg.foi=[10 16]; cfg.clim = [];
cfg.type='pow';  cfg.contrast={'Rem' 'For'}; 
cfg.L=L; cfg.ga=ga;
cfg.template=template.topo;
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet');
cfg.perspective = {'left'}; jp_topoplotFT_v3(cfg); colormap('jet');


cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el({'36' '38'},LABEL); 
cfg.time=[-.4 1.2]; cfg.foi=[10 16]; cfg.test=[.7 1];
cfg.contrast={'Rem' 'For'};
cfg.type='logpower'; tstruct.legend={'Remembered' 'Forgotten'};
cfg.tstruct=tstruct; frq_ss_plot_v2(cfg);


cfg=[]; cfg.template=template.erp;
cfg.topo_lim=[-4 4];
ft_plot_freq('Rem','For',ga,cfg); colormap('jet');













