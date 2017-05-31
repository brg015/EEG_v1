%=========================================================================%
% Scratch Space
%=========================================================================%
cod_dir='C:\Users\brg13\Desktop\MatlabCode\';
% Specify where your function files are
code_dir{1}=fullfile(cod_dir,'EEG_v1-master\ft_functions');
code_dir{2}=fullfile(cod_dir,'EEG_v1-master\functions');
code_dir{3}=fullfile(cod_dir,'funtion_files-master');
for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end

template.erp=fullfile(cod_dir,'\EEG_v1-master\templates\layoutMartyFixed.mat');
template.topo=fullfile(cod_dir,'\EEG_v1-master\templates\mw64_withMastoids_fixed.ced');
%=========================================================================%
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_ERP.mat'),'ga'); % Has 28
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_Z.mat'),'ga'); 
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_LP.mat'),'ga'); 

L=true(1,28); LABEL=ga.Rem.label;
%=========================================================================%
ROI1={'55','57','53','25','63','47','45'};
ROI2={'55' '53'};
ROI3={'53'};
ROI4={'36' '45' '37' '46'}; % Occip
ROI5={'1' '3' '5' '6' '9' '10'}; %midfront
ROI6={'38','49','47','57'}; % Lpar
ROI7={'Cz' '3' '14' '50' '38' '49' '13'}; % Cz cluster
ROI8={'23' '25' '27'};
ROI9={'37' '38' '48' '46' '36' '43' '45'}; % Pz cluster
%=========================================================================%
% ERP
%=========================================================================%
% cfg=[]; cfg.scale=[-15 15]; cfg.filt=0; 
% cfg.subj=L; cfg.ga=ga; cfg.template=template.erp;
% ft_plot_v3({'Rem' 'For'},cfg);

% DM TOPO
cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.2; cfg.t=[0 .8]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'Rem' 'For'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
% N2pc TOPO
cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 .4]; cfg.clim = [-.5 .5]; cfg.ga=ga;
cfg.contrast={'Rem_CI' 'For_CI'}; cfg.L=L; 
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

% DM Plot
tstruct=[]; tstruct.legend={'Remembered' 'Forgotten' 'Difference'}; 
tstruct.window=[-.2 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.7 .9]; tstruct.fig=[0 1 0 1]; tstruct.delta=1;
tstruct.filter.on=1; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
O1=erp_ss_plot_v3({'Rem' 'For'},get_el(ROI7,LABEL),tstruct);

tstruct=[];
tstruct.legend={'Remembered' 'Forgotten'}; 
tstruct.window=[-.5 .5]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.25 .35]; tstruct.fig=[0 1 0 1]; tstruct.delta=1;
tstruct.filter.on=1; tstruct.filter.lpfilter='yes'; tstruct.filter.lpfreq=20;
O1=erp_ss_plot_v3({'Rem_CI' 'For_CI'},get_el(ROI2,LABEL),tstruct);

%=========================================================================%
% FRQ
%=========================================================================%
addpath(fullfile(cod_dir,'EEG_toolboxes\fieldtrip-20170430\')); ft_defaults;
addpath(fullfile(cod_dir,'EEG_toolboxes\eeglab14_1_0b\'));

% Top Topo REM-FOR
cfg=[]; cfg.Inc=0.25; cfg.t=[-.5 .75]; cfg.foi=[4 7]; cfg.clim = [-.025 .025];
cfg.type='pow';  cfg.contrast={'Rem' 'For'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet');
% Top Topo Targets
cfg=[]; cfg.Inc=0.25; cfg.t=[-.5 .75]; cfg.foi=[4 7]; cfg.clim = [0 .8];
cfg.type='pow';  cfg.contrast={'Targets'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet');
% Back Topo REM-FOR
L=true(28,1);
cfg=[]; cfg.Inc=0.2; cfg.t=[-.2 1]; cfg.foi=[8 14]; cfg.clim = [-.1 .1];
cfg.type='pow';  cfg.contrast={'Rem_CI' 'For_CI'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet');
% Back Topo REM-FOR
L=true(28,1);
cfg=[]; cfg.Inc=0.2; cfg.t=[-.2 1]; cfg.foi=[8 14]; cfg.clim = [-.1 .1];
cfg.type='pow';  cfg.contrast={'Targets_CI'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet');

% Midfrontal Theta
cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el(ROI5,LABEL); 
cfg.time=[-.5 1]; cfg.foi=[3.5 7]; cfg.test=[0 .4];
cfg.contrast={'Rem' 'For'};
cfg.type='logpower'; tstruct.legend={'Remember' 'Forgotten'};
cfg.tstruct=tstruct; O1=frq_ss_plot_v2(cfg);

% Alpha CI
cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el(ROI2,LABEL); 
cfg.time=[-.5 1]; cfg.foi=[8 14]; cfg.test=[.5 1];
cfg.contrast={'Targets'};
cfg.type='logpower'; tstruct.legend={'Remember' 'Forgotten'};
cfg.tstruct=tstruct; O1=frq_ss_plot_v2(cfg);


cfg=[]; cfg.template=template.erp;
cfg.topo_lim=[-4 4];
ft_plot_freq('Rem','For',ga,cfg); colormap('jet');