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

load(template.erp);
ga.('RemHi_CI')  = contraMinusIpsi(layoutmw64,ga.('Left_RemHi'),ga.('Right_RemHi'),'erp'); 
ga.('ForHi_CI')  = contraMinusIpsi(layoutmw64,ga.('Left_ForHi'),ga.('Right_ForHi'),'erp'); 
ga.('Fast_CI')  = contraMinusIpsi(layoutmw64,ga.('Left_Fast'),ga.('Right_Fast'),'erp'); 
ga.('Slow_CI')  = contraMinusIpsi(layoutmw64,ga.('Left_Slow'),ga.('Right_Slow'),'erp'); 
%=========================================================================%
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_ERP_beta.mat'),'ga'); % Has 28
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_Z_RT.mat'),'ga'); 
load(fullfile('C:\Users\brg13\Desktop\SEA\Writing\','GA_28_LP_RT.mat'),'ga'); 

L=true(1,28); LABEL=ga.Rem.label;
%=========================================================================%
% ERP
%=========================================================================%
cfg=[]; cfg.scale=[-15 15]; cfg.filt=0; 
cfg.subj=L; cfg.ga=ga; cfg.template=template.erp;
ft_plot_v3({'RemHi' 'ForHi'},cfg);
%-------------------------------------------------------------------------%
% Basic Contrast
%-------------------------------------------------------------------------%
cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-5 5]; cfg.ga=ga;
cfg.contrast={'mean'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-5 5]; cfg.ga=ga;
cfg.contrast={'Living' 'NonLiving'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-3 3]; cfg.ga=ga;
cfg.contrast={'Fast' 'Slow'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-3 3]; cfg.ga=ga;
cfg.contrast={'Left' 'Right'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-5 5]; cfg.ga=ga;
cfg.contrast={'Catch' 'Target'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')
%-------------------------------------------------------------------------%
% Memory
%-------------------------------------------------------------------------%
cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-3 3]; cfg.ga=ga;
cfg.contrast={'RemHi' 'ForHi'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')
%-------------------------------------------------------------------------%
% Interaxns
%-------------------------------------------------------------------------%
cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'Memory_x_Left_pos' 'Memory_x_Left_neg'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'invRT_x_Left_pos' 'invRT_x_Left_neg'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'Living_x_Left_pos' 'Living_x_Left_neg'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')

cfg=[]; cfg.template=template.topo; 
cfg.Inc=0.1; cfg.t=[-.1 1]; cfg.clim = [-1 1]; cfg.ga=ga;
cfg.contrast={'invRT_x_Memory_pos' 'invRT_x_Memory_neg'}; cfg.L=L; 
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet')
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet')
%-------------------------------------------------------------------------%
% ROIs
%-------------------------------------------------------------------------%
ROI1={'55','57','53','25','63','47','45'};
ROI2={'55' '53'};
ROI3={'53'};
ROI4={'36' '45' '37' '46'}; % Occip
ROI5={'1' '3' '5' '6' '9' '10'}; %midfront
ROI6={'38','49','47','57'}; % Lpar
ROI7={'Cz' '3' '14' '50' '38' '49' '13'}; % Cz cluster
ROI8={'23' '25' '27'};
ROI9={'37' '38' '48' '46' '36' '43' '45'}; % Pz cluster

% N2pc
tstruct=[]; tstruct.window=[-.2 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.25 .35]; tstruct.fig=[1 1 1 1]; tstruct.delta=1;
O1=erp_ss_plot_v3({'Memory_x_Left_pos' 'Memory_x_Left_neg'},get_el({'56' '54'},LABEL),tstruct);

tstruct=[]; tstruct.window=[-.2 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.25 .35]; tstruct.fig=[1 1 1 1]; tstruct.delta=0;
O1=erp_ss_plot_v3({'Right' 'mean' 'RemHi' 'Memory_x_Left_pos' 'invRT_x_Left_pos' 'Fast'},get_el(ROI2,LABEL),tstruct);

tstruct=[]; tstruct.window=[-.2 1]; tstruct.subj=L; tstruct.ga=ga;
tstruct.time=[.2 .4]; tstruct.fig=[1 1 1 1]; tstruct.delta=1;
O1=erp_ss_plot_v3({'Living' 'NonLiving'},get_el(ROI7,LABEL),tstruct);

%=========================================================================%
% FRQ
%=========================================================================%
addpath(fullfile(cod_dir,'EEG_toolboxes\fieldtrip-20170430\')); ft_defaults;
addpath(fullfile(cod_dir,'EEG_toolboxes\eeglab14_1_0b\'));

% Top Topo REM-FOR
cfg=[]; cfg.Inc=0.25; cfg.t=[-.5 .75]; cfg.foi=[4 7]; cfg.clim = [-.025 .025];
cfg.type='pow';  cfg.contrast={'slow' 'fast'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'top'}; jp_topoplotFT_v3(cfg); colormap('jet');

cfg=[]; cfg.Inc=0.2; cfg.t=[-.4 1]; cfg.foi=[8 14]; cfg.clim = [-.025 .025];
cfg.type='pow';  cfg.contrast={'slow' 'fast'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet');


% Back Topo REM-FOR
cfg=[]; cfg.Inc=0.2; cfg.t=[-.2 1]; cfg.foi=[8 14]; cfg.clim = [-.1 .1];
cfg.type='pow';  cfg.contrast={'fast_CI'}; 
cfg.L=L; cfg.ga=ga; cfg.template=template.topo;
cfg.perspective = {'back'}; jp_topoplotFT_v3(cfg); colormap('jet');


% Midfrontal Theta
cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el(ROI5,LABEL); 
cfg.time=[-.5 1]; cfg.foi=[3.5 7]; cfg.test=[0 .4];
cfg.contrast={'slow' 'fast'};
cfg.type='logpower'; tstruct.legend={'Remember' 'Forgotten'};
cfg.tstruct=tstruct; O1=frq_ss_plot_v2(cfg);

% Alpha CI
cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el(ROI2,LABEL); 
cfg.time=[-.5 1]; cfg.foi=[8 14]; cfg.test=[.5 1];
cfg.contrast={'slow_CI' 'fast_CI'};
cfg.type='logpower'; tstruct.legend={'Fast' 'Slow'};
cfg.tstruct=tstruct; O1=frq_ss_plot_v2(cfg);


% Alpha CI
cfg=[]; cfg.L=L; cfg.ga=ga;
cfg.elec=get_el({'36' '37'},LABEL); 
cfg.time=[-.6 1]; cfg.foi=[8 14]; cfg.test=[-.4 -.2];
cfg.contrast={'slow' 'fast'};
cfg.type='logpower'; tstruct.legend={'Fast' 'Slow'};
cfg.tstruct=tstruct; O1=frq_ss_plot_v2(cfg);




cfg=[]; cfg.template=template.erp;
cfg.topo_lim=[-4 4];
ft_plot_freq('slow','fast',ga,cfg); colormap('jet');