%=========================================================================%
%% Initialize
%=========================================================================%
% restoredefaultpath;matlabrc;
clc; close all; clear all;
global RUN;

serv='local';
switch serv
    case 'local', root_dir='C:\Users\brg13\Desktop\MatlabCode\';
    case 'serv1', root_dir='D:\Data\';
end

wrk_dir=fullfile('F:\SEA\Enc_EEG\');

% Specify where your function files are
code_dir{1}=fullfile(root_dir,'function_files-master');
code_dir{2}=fullfile(root_dir,'EEG_v1-master\functions\');
code_dir{3}=fullfile(root_dir,'EEG_v1-master\ft_functions\');
code_dir{4}=fullfile(root_dir,'EEG_v1-master\templates\');
for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end
addpath(fullfile(root_dir,'EEG_toolboxes\fieldtrip-20170430\'));
addpath('C:\Users\brg13\Desktop\SEA\SEA_scripts\');
rmpath(fullfile(root_dir,'EEG_toolboxes\eeglab14_1_0b\'));
%=========================================================================%
%% Initialize
%=========================================================================%
% 3446 too many lat eye movements
%
RUN.dir.plot=true(1,38);
RUN.dir.subjects={'3446' '3447','3448','3449','3450',...
    '3451','3452','3453','3454','3455',...
    '3456','3457','3458','3480','3483',...
    '3484','3486','3487','3488','3490',...
    '3491','3492','3500','3501','3507',...
    '3509' '3510','3511','3514','3517',...
    '3519','3520','3525','3526','3528',...
    '3540','3545','3547'}; 
RUN.dir.plot(32)=0;
RUN.dir.plot(25)=0;
RUN.dir.plot(23)=0;
RUN.dir.plot(10)=0;
RUN.dir.plot(9)=0;
RUN.dir.plot(6)=0;
RUN.dir.plot(3)=0;
% RUN.dir.plot(1:32)=0;
RUN.dir.plot([34 37 38])=0;
RUN.dir.ver='v0';
RUN.dir.pro=fullfile(wrk_dir,'pro',RUN.dir.ver);
RUN.dir.sav=fullfile(wrk_dir,'group_freq',RUN.dir.ver);

RUN.dir.plot(1:38)=0;
RUN.dir.plot([35 36])=1;

eeg_sea_enc_setup();

eeg_plot_ERP_RT();
eeg_plot_FREQ_RT();
