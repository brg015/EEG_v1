function eeg_sefer_setup_student()
global RUN;
% Personal Notes:
% Interperlation is done at the very end. Thus, interpolerated channels
% greater than 31 should be referenced as the channel name minus 1
%=========================================================================%
% 2771
%=========================================================================%
s=1;
subj{s}.name='2771';
%no bad channels, 31 looked a little bad
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[1 2 3];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%

% 2781
%=========================================================================%
subj{s}.name='2781';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[1 2];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1 2];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2906
%=========================================================================%
subj{s}.name='2906';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[41];
subj{s}.enc.pre.ICA=[1 2];
%Channel 42 very noisy and all over the place //redo 2906 pre, did not
%remove IC's 1, 2
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[33 39 51];
subj{s}.ret.pre.ICA=[1 2 3]; %seems to all be eyeblinks
s=s+1;
%=========================================================================%
%=========================================================================%
% 2928
%=========================================================================%
subj{s}.name='2928';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[1 2];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1 4];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2960
%=========================================================================%
subj{s}.name='2960';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[32 51 35];
subj{s}.enc.pre.ICA=[1 9];
%ch36 eye movement on back of head?
%ic 9 heart beat removed
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[38 39 50 51];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 3013
%=========================================================================%
subj{s}.name='3013';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[5];
subj{s}.ret.pre.ICA=[6];
s=s+1;
% Amazing looking data!
%=========================================================================%
%=========================================================================%
% 3010
%=========================================================================%
subj{s}.name='3010';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[2];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2997
%=========================================================================%
subj{s}.name='2997';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2991
%=========================================================================%
subj{s}.name='2991';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2983
%=========================================================================%
subj{s}.name='2983';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2926
%=========================================================================%
subj{s}.name='2926';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[45];
subj{s}.ret.pre.ICA=[2];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2918
%=========================================================================%
subj{s}.name='2918';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1 5];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2912
%=========================================================================%
subj{s}.name='2912';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[33 41 60];
subj{s}.ret.pre.ICA=[1 4];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2870
%=========================================================================%
subj{s}.name='2870';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1]; %eyeblink
s=s+1;
%=========================================================================%
%=========================================================================%
% 2862
%=========================================================================%
subj{s}.name='2862';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2836
%=========================================================================%
subj{s}.name='2836';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2793
%=========================================================================%
subj{s}.name='2793';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1 2];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2789
%=========================================================================%
subj{s}.name='2789';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[];
subj{s}.ret.pre.ICA=[1 3];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2779
%=========================================================================%
subj{s}.name='2779';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[34];
subj{s}.ret.pre.ICA=[];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2772
%=========================================================================%
subj{s}.name='2772';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[29 32];
subj{s}.ret.pre.ICA=[1];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2764
%=========================================================================%
subj{s}.name='2764';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[8];
subj{s}.ret.pre.ICA=[3];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2735
%=========================================================================%
subj{s}.name='2735';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[21 34 50];
subj{s}.ret.pre.ICA=[];
s=s+1;
%=========================================================================%
%=========================================================================%
% 2954
%=========================================================================%
subj{s}.name='2954';
subj{s}.enc.pre.interp=[];
subj{s}.enc.pre.interp_ica=[];
subj{s}.enc.pre.ICA=[];
subj{s}.ret.pre.interp=[];
subj{s}.ret.pre.interp_ica=[24 29];
subj{s}.ret.pre.ICA=[];
s=s+1;
%=========================================================================%

%% Set Processing Defualts
%=========================================================================%
RUN.dir.analyze=ones(1,length(RUN.dir.subjects));
RUN.set.ica=1;
RUN.err=0;
RUN.dir.check='natural';

for ii=1:length(subj), subject_list{ii}=subj{ii}.name; end

c=1;
for ii=1:length(RUN.dir.subjects)
    I=find(strcmp(RUN.dir.subjects{ii},subject_list));
    if ~isempty(I)
        
        RUN.enc.pre.interp{ii}=subj{I}.enc.pre.interp;
        RUN.enc.pre.interp_ica{ii}=subj{I}.enc.pre.interp_ica;
        RUN.enc.pre.ICA{ii}=subj{I}.enc.pre.ICA;
        
        RUN.ret.pre.interp{ii}=subj{I}.ret.pre.interp;
        RUN.ret.pre.interp_ica{ii}=subj{I}.ret.pre.interp_ica;
        RUN.ret.pre.ICA{ii}=subj{I}.ret.pre.ICA;
        
        if RUN.dir.analyze(ii)==1
            RUN.dir.inc_subjects{c}=RUN.dir.subjects{ii}; 
            c=c+1;
        end
        switch RUN.dir.ver
            case 'v0_20140708' 
                RUN.enc.pre.ref{ii}=setdiff([32,64],RUN.enc.pre.interp{ii});
                RUN.ret.pre.ref{ii}=setdiff([32,64],RUN.ret.pre.interp{ii});
            case 'v1_20140807'
                RUN.enc.pre.ref{ii}=[32 64];
                RUN.ret.pre.ref{ii}=[32 64];
        end
    else
        display([RUN.dir.subjects{ii} ' DNE']);
        error('Kindly select valid subjects');
    end
end
%=========================================================================%
RUN.template.ced='mw64_withMastoids_fixed.ced';
RUN.template.ica='mw64.ced';
RUN.template.elp='/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp';
RUN.template.plt='D:\Data\SEFER\EEG\template\layoutmw64_ActiChamp.mat';
RUN.template.plt='D:\Data\SEFER\EEG\template\layoutmw64.mat';
RUN.template.plt2='D:\Data\SEFER\EEG\template\layoutmw64_martyPlot.mat';
RUN.template.COCA='D:\Data\SEFER\EEG\template\COCA.csv';
RUN.template.dprime='D:\Data\SEFER\EEG\template\dprime.csv';
        
switch RUN.dir.ver
    case 'v2'
        % Lengthened these on 10/29/14
        RUN.enc.pre.epoch=[-1 2.5];
        RUN.enc.pre.baseline=[-400 -200];
        RUN.ret.pre.epoch=[-1 2.5];
        RUN.ret.pre.baseline=[-200 0];

        RUN.set.filter_ica=[1 70];
        RUN.set.filter_dta=[0.05 70];

        RUN.set.sf=250;
        RUN.set.thresh_ica=[-500 1500];
        % Changed from -150 as it was too sensitive
        
        RUN.LS='';
        switch RUN.LS
            case ''
                RUN.set.thresh_epc=[-.4 1.1];
                RUN.set.thresh_dta=[-120 120]; LS='';
            case '_strict'
                RUN.set.thresh_epc=[-.5 1.5];
                RUN.set.thresh_dta=[-100 100]; LS='_strict';
            case '_trend'
                RUN.set.thresh_epc=[-.4 1.5];
                RUN.set.thresh_dta=[-120 120]; LS='_trend';
            case '_longest'
                RUN.set.thresh_epc=[-.5 2];
                RUN.set.thresh_dta=[-100 100]; LS='_longest';
        end
        RUN.set.reref='mastoid';
        
        switch RUN.dir.check
            case 'natural', RUN.set.ch_interp=1; % RUN.set.ica=1;
            case 'channel', RUN.set.ch_interp=0; % RUN.set.ica=0;
        end
    otherwise
end

% Let's setup some save strings for posterity:
switch RUN.dir.hc
    case 0, hc_str='';
    case 1, hc_str='_hc';
    case 2, hc_str='_hc_only';
    case 3, hc_str='_hc_toss';
end
RUN.dir.beh=fullfile('D:\Data\','SEFER','behav','eeg',['pre' hc_str]);
        
sess_str=['_' RUN.dir.sess];

ref_str=['_' RUN.set.reref];
switch RUN.set.ica
    case 0, ica_str='';
    case 1, ica_str='_ica';
end

% ICrun files
RUN.save_str{1}=['ICrun' sess_str '.set'];
RUN.save_str{2}=['eegbasic' sess_str '_f' num2str(RUN.set.filter_ica(1)*100) '.set'];
% Pre files
RUN.save_str{3}=['eegbasic' sess_str '_f' num2str(RUN.set.filter_dta(1)*100) '.set'];
% Post-pre-long

epoch_length='long';
switch epoch_length
    case 'short'
        RUN.save_str{4}=['eegepoch' sess_str ica_str '_f' num2str(RUN.set.filter_dta(1)*100) '.set'];
        RUN.save_str{5}=['eegepoch' sess_str ica_str ref_str '_f' num2str(RUN.set.filter_dta(1)*100) hc_str '.mat'];
    case 'long'
        RUN.save_str{4}=['eegepoch_long' LS sess_str ica_str '_f' num2str(RUN.set.filter_dta(1)*100) '.set'];
        RUN.save_str{5}=['eegepoch_long' LS sess_str ica_str ref_str '_f' num2str(RUN.set.filter_dta(1)*100) hc_str '.mat'];
end















