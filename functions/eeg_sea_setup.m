function eeg_sea_setup()
global RUN;
%=========================================================================%
% 3446
%=========================================================================%
s=1;
subj{s}.name='2735';
subj{s}.pre.interp=[];
subj{s}.pre.interp_ica=[];
subj{s}.pre.ICA=[];

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
        
        % Limit
        % Reg
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















