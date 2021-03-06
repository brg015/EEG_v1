function eeg_sea_enc_setup()
global RUN;
%=========================================================================%
% 3446
%=========================================================================%
% Encoding data looked excellent for this subject, all impedences were
% great and overall the data looked excellent
s=1;
subj{s}.name='3446';
subj{s}.pre.interp=[51];
subj{s}.pre.ICA=[1];
% Electrode 52 has HF noise (index is 52-1=51)
% Chopped filtered data looks fine
% Reject 467 for lat eye movements
s=s+1;
subj{s}.name='3447';
subj{s}.pre.interp=[28];
subj{s}.pre.ICA=[11];
% 10 looks odd on the ICAs, spatially it looks like blinks, but is much to
% slow
% Reject 16 for lat eye movements

s=s+1;
subj{s}.name='3448';
subj{s}.pre.interp=[47];
subj{s}.pre.ICA=[1 2 3 4 5];
% Not sure how these got all spread out, the data looks nosier than I'd
% like
% Reject 58 for lat eye movements

s=s+1;
subj{s}.name='3449';
subj{s}.pre.interp=[8,47];
subj{s}.pre.ICA=[2 3];
% ICA 2 = heartbeat
% ICA 3 = blinks (might have lateral eyemovent as well)
% Reject 18 for lat eye movements

s=s+1;
subj{s}.name='3450';
subj{s}.pre.interp=[52];
subj{s}.pre.ICA=[1];
% Reject 11 for lat eye movements

s=s+1;
subj{s}.name='3451';
subj{s}.pre.interp=[47];
subj{s}.pre.ICA=[5];
% 48 is broken
% Pretty sure said subject has shitty data - RTs are absurd

s=s+1;
subj{s}.name='3452';
subj{s}.pre.interp=[47 63];
subj{s}.pre.ICA=[1 2];
% 48 is broken, 64 has HF noise (muscle?)

s=s+1;
subj{s}.name='3453';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[2];
% Might not want to reference mastoid...

s=s+1;
subj{s}.name='3454';
subj{s}.pre.interp=[23 49];
subj{s}.pre.ICA=[9];
% 23 and 50 are noisy interp them

s=s+1;
subj{s}.name='3455';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[2];
% all channels look okay
% 2 is blinks, there are many

s=s+1;
subj{s}.name='3456';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1 4];
% all channels look okay
% 4 is HR

s=s+1;
subj{s}.name='3457';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[5];
% 17 and 19 have HF, but leaving for now.

s=s+1;
subj{s}.name='3458';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1 2];
% manually added trial drop code at 117827

s=s+1;
subj{s}.name='3480';
subj{s}.pre.interp=[51];
subj{s}.pre.ICA=[1 2];
% channel 52 is choppy

s=s+1;
subj{s}.name='3483';
subj{s}.pre.interp=[51];
subj{s}.pre.ICA=[1 2];

s=s+1;
subj{s}.name='3484';
subj{s}.pre.interp=[51];
subj{s}.pre.ICA=[];
% Remove no components

s=s+1;
subj{s}.name='3486';
subj{s}.pre.interp=[8,14];
subj{s}.pre.ICA=[1,4];

s=s+1;
subj{s}.name='3487';
subj{s}.pre.interp=[33, 40];
subj{s}.pre.ICA=[1,2,3,4];
% interp 34 and 41

s=s+1;
subj{s}.name='3488';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3490';
subj{s}.pre.interp=[28];
subj{s}.pre.ICA=[1];
% interp none

s=s+1;
subj{s}.name='3491';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[];
% 104/640 rejected trials seemms high

s=s+1;
subj{s}.name='3492';
subj{s}.pre.interp=[17];
subj{s}.pre.ICA=[1 2 3];

s=s+1;
subj{s}.name='3500';
subj{s}.pre.interp=[28 29 36 51];
subj{s}.pre.ICA=[4];
% lol 366/640 trials rejected in 'ica' condition
% 31/640 marked for rejection in 'pre'

s=s+1;
subj{s}.name='3501';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
% 1/640 marked in 'ica'
% 9/640 marked in 'pre'

s=s+1;
subj{s}.name='3507';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3509';
subj{s}.pre.interp=[36];
subj{s}.pre.ICA=[1];


s=s+1;
subj{s}.name='3510';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3511';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1,2,7];

s=s+1;
subj{s}.name='3514';
subj{s}.pre.interp=[38];
subj{s}.pre.ICA=[1,6];

% Interp LM
s=s+1;
subj{s}.name='3517';
subj{s}.pre.interp=[28,33];
subj{s}.pre.ICA=[2];

% Interp LM
s=s+1;
subj{s}.name='3519';
subj{s}.pre.interp=[6,33];
subj{s}.pre.ICA=[1,2];

s=s+1;
subj{s}.name='3520';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1,2,4,5,6,7,8,9,10,11,12];

s=s+1;
subj{s}.name='3525';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[];

s=s+1;
subj{s}.name='3526';
subj{s}.pre.interp=[47,27];
subj{s}.pre.ICA=[];
RUN.dir.plot(s)=0;

s=s+1;
subj{s}.name='3528';
subj{s}.pre.interp=[28];
subj{s}.pre.ICA=[];

s=s+1;
subj{s}.name='3540';
subj{s}.pre.interp=[3];
subj{s}.pre.ICA=[];

s=s+1;
subj{s}.name='3545';
subj{s}.pre.interp=[3];
subj{s}.pre.ICA=[];

s=s+1;
subj{s}.name='3547';
subj{s}.pre.interp=[47,27];
subj{s}.pre.ICA=[];

for ii=1:length(RUN.dir.subjects)
    if ~strcmp(RUN.dir.subjects{ii},subj{ii}.name)
        error('Subject IDs defined wrong');
    end
end
%=========================================================================%
%% Add behave data
%=========================================================================%
for ii=1:length(RUN.dir.subjects)
    dat=sea_enc_beh(ii);
    subj{ii}.beh=dat;
end

RUN.subj=subj;
%=========================================================================%
%% Set Processing Defualts
%=========================================================================%
switch RUN.dir.ver
    case 'v0'
        RUN.pre.ref=[32 64];
        
        RUN.pre.epoch=[-1 1.8];
        RUN.pre.baseline=[-100 0];

        RUN.set.filter_ica=[1 70];
        RUN.set.filter_dta=[0.05 70];
        RUN.set.sf=250;
        RUN.set.thresh_ica=[-500 1500];
        
        RUN.set.thresh_epc=[-.5 1];
        RUN.set.thresh_dta=[-100 100]; 

        RUN.set.reref='mastoid';
end
%=========================================================================%
RUN.template.ced='mw64_withMastoids_fixed.ced';
RUN.template.ica='mw64.ced';
RUN.template.elp='/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp';
RUN.template.plt='F:\SEA\Enc_EEG\EEG_v1-master\templates\layoutmw64_ActiChamp.mat';
RUN.template.plt='F:\SEA\Enc_EEG\EEG_v1-master\templates\layoutmw64.mat';
RUN.template.plt2='F:\SEA\Enc_EEG\EEG_v1-master\templates\layoutmw64_martyPlot.mat';
RUN.template.plt3='C:\Users\brg13\Desktop\SEA\layoutmw64.mat';

RUN.template.COCA='D:\Data\SEFER\EEG\template\COCA.csv';
RUN.template.dprime='D:\Data\SEFER\EEG\template\dprime.csv';

        















