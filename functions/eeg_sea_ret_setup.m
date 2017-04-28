function eeg_sea_ret_setup()
global RUN;
% If channel is >32 must subtract by one as channel 32 does not exist
% (reference)
%=========================================================================%
% 3446
%=========================================================================%
s=1;
subj{s}.name='3446';
subj{s}.pre.interp=[51];
subj{s}.pre.ICA=[2];
s=s+1;
% channel 52 has a great deal of HF noise and responds strongly
% (positively) to eye blinks

subj{s}.name='3447';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[];
s=s+1;

subj{s}.name='3448';
subj{s}.pre.interp=[47];
subj{s}.pre.ICA=[1 2 3 4];
s=s+1;
% removed channel 48 which is index 47
% UPDATE: chopped out the weird spikes along with the nearest '5' event
% code
% chopped out sections' latencies in seconds (event codes, if applicable,
% are indicated in parentheses):
% 811.29 - 811.58; 814.45 - 814.69; 820.86 - 821.41 (4); 827.50; 830.70;
% 832.97; 877.78 - 878.99 (5); 884.48; 944.83 - 945.00; 947.62;
% 1022.32 - 1022.78; 1033.19 - 1033.74; 1044.83 - 1047.47 (2,5);
% 1122.54; 1125.00; 1127.34 - 1127.46; 1128.81; 1212.15; 2003.98 - 2004.79 (4)
% 2047.32 -2047.44; 2049.08 - 2049.13; 2056.36 - 2057.28 (4); 
% 2058.48 - 2058.70; 2065.01 - 2066.45 (5);

subj{s}.name='3449';
subj{s}.pre.interp=[8 47];
subj{s}.pre.ICA=[1 2 3];
s=s+1;
% removed channels 8 and 48 (index 47)
% considerable alpha

subj{s}.name='3450';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1 6];
s=s+1;
% All channels look good. Minimal weirdness

subj{s}.name='3451';
subj{s}.pre.interp=[47];
subj{s}.pre.ICA=[1 4];
s=s+1;
% some strange drifts lasting 2-4 epochs

subj{s}.name='3452';
subj{s}.pre.interp=[47];
subj{s}.pre.ICA=[1 2];
s=s+1;
% tons of blinks (just less than once every epoch)

subj{s}.name='3453';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[3 7 8];
s=s+1;
% All channels look good

subj{s}.name='3454';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[6 8];
s=s+1;

subj{s}.name='3455';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3456';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3457';
subj{s}.pre.interp=[17 19];
subj{s}.pre.ICA=[3];
s=s+1;
% Tons of noise at times, especially following breaks. Looks like muscular
% activity. Smooths out after 10-20 seconds. 17 and 19 are particulary
% noisy - sharp spikes at ~11 Hz.

subj{s}.name='3458';
subj{s}.pre.interp=[17];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3480';
subj{s}.pre.interp=[26];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3483';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1 4];
s=s+1;
% Lots of intermittent HF noise. A bunch of electrodes are at times noisy
% then normal. Quite a bit of drift too.

subj{s}.name='3484';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[];
s=s+1;
% Loads of weird identical drift on greatly disparate channels (on 
% different bundles)
% ICA can't find a blink component

subj{s}.name='3486';
subj{s}.pre.interp=[14];
subj{s}.pre.ICA=[3];
s=s+1;
% Troublesome trial. Channel 14 acted up throughout, and all channels
% flatlined towards the end.
% ICA plots look very funny - no clear eyeblink component

subj{s}.name='3487';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3488';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3490';
subj{s}.pre.interp=[6];
subj{s}.pre.ICA=[1];
s=s+1;

subj{s}.name='3491';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[2];
s=s+1;
% muscle activity on channels 17 and 19 comes and goes

subj{s}.name='3492';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3500';
subj{s}.pre.interp=[28 29 36 51];
subj{s}.pre.ICA=[2];
% 111/768 trials marked for rejection in 'pre'

s=s+1;
subj{s}.name='3501';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3507';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[8];

s=s+1;
subj{s}.name='3509';
subj{s}.pre.interp=[36];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3510';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1];
% running ICA scenario for this yields an error: "Cell contents assignment
% to a non-cell array object"

s=s+1;
subj{s}.name='3511';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[1,2,8];
% 8 is heartrate
% a great deal of identical drift is present across multiple channels. Due
% to the bad electrode which has now been fixed

s=s+1;
subj{s}.name='3514';
subj{s}.pre.interp=[38];
subj{s}.pre.ICA=[2];
% Component 1 is bad channel at start

s=s+1;
subj{s}.name='3517';
subj{s}.pre.interp=[];
subj{s}.pre.ICA=[3];
% Component 1 is reference noise - its bad

s=s+1;
subj{s}.name='3519';
subj{s}.pre.interp=[61];
subj{s}.pre.ICA=[1];

s=s+1;
subj{s}.name='3520';
subj{s}.pre.interp=[6,28];
subj{s}.pre.ICA=[1,2,4];
% Tourette's - gonna be booted anyways

for ii=1:length(RUN.dir.subjects)
    if ~strcmp(RUN.dir.subjects{ii},subj{ii}.name)
        error('Subject IDs defined wrong');
    end
end
%=========================================================================%
%% Add behave data
%=========================================================================%
for ii=1:length(RUN.dir.subjects)
    dat=sea_ret_beh(ii);
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
RUN.template.plt='D:\Data\SEFER\EEG\template\layoutmw64_ActiChamp.mat';
RUN.template.plt='D:\Data\SEFER\EEG\template\layoutmw64.mat';
RUN.template.plt2='D:\Data\SEFER\EEG\template\layoutmw64_martyPlot.mat';
RUN.template.plt3='F:\Data2\SEA\layoutmw64.mat';

RUN.template.COCA='D:\Data\SEFER\EEG\template\COCA.csv';
RUN.template.dprime='D:\Data\SEFER\EEG\template\dprime.csv';

        















