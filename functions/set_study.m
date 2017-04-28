%% set the working directory
clear all
restoredefaultpath
cd \\ccn-woldorffser.win.duke.edu\data1\Projects\April\cuedAttIntf;
dataPath = 'C:\Users\adr25\Desktop\data\';
eeglabDir = 'matlab_apps\eeglab13_2_2b\';


subjectID = {'2700' '2711' '2713' '2717' '2718' '2729'}; 

%establish bad channels to be interpolated in preprocessing.m

interp = {[43 44 52] [8 15] [7 15 16 43 44] [16 52] [63 59 26 41] [18 28]}; 

%establish ICA components (eye blinks) to be subtracted
ICA = {[3] [1] [4] [1] [1] [11]};





%%

open analysis\eegAnalysis\preprocessing.m